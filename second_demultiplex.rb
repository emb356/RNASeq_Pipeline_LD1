# file name: second_demultiplex.rb
# performs plateseq demultiplexing
# splits fastq files based on barcode
# barcode matching using hamming distance (<=1 to real and >1 to the other)
# generates a set of possible codes with hd <=1 as a hash

require 'zlib'
require 'multiple_files_gzip_reader'

def main
  inputdir = ARGV[0]
  outdir = ARGV[1]
  prefix = ARGV[2]
  barcode = ARGV[3] # comma-delimited csv file
  
  nt = ARGV[4]  # max threads  
  if nt == nil
    nt = 4
  else
    nt = nt.to_i
  end

  barcodesize = 9  # barcode size

  coding = readBar(barcode) # return a hash {lane => { code => sample  }}
 
  multiplex = {}
  $reads_lane ={}

  coding.each do |lane, ch|
    multiplex[lane] = {}
    ch.each do |str,sampleID|
      mutate1(str).each do |strmut|   # only allow hamming distance of 1
        multiplex[lane][strmut] = sampleID + "_" + str
      end
    end
  end
  $stderr.puts "multiplex mapping: #{multiplex}"

  outprefix = outdir + "/" + prefix 

  assignment = decode(inputdir, multiplex, outprefix, barcodesize, nt)
  
  # print summary on reads - mean and stdev per lane
  $reads_lane.each do |lane,arr_samples| 
         $stderr.puts "Lane #{lane} \t" + arr_samples.join("\t") 
  end
  $stderr.puts "\n"
  $reads_lane.each do |lane,arr_samples|
         mean = arr_samples.inject{ |sum, el| sum + el }.to_f / arr_samples.size
         variance = arr_samples.inject(0) { |variance, x| variance += (x - mean) ** 2 }
         stddev = Math.sqrt(variance/(arr_samples.size-1))
         $stderr.puts "Lane #{lane}\tMean: #{mean.to_i}\tStdev: #{stddev.to_i}"
  end
end


def decode(inputdir, multiplex, outprefix, barcodesize, nt )
  
  barcodefq = Dir.new(inputdir).select {|a| a.match(/\w+d*w*\_\w+\d+\_\L00\d\_R1_\d+\.fastq.gz$/) }  
  $stderr.puts "barcode files: \n#{barcodefq.join("\n")}"
  nprocess = 0

  barcodefq.sort.each do |bfq|

    if bfq.match(/(\w+d*w*)\_\w+(\d+)\_\L00(\d)\_R1_(\d+)\.(\S+)/)
      sampleID = "#{$1}"
      sampleOrder= "#{$2}"
      lane = "#{$3}"
      piece = "#{$4}"
      $stderr.puts "sampleID:#{sampleID}"
      $stderr.puts "Piece:#{piece}"
      $stderr.puts "Lane:#{lane}"
      targetfq1 = "#{inputdir}/#{sampleID}\_S#{sampleOrder}_\L00#{lane}\_R2_#{piece}\.#{$5}"
      $stderr.puts "targetFile:#{targetfq1}"

      
      bfql = "#{inputdir}/#{bfq}"
      next unless multiplex.key?(lane)
      $stderr.puts "working on lane #{lane}"
      Process.fork do
        doSplit(sampleID, bfql, targetfq1, lane, piece, multiplex, outprefix, barcodesize) if File.exist?(targetfq1)
      end

      nprocess += 1

      if nprocess >= nt # need wait
        Process.waitall
        nprocess = 0 
      end
    end
  end
  Process.waitall
end


def doSplit(sID, bfq, targetfq, lane, piece, multiplex, outprefix, barcodesize)
  
  ndecode = 0
  ndiscard = 0


  ends = "1"
  if targetfq =~ /\w+d*w*\_\w+\d+\_\L00(\d)\_R2_(\d+)\.f(\S+)$/  
    ends = $2
  end


  outio = {}
  reads_sample ={}
  multiplex[lane].values.sort.each do |sampleID|
  
    outName = "#{sID}#{sampleID}_L00#{lane}_R1_#{piece}.fastq.gz"
  
    outio[sampleID] = Zlib::GzipWriter.new(File.open(outName,'w')) 
  
    reads_sample[sampleID] =0;
 end
  
  outio[:discard] = Zlib::GzipWriter.new(File.open(outprefix + "_L00#{lane}_#{ends}.unknown.fastq.gz", 'w'))    

  print "#{bfq}"
  print "#{targetfq}"

  if bfq.match(/.gz$/)
    biof = File.new(bfq, 'r')
    bio = MultipleFilesGzipReader.new(biof).to_enum
  else
    bio = File.new(bfq, 'r')
  end

  if targetfq.match(/.gz$/)
    tiof = File.new(targetfq, 'r')
    tio = MultipleFilesGzipReader.new(tiof).to_enum
  else
    tio = File.new(targetfq, 'r')
  end

  
  loop do
    bunit = []
    tunit = []
    bline = bio.next()
    bunit << bline
    tline = tio.next()
    tunit << tline
    bline = bio.next()
    bunit << bline
    tline = tio.next()
    tunit << tline
    bline = bio.next()
    bunit << bline
    tline = tio.next()
    tunit << tline
    bline = bio.next()
    bunit << bline
    tline = tio.next()
    tunit << tline

    bc = bunit[1][0..barcodesize-1]
    
    if multiplex[lane].key?(bc)  
      
      sampleID = multiplex[lane][bc]
      
      outio[sampleID].puts tunit
      reads_sample[sampleID] += 1
      ndecode += 1
    else # discarded
      outio[:discard].puts tunit
      ndiscard += 1
    end
  end

 # Output summary of reads per sample and reads per lane etc.

  $stderr.puts "#{targetfq}:\t#decoded=#{ndecode}\t#unknown=#{ndiscard}"
  read_summary_file =  File.open(outprefix + "_summary_lane#{lane}.stats",'a+')
  lane_out="lane"
  temp_out="#{lane}"
  $reads_lane[lane] = []

  reads_sample.each do |sample,reads|
	lane_out << "\t#{sample}"
	temp_out << "\t#{reads}"
	$reads_lane[lane].push(reads)	
  end
  lane_out << "\tDecoded"
  temp_out << "\t#{ndecode}"
  lane_out << "\tUnknown"
  temp_out << "\t#{ndiscard}"
  read_summary_file.puts lane_out
  read_summary_file.puts temp_out
  read_summary_file.close

  outio.values.each {|oio| oio.close}
  bio.close
  tio.close

end

# create hash for hamming distance of 1
def mutate1(string)
  l = string.length
  hd1 = []
  0.upto(l-1) do |i|
    cp = String.new(string)
    ["A", "T", "G", "C", "."].each do |nt|
      cp[i] = nt
      hd1 << String.new(cp)
    end
  end
  return hd1
end


def readBar(b)
  # parses the barcode csv file (modified from ~/Pipeline_LD1/References)
  coding = {}
  File.new(b, 'r').each do |line|
    next if line.match(/^#/) # header line
    
    cols = line.chomp.split(/,/)  
    lane, sampleID, code = cols[0].strip, cols[1].strip, cols[2].strip

    if !coding.key?(lane)
      coding[lane] = {}
    end

	if coding[lane].key?(code)
                $stderr.puts "Duplicate barcodes in lane : #{lane}  for #{code} from sample #{coding[lane][code]}"
		exit
	end

    coding[lane][code] = sampleID.tr("/", "_").tr(" ","_")
  end
  return coding
end


def sanity_check(inputdir, outprefix, coding)
 flag=1
   coding.each do |lane,junk|	
		num_samples_this_lane = coding[lane].length
		count = 0
		all_codes=""
		bfq="s_#{lane}_2.fastq.barcode-stats"
		f = File.new("#{inputdir}/#{bfq}","r")
		i = 0
		begin
		    while (line = f.readline)
		        line.chomp
			break if i == 10
			cols = line.chomp.split(/\t/)
			all_codes <<  cols[0].strip 
			all_codes << "\t"
			count += 1
			i += 1
		    end
		rescue EOFError
		    f.close
		end
			sanity_flag=1
			coding[lane].each_key do |code|
		if all_codes.scan( "#{code}A" ).empty? and all_codes.scan( "#{code}C" ).empty? and all_codes.scan( "#{code}G" ).empty? and all_codes.scan( "#{code}T" ).empty?
					sanity_flag=0
					break
				else
				end
			end
			if sanity_flag == 1
				$stderr.puts "Sanity Check - Lane #{lane} : Successful"
			else
                        	$stderr.puts "Sanity Check - Lane #{lane} : Failed"
	                        flag=0
                	end
 end
 return flag
end

main()


