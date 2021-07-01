#!/usr/bin/ruby
# filename: SCC_NovaSeqSampleSheet.rb
# creates sample sheet formatted for cellRanger
# executed in cronjob_Illumina_Novaseq.sh

def main()
 input = ARGV[0] # original sample sheet
 sccout = ARGV[1] # SCC formatted sample sheet

  infile = File.new(input, "r")
  sccfile = File.new(sccout, "w") # SCC samplesheet for cellRanger

  # counter for headers
  flag = 1

  infile.each {
    |line|

    cols = line.chomp.split(/,/)
    if flag ==2
 
      sccfile.write "[Header]" + "\n"
      sccfile.write "IEMFileVersion,4" +"\n\n"
      sccfile.write "[Reads]" + "\n"
      sccfile.write "26" +"\n" + "8" + "\n\n"
      sccfile.write "[Data]" + "\n"
      sccfile.write "Lane,Sample_ID,index,Sample_Project" +"\n"
    end
 
    
    if (flag >= 2)

      sequenceLength = cols[9].split("_")
      arrlength=sequenceLength.length()                                                                                                                            
      num = sequenceLength[arrlength-2].split("E")      
      # identify single cell projects
      if cols[9].include? "_SCC_"
        sccfile.write cols[1]+ "," + cols[2] + "," + cols[4] + "," + cols[9] + "\n"
      end
    end
    flag = flag + 1
   
  }

  infile.close()
  sccfile.close()
end


main()
