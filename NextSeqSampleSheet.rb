#!/usr/bin/ruby
# filename: NextSeqSampleSheet.rb
# creates sample sheet formatted for bcl2fastq
# executed in cronjob_Illumina_Novaseq.sh

def main()
  input = ARGV[0] # origial sample sheet
  output = ARGV[1] # formated sample sheet

  infile = File.new(input, "r")
  outfile = File.new(output, "w")
  flag = 1 # counter denoting first line
  infile.each {
    |line|

    cols = line.chomp.split(/,/)

    if flag == 2 # start writing

      outfile.write "[Header],," + "\n"
      outfile.write "IEMFileVersion,5,," +"\n"
      temp = `date`
      outfile.write "Date," + temp + "\n"
      outfile.write "Workflow,GenerateFASTQ" + "\n" 
      outfile.write "Application,NextSeq FASTQ Only" +"\n"
      outfile.write "Instrument Type,NextSeq,,,,,,,,," +"\n"
      outfile.write "Assay,TruSeq Stranded Total RNA,,,,,,,,," +"\n"
      outfile.write "Index Adapters,IDT-ILMN TruSeq RNA UD Indexes (96 Indexes),,,,,,,,," +"\n"
      outfile.write ",,,,,,,,,," +"\n"


      outfile.write "Description,," + "\n" + "Chemistry,Amplicon,," + "\n" + ",," + "\n\n" + "[Reads],," + "\n"

      outfile.write "[Reads],," + "\n"
      sequenceLength = cols[9].split("_")
      arrlength=sequenceLength.length()
      num = sequenceLength[arrlength-2].split("E")
      outfile.write num[1].to_s + ",," +"\n"
  
      if num[0]== "P"
        outfile.write num[1] + ",," +"\n"
      end
 
      outfile.write "\n"

      outfile.write "[Settings],," + "\n"
      # adapters for both Truseq and Nextera are included here for trimming
      outfile.write "Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA+CTGTCTCTTATACACATCT,,,,,,,,," +"\n"
      outfile.write "AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT+CTGTCTCTTATACACATCT,,,,,,,,," +"\n"

      outfile.write ",,\n[Data]" +"\n"
      outfile.write "Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description" +"\n"
    end

# start the sample portion of samplesheet after creating header using first line (flag = 1)
    if (flag >= 2)

      sequenceLength = cols[9].split("_")
      arrlength=sequenceLength.length()                                                                                                                            
      num = sequenceLength[arrlength-2].split("E")      

# check if dual index      
      if cols[4].length > 10
        temp1 = cols[4].split("-")
      # NextSeq needs reverse comp of second index
        reverseTemp = reverseComp(temp1[1])
        outfile.write cols[1]+ "," + cols[2] + "," + cols[2] + "," + "," + "," + "," + temp1[0] + "," + "," + reverseTemp + "," + cols[9] + "\n"
      else
        outfile.write cols[1]+ "," + cols[2] +"," + cols[2] +"," + "," + "," + "," + cols[4]+ "," + "," + "," + cols[9] + "\n" 
      end
    end
    flag = flag + 1
}
  
  infile.close()
  outfile.close()
end


def reverseComp(index)
  temp = ""
  count = index.length-1
  
  for i in 0..index.length-1
    if index[i] == "T"
      temp = "A".concat(temp)
    elsif index[i] == "A"
      temp = "T".concat(temp)
    elsif index[i] == "C"
      temp = "G".concat(temp)
    elsif index[i] == "G"
      temp = "C".concat(temp)
    end
    count = count -1
  end  
  return temp   
end

main()
