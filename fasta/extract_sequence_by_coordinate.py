from Bio import SeqIO
import sys

#use conda of build_tree

inFile = open(sys.argv[1],'r') #input fasta file
fw=open(sys.argv[3],'w') #output fasta file

for record in SeqIO.parse(inFile,'fasta'):
	
	for line in open(sys.argv[2]): #input coordinate file form blastout

		line= line.rstrip()
		col=line.split('\t')
		seq_name=col[1]
		start=int(col[6])-1
		end=int(col[7])

		if str(record.description) == seq_name:
			
			fw.write(">" + record.description + "\n")
			fww = (str(record.seq[start:end]) + '\n')
			fw.write(fww)
    
fw.close()
inFile.close()
#usage
#python extract_sequence_by_coordinate.py rename_pass_glassgo2copra.fasta SR1P.blastout rename.out
