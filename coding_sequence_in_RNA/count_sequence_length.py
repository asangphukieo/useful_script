from Bio import SeqIO
import sys

for seq_rec in SeqIO.parse(sys.argv[1],"fasta"):
	header=seq_rec.description.replace('   ','  ')
	#header=header.split('  ')[1].split('_')[0].replace(' ','')
	header=header.split('  ')[1]
	print(header,len(str(seq_rec.seq)))

#to count uniq genome
#python count_sequence_length.py renamed_SR1P_blastout.fasta |cut -f1 -d' '|sort -u|wc -l

