import sys
from Bio import SeqIO

#>gi|1143527503|dbj|AP017959.1|:c439870-439626 Synechococcus sp. NIES-970 DNA, complete genome-p.c.VAL:55.97%-taxID:1827144
#>AP017959.1:c439870-439626 Synechococcus sp. NIES-970 DNA, complete genome-p.c.VAL:55.97%-taxID:1827144

def fix_header_glassgo2copra(filein,fileout):#query file contains many more sequences than the sequences in target file, (target file is synteny table file, query_file is fasta file)
	#skip query sequence
	out_file=open(fileout,'w')
	for seq_rec in SeqIO.parse(filein,"fasta"):
		if 'gi|' in seq_rec.description:
			header=str(seq_rec.description)
			header=header.split('|')[3]+header.split('|')[4]

			out_file.write('>'+str(header)+'\n'+str(seq_rec.seq)+'\n')	

	out_file.close()
	#print(dict_target)

fix_header_glassgo2copra(sys.argv[1],sys.argv[2])
