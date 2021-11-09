import sys
from Bio import SeqIO

def replace_fasta_header(fasta,outfile):

	out_file=open(outfile,'w')

	for seq_rec in SeqIO.parse(fasta,"fasta"):
		header_ori=str(seq_rec.description).split()
		header_replace=''
		for i in header_ori:
			if 'locus_tag=' in i:
				header_replace=i.replace('[locus_tag=','').replace(']','')

		out_file.write('>'+str(header_replace)+'\n'+str(seq_rec.seq)+'\n')
	out_file.close()

replace_fasta_header(sys.argv[1],sys.argv[2])
