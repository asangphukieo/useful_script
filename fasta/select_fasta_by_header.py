import sys
from Bio import SeqIO
#require biopython

def select_fasta_from_name(name_list,fasta,outfile): #name_list,fasta,outfile
	list_fasta=[]
	out_file=open(outfile,'w')	
	for seq_rec in SeqIO.parse(fasta,"fasta"):
		list_fasta.append('>'+str(seq_rec.description)+'\n'+str(seq_rec.seq)+'\n')

	for code in open(name_list):
		code=code.rstrip()
		#print(code)
		code=code.split('\t') #split tab to select the seconf colund which is accession id
		
		for fasta_ind in list_fasta:

			if code[0] in fasta_ind: #Remove first sequence as GLASSgo query	
				out_file.write(fasta_ind)
				
	out_file.close()

select_fasta_from_name(sys.argv[1],sys.argv[2],sys.argv[3])

#usage:
#python select_fasta_by_header.py selected_seq.txt target.fasta output_selected.fasta
#input selected_seq.txt >> file containing header of selected sequencuences 
