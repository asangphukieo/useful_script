import sys
from Bio import SeqIO
#to select the most similar sequence from GLASSgo output

input_file1=sys.argv[1] #GLASSgo output 1 as first round
input_file2=sys.argv[2] #GLASSgo output 2 from repeting GLASSgo by using closest sequence from glassgo round 1
outfile=sys.argv[3] 

id_list={}
out_file=open(outfile,'w')
#>gi|939195038|gb|CP012832.1|:c3161297-3161230 Synechocystis sp. PCC 6803 substrain GT-G, complete genome-p.c.VAL:100.0%-taxID:1148
for seq_rec in SeqIO.parse(input_file1,"fasta"):
	id_list[str(seq_rec.id)]='found'
	out_file.write('>'+str(seq_rec.description)+'\n'+str(seq_rec.seq)+'\n')


i=1
for seq_rec in SeqIO.parse(input_file2,"fasta"):
	if i >= 2: #Remove first sequence as GLASSgo query
		if seq_rec.id not in id_list:
			out_file.write('>'+str(seq_rec.description)+'-(rep_1)\n'+str(seq_rec.seq)+'\n')
		
	i+=1

out_file.close()

#python combine_seq.py ./../OUTPUT/lastest_14/01_GLASSgo_Results/GLASSgo_output_6.fasta ./../OUTPUT/lastest_14/01_GLASSgo_Results/GLASSgo_output_6.fasta ./test.out