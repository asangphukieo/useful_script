import sys
from Bio import SeqIO
#to select the most similar sequence from GLASSgo output

input_file=sys.argv[1] #GLASSgo output
outfile=sys.argv[2] 

out_file=open(outfile,'w')
i=1
for seq_rec in SeqIO.parse(input_file,"fasta"):
	max_ident=0.0
	closest_seq=''
	if i >= 2: #Remove first sequence as GLASSgo query
		identity=seq_rec.description.split('p.c.VAL:')[1].split('%-')[0]
		if float(identity) > max_ident:
			max_ident=float(identity)
			closest_seq='>'+str(seq_rec.id.replace('|','_'))+'\n'+str(seq_rec.seq)+'\n'

		
	i+=1
out_file.write(closest_seq)
out_file.close()

#python select_closest_seq.py ./../OUTPUT/lastest_14/01_GLASSgo_Results/GLASSgo_output_6.fasta ./test.out