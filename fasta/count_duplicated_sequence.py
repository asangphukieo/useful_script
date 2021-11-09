
#build_tree dependency
from Bio import SeqIO
import sys

seq_dict={}

for seq_rec in SeqIO.parse(sys.argv[1],"fasta"):
	#print(seq_rec.id)
	pos=seq_rec.id.split(':')
	acc=pos[0]
	#print(len(seq_rec))

	if 'c' in pos[1]: #reverse transcipt
		col=pos[1].replace('c','')
		col=col.split('-')
		str_srna=col[1]
		stop_srna=col[0]
		#print(pos[1],str_srna,stop_srna)

	else:
		col=pos[1].split('-')
		str_srna=col[0]
		stop_srna=col[1]
		#print(pos[1],str_srna,stop_srna)

	if acc not in seq_dict:
		seq_dict[acc]=[[str_srna,stop_srna,seq_rec.description,seq_rec.seq,len(seq_rec)],]

		#print('New data add',acc,seq_dict[acc])
	
	else:
		#print('Redundant',acc,len(seq_dict[acc]))
		#print(seq_dict[acc])
		collect_i=[]
		#print(acc,seq_dict[acc])
		for i in seq_dict[acc]:
			#collect_i.append(i)

			collect_i.append(i) #add existing genes to db
		
		collect_i.append([str_srna,stop_srna,seq_rec.description,seq_rec.seq,len(seq_rec)])
		del seq_dict[acc]
		seq_dict[acc]=collect_i

		#print(seq_dict[acc])

		#print('Update Redundant',acc,len(seq_dict[acc]))
	#seq_dict[seq_rec.id]=seq_rec.seq
	#key.append(seq_rec.id)

w2stat=open('gene_duplicate_stat.temp','w')
for k in seq_dict:
	
	w2stat.write(seq_dict[k][0][2].split(' ')[1]+' '+seq_dict[k][0][2].split(' ')[2]+'\t'+k+'\t'+str(len(seq_dict[k]))+'\n')
w2stat.close()
print('Done!!')

#usage
#python count_duplicated_sequence.py ../test_dup.fasta


