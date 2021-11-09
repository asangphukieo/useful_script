
#build_tree dependency
from Bio import SeqIO
import sys

#remove redundant sequences but keep paralogs 
#keep longer sequence 

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
		collect_shorter=[]
		status='accept'
		#print(acc,seq_dict[acc])
		for i in seq_dict[acc]:
			#collect_i.append(i)

			db_start=i[0]
			db_stop=i[1]
			db_len=i[4]

			if ( int(str_srna) <= int(db_start) < int(stop_srna))  or ( int(str_srna) < int(db_stop) <= int(stop_srna)) :
				if len(seq_rec) > db_len:
					status='accept'
				else:
					collect_i.append(i)
					status='Unaccept'
			else:
				collect_i.append(i) #add existing genes to db
				
		if status=='accept' :
			collect_i.append([str_srna,stop_srna,seq_rec.description,seq_rec.seq,len(seq_rec)])
			del seq_dict[acc]
			seq_dict[acc]=collect_i

		#print(seq_dict[acc])

		#print('Update Redundant',acc,len(seq_dict[acc]))
	#seq_dict[seq_rec.id]=seq_rec.seq
	#key.append(seq_rec.id)

w2file=open(sys.argv[2],'w')
for k in seq_dict:
	for j in seq_dict[k]:
		w2file.write('>'+j[2]+'\n'+str(j[3])+'\n')
w2file.close()

#usage
#python remove_redundant_seq.py ../test_dup.fasta outtest.fasta
