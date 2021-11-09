#written by Apiwat Sangphukieo , apiwat.bif@mail.kmutt.ac.th

import sys
w2file=open(sys.argv[2],'w')

for i in open(sys.argv[1]):
	i=i.replace('\n','')
	col=i.split()
	gos=col[1].split(';')
	for j in gos:
		w2file.write(col[0]+'\t'+j+'\n')
		#print (col[0],j)

w2file.close()	

#Usage python split_GO.py input_ipr2goatools.tsv output_ipr2go.tsv
