import commands

list_uniprot=commands.getoutput('cut -f1 uniq_all_specific_photo.out')
list_uniprot=list_uniprot.split('\n')

for i in list_uniprot:
	check=commands.getoutput('grep '+i+' all_seq4.faa')
	#print i
	if check == '':
		print i,' No in fasta'
	
