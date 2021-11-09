from Bio import SeqIO
import sys
gb_filename = str(sys.argv[1])
out_filename = str(sys.argv[2])

#usage: python gbk_faa.py GCF_000009725.1_ASM972v1_genomic.gbff NC_000911.faa
w2file=open(out_filename,'w')

for seq_record in SeqIO.parse(gb_filename, "genbank") :
	#print("Dealing with GenBank record %s") % seq_record.id
	for seq_feature in seq_record.features :
		if seq_feature.type=="CDS" :
			try:
				if 'join' in str(seq_feature.location):
					if seq_feature.location.strand == 1: #fwd strand
						if 0 in (seq_feature.location.parts[0].start,seq_feature.location.parts[1].start,seq_feature.location.parts[0].end,seq_feature.location.parts[1].end):
							w2file.write('>'+seq_feature.qualifiers['locus_tag'][0]+' coordinate=ooi:'+str(seq_feature.location.parts[0].start)+'_'+str(seq_feature.location.parts[1].end)+'_fwd'+'\n')
							w2file.write(seq_feature.qualifiers['translation'][0]+'\n')
							#print(seq_feature.location.parts[0].start,seq_feature.location.parts[1].start)
							#print(seq_feature.location.parts[0].end,seq_feature.location.parts[1].end)
							#print(seq_feature.qualifiers['locus_tag'][0],seq_feature.location)
						else:
							w2file.write('>'+seq_feature.qualifiers['locus_tag'][0]+' coordinate=ooi:'+str(seq_feature.location.start)+'_'+str(seq_feature.location.end)+'_fwd'+'\n')
							w2file.write(seq_feature.qualifiers['translation'][0]+'\n')
					else: #rev strand
						if 0 in (seq_feature.location.parts[0].start,seq_feature.location.parts[1].start,seq_feature.location.parts[0].end,seq_feature.location.parts[1].end):
							w2file.write('>'+seq_feature.qualifiers['locus_tag'][0]+' coordinate=ooi:'+str(seq_feature.location.parts[0].end)+'_'+str(seq_feature.location.parts[1].end)+'_rev'+'\n')
							w2file.write(seq_feature.qualifiers['translation'][0]+'\n')
						else:
							w2file.write('>'+seq_feature.qualifiers['locus_tag'][0]+' coordinate=ooi:'+str(seq_feature.location.start)+'_'+str(seq_feature.location.end)+'_rev'+'\n')
							w2file.write(seq_feature.qualifiers['translation'][0]+'\n')
				#w2file.write('>'+seq_feature.qualifiers['locus_tag'][0]+'\n')
				#print(seq_feature.location)
				#print(seq_feature.qualifiers['locus_tag'][0])
				#print(seq_feature.qualifiers['translation'][0])
				else:
					#print(seq_feature.location)
					if seq_feature.location.strand == 1:
						w2file.write('>'+seq_feature.qualifiers['locus_tag'][0]+' coordinate=ooi:'+str(seq_feature.location.start)+'_'+str(seq_feature.location.end)+'_fwd'+'\n')
						w2file.write(seq_feature.qualifiers['translation'][0]+'\n')
					else:
						w2file.write('>'+seq_feature.qualifiers['locus_tag'][0]+' coordinate=ooi:'+str(seq_feature.location.start)+'_'+str(seq_feature.location.end)+'_rev'+'\n')
						w2file.write(seq_feature.qualifiers['translation'][0]+'\n')
			except KeyError:
				print(seq_feature.qualifiers['locus_tag'][0])
				print(seq_feature.qualifiers['note'][0])
				
w2file.close()
print('Finish!!!')
