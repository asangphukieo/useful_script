from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation
import sys

start=sys.argv[1]
stop=sys.argv[2]
strand=sys.argv[3]
genome_record = SeqIO.read("Chromosome_with_new_ORFs.genbank", "genbank")
#print(genome_record.seq)
feature_loc = FeatureLocation(int(start), int(stop), strand=int(strand))


print(feature_loc.extract(genome_record.seq))

#ncl0310 complement(764097..764234)
#sml0013 complement(3188105..3188227)
