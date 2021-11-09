
import subprocess 
from subprocess import PIPE,STDOUT,Popen
from ete3 import NCBITaxa
ncbi = NCBITaxa()

#need to install 'sudo apt-get install ncbi-entrez-direct'
#or 'conda install -c bioconda entrez-direct'

#use build_tree conda

def call_taxonomy(accession,tax_level): #accession e.g. NC_00911
	cmd='esearch -db nucleotide -query "'+str(accession)+'"|esummary|xtract -pattern TaxId -element TaxId'
	print(cmd)
	ps=Popen(cmd,shell=True, stdout=PIPE, encoding='utf8') #use utf8 to remove b'..'
	output = ps.stdout.read().rstrip()
	print(output)

	if output != '':
		tax_list=[]
		tax_level_found='no'
		for lineage in ncbi.get_lineage(output) : #lineage is taxonomy id for each level

			if tax_level == ncbi.get_rank([lineage]).get(lineage):
				tax_level_found='yes'
				ranking=ncbi.get_rank([lineage]).get(lineage)
				name_rank=str(ncbi.get_taxid_translator([lineage]).get(lineage))
				
				tax_ranking=ranking+' : '+name_rank+' ('+str(lineage)+')'
		if tax_level_found=='no':
			tax_ranking='N/A'
		print(tax_ranking)
	else:
		print('Query ID',accession,' not found!!')
call_taxonomy('NC_000912','species') #tax_level follows ncbi taxonomy: e.g. superkingdom, phylum, order, family, genus, species



