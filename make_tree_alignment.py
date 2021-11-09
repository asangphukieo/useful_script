# Import Python libraries
from ete3 import PhyloTree
import sys

#input1="SR1P/clustalo_default-none-none-fasttree_full/SR1P_blastout.fasta.final_tree.nw"
#input2="SR1P/clustalo_default-none-none-fasttree_full/SR1P_blastout.fasta.final_tree.used_alg.fa"
#output="rdp_alignment_tree.pdf"

#why does it not work?
input1=str(sys.argv[1])
input2=str(sys.argv[2])
output=str(sys.argv[3])

# Import the tree using PhyloTree class
tree1 = PhyloTree(input1)

# Add alignment
tree1.link_to_alignment(input2)
tree1.render(file_name=output )
#tree1.render("%%inline")
# Render in jupyter notebook
#tree1.render("%%inline", h=150, units="mm", dpi=100)

##usage
#python make_tree_alignment.py SR1P_nt/clustalo_default-none-none-fasttree_full/coding_sequence.fasta.final_tree.nw SR1P_nt/clustalo_default-none-none-fasttree_full/coding_sequence.fasta.final_tree.used_alg.fa coding_SR1P.pdf
