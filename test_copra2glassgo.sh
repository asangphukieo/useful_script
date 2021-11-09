source ./activate r342
GLASSgo_OUT_FOLDER='01_GLASSgo_Results'
GLASSgo2CopraRNA="02_GLASSgo2CopraRNA"


echo '#Input\tGLASSgo2CopraRNA_status' > GLASSgo2CopraRNA_status.out

for f in `ls $GLASSgo_OUT_FOLDER`; do
    if [[ $f =~ .*\.(fasta|fa|fas) ]]; then 
        echo "Running GLASSgo2CopraRNA for sRNA: " $f
        fullpath=$GLASSgo_OUT_FOLDER/$f
        filestem=${f%.*}
        ID=${filestem#GLASSgo_output_}

	if R --slave -f  GLASSgo_postprocessing_7.r --args filename=$fullpath duplicates_allowed=FALSE synteny_window=3000 name=${ID} coprarna_compatible=TRUE ooi=$ORG_REF ; then # TRY
            echo $f'\tPass' >> GLASSgo2CopraRNA_status.out 
	else # CATCH
	    echo 'Number of organisms is too low'
            echo $f'\tNot pass' >> GLASSgo2CopraRNA_status.out
        fi
        mkdir -p $GLASSgo2CopraRNA/$ID
        mv ${ID}_* $GLASSgo2CopraRNA/$ID

        rm full_tree.Rdata distances.Rdata coor.Rdata 16s_subset.fasta 16s_aligned2.fasta 16s_aligned.fasta psi_out -r
        mv 16s.fasta $GLASSgo2CopraRNA/$ID 

    fi
done
conda deactivate
