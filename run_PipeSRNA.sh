#!/bin/bash


# -----------
# parameter file
# -----------
source ./PARAMETER_FILE/parameters.conf

GLASSgo_OUT_FOLDER='01_GLASSgo_Results'
GLASSgo2CopraRNA="02_GLASSgo2CopraRNA"
CopraRNA_OUT_FOLDER="03_CopraRNA_Results"
MAFFT_OUT_FOLDER="04_Mafft_Results"
RNAalifold_OUT_FOLDER="05_RNAalifold_Results"
Rscape_OUT_FOLDER="06_Rscape_Results"
RNAcode_OUT_FOLDER="07_RNAcode_Results"
Synteny_OUT_FOLDER="08_Synteny_Results"
Promoter_OUT_FOLDER="09_Promoter_Results"

# ------------
# Split multi fasta file into single sequence fasta file
# ------------


if  [ $DO_SPLIT_INPUT == 'Y' ]; then
    if [ -d "./SRNA_INPUT/SPLIT_FASTA" ] ; then
         echo 'WARNING: Folder SPLIT_FASTA exists! , Please remove before start running, otherwise the old sequences will be used'
         INPUT_FOLDER="./SRNA_INPUT/SPLIT_FASTA"
    else
         python split_multifasta.py ./SRNA_INPUT/${INPUT} ./SRNA_INPUT/SPLIT_FASTA
         INPUT_FOLDER="./SRNA_INPUT/SPLIT_FASTA" 
    fi
else
    INPUT_FOLDER="./SRNA_INPUT"
fi

# -----------
# run GLASSgo
# -----------

if  [ $DO_GLASSgo == 'Y' ]; then
    mkdir -p $GLASSgo_OUT_FOLDER
    source ./activate aligners

    for f in `ls $INPUT_FOLDER`; do
        if [[ $f =~ .*\.(fasta|fa|fas) ]]; then
        echo "Running GLASSgo for sRNA: " $f
        fullpath=${INPUT_FOLDER}/$f
        output_sRNAs=GLASSgo_output_${f}
        output_json=${output_sRNAs%.*}.json

        if  [ $GI_MODE == 'Y' ]; then   
            gi='-g '$GI_FILE
        else
            gi=''
        fi

        ./GLASSgo.py -d $nt_DATABASE_PATH -l $LONDEN $gi -p $MIN_PIDENT -e $E_value_GG -i $fullpath -t $CPU_CORE_GG -o $output_sRNAs  
        ./createJSON-1_0.py -g $output_sRNAs -o $output_json
        ./generate_GLASSgo_html.py $output_json
        fi
    done
    mv GLASSgo_output_* $GLASSgo_OUT_FOLDER

    conda deactivate
fi


if  [ $DO_REPETE_GLASSgo == 'Y' ]; then
    
    source ./activate aligners

    for f in `ls $GLASSgo_OUT_FOLDER`; do
        if [[ $f =~ .*\.(fasta|fa|fas) ]]; then
            echo "REPETE: Running GLASSgo for the most similar sRNA: " $f
            fullpath=${GLASSgo_OUT_FOLDER}/$f

            output_sRNAs=${f}
            output_json=${output_sRNAs%.*}.json

            #select the most similar sequence
            python ./select_closest_seq.py $fullpath GLASS_go_query_repete.fasta

            if  [ $GI_MODE == 'Y' ]; then   
                gi='-g '$GI_FILE
            else
                gi=''
            fi

            ./GLASSgo.py -d $nt_DATABASE_PATH -l $LONDEN $gi -p $MIN_PIDENT -e $E_value_GG -i GLASS_go_query_repete.fasta -t $CPU_CORE_GG -o rep1_$output_sRNAs
            rm GLASS_go_query_repete.fasta
            python ./combine_seq.py $fullpath rep1_$output_sRNAs $output_sRNAs
            ./createJSON-1_0.py -g $output_sRNAs -o $output_json
            ./generate_GLASSgo_html.py $output_json
            rm rep1_$output_sRNAs
        fi
    done
    mkdir ${GLASSgo_OUT_FOLDER}/ALL_SEQ_NO_REPETE
    mv ${GLASSgo_OUT_FOLDER}/* ${GLASSgo_OUT_FOLDER}/ALL_SEQ_NO_REPETE/
    mv GLASSgo_output_* $GLASSgo_OUT_FOLDER

    conda deactivate
fi


if  [ $GI_FILTER == 'Y' ]; then     
    source ./activate aligners       
    for k in `ls $GLASSgo_OUT_FOLDER`; do
        path2GG=$GLASSgo_OUT_FOLDER/$k
        if [[ $k =~ .*\.(fasta|fa|fas) ]]; then 
            python3 filter_taxonomy.py $path2GG ${GI_ID}
            ./createJSON-1_0.py -g $path2GG -o ${k%.*}.json
            ./generate_GLASSgo_html.py ${k%.*}.json

        fi
    done
    mkdir ${GLASSgo_OUT_FOLDER}/ALL_SEQ
    mv ${GLASSgo_OUT_FOLDER}/*.fasta_x ${GLASSgo_OUT_FOLDER}/ALL_SEQ/
    mv ${GLASSgo_OUT_FOLDER}/*.json ${GLASSgo_OUT_FOLDER}/ALL_SEQ/
    mv ${GLASSgo_OUT_FOLDER}/*.html ${GLASSgo_OUT_FOLDER}/ALL_SEQ/

    mv GLASSgo_output_* $GLASSgo_OUT_FOLDER
    conda deactivate
fi

if  [ $GI_FIXHEADER == 'Y' ]; then     
    source ./activate r342
    for x in $GLASSgo_OUT_FOLDER/*.fasta ;do

         python fix_header_glassgo2copra.py $x ${x}2
    done
    mkdir ${GLASSgo_OUT_FOLDER}/ALL_SEQ_GI_HEADER
    mv ${GLASSgo_OUT_FOLDER}/*.fasta ${GLASSgo_OUT_FOLDER}/ALL_SEQ_GI_HEADER/
    for x in $GLASSgo_OUT_FOLDER/*.fasta2 ; do mv $x ${x%.*}.fasta;done 
    conda deactivate 
fi

# --------------------
# run GLASSgo2CopraRNA 
# --------------------

if  [ $DO_GLASSgo2CopraRNA == 'Y' ]; then
    mkdir -p $GLASSgo2CopraRNA

    source ./activate r342

    # update CopraRNA organism list
    if  [ $UPDATE_CP_ORGANISMS == 'True' ]; then
        echo " Update organism list in CopraRNA ... "
        mkdir update_org && cd update_org && perl ../build_kegg2refseq.pl && cd ..
        
        mv ~/anaconda2/envs/JensCopraRNA2/bin/coprarna_aux/CopraRNA_available_organisms.txt ~/anaconda2/envs/JensCopraRNA2/bin/coprarna_aux/CopraRNA_available_organisms.txt.previous
        mv ~/anaconda2/envs/JensCopraRNA2/bin/coprarna_aux/kegg2refseqnew.csv ~/anaconda2/envs/JensCopraRNA2/bin/coprarna_aux/kegg2refseqnew.csv.previous
        mv ./update_org/CopraRNA_available_organisms.txt ~/anaconda2/envs/JensCopraRNA2/bin/coprarna_aux/
        mv ./update_org/kegg2refseqnew.csv ~/anaconda2/envs/JensCopraRNA2/bin/coprarna_aux/
        cp ~/anaconda2/envs/JensCopraRNA2/bin/coprarna_aux/CopraRNA_available_organisms.txt ./
        rm -r update_org

    fi

    echo '#Input   GLASSgo2CopraRNA_status' > GLASSgo2CopraRNA_status.out
    for f in `ls $GLASSgo_OUT_FOLDER`; do
        if [[ $f =~ .*\.(fasta|fa|fas) ]]; then 

            echo "Running GLASSgo2CopraRNA for sRNA: " $f
            fullpath=$GLASSgo_OUT_FOLDER/$f
            filestem=${f%.*}
            ID=${filestem#GLASSgo_output_}

            if [ -d $GLASSgo2CopraRNA/$ID ]; then    
                echo "File found! >> " $f
            else

                #updated GLASSgo_postprocessing 7.0
                if R --slave -f  GLASSgo_postprocessing_8_modi.r --args filename=$fullpath duplicates_allowed=FALSE synteny_window=10000 name=${ID} coprarna_compatible=TRUE ooi=$ORG_REF ; then # TRY
                    echo $f'    Pass' >> GLASSgo2CopraRNA_status.out
                else # CATCH
                    echo 'Number of organisms is too low'
                    echo $f'    Not pass' >> GLASSgo2CopraRNA_status.out

                fi            
                
              
                # fasta generation
                #R --slave -f  GLASSgo2CopraRNA_fasta_generation_by_accession.r --args input_file=$fullpath refpath=taxid_to_refseq cop_path=CopraRNA_available_organisms.txt output_file=${ID}_coprarna_candidates.txt
                # exclusion 
                #R --slave -f  GLASSgo2CopraRNA_exclusion_script.r --args datapath=full_GLASSgo_table.Rdata wildcard=$WILDCARD ooi=$ORG_REF
                # balanced selection
                #R --slave -f  GLASSgo2CopraRNA_balanced_ooi_selection.r --args wildcard=$WILDCARD max_number=20 outfile_prefix=${ID}_sRNA ooi=$ORG_REF
                # move results to its folder
                mkdir -p $GLASSgo2CopraRNA/$ID
                
                #remove sequences of not available in CopraRNA
                #python2 Remove_unmatch_organism.py $GLASSgo2CopraRNA/${ID}/*.fasta CopraRNA_available_organisms.txt
                rm full_tree.Rdata distances.Rdata coor.Rdata 16s_subset.fasta 16s_aligned2.fasta 16s_aligned.fasta psi_out -r
                mv 16s.fasta $GLASSgo2CopraRNA/$ID
                mv ${ID}_* $GLASSgo2CopraRNA/$ID
		mv GLASSgo2CopraRNA_status.out $GLASSgo2CopraRNA/$ID
		mv pass_glassgo2copra.fasta $GLASSgo2CopraRNA/$ID
		mv aligned_sRNA.fasta $GLASSgo2CopraRNA/$ID


            fi

        fi
    done

    #source deactivate

    conda deactivate
fi

if  [ $DO_CHECK_AND_REMOVE_LOW_MATCH == 'Y' ]; then

for file in ./*.pdf####################
do
    if [ -f "${file}" ]; then
    echo 'true';
    break
    fi
done

fi

if  [ $DO_CHECK_REPETE_GLASSgo == 'Y' ]; then
    
    source ./activate aligners

    for f in `ls $GLASSgo_OUT_FOLDER`; do
        if [[ $f =~ .*\.(fasta|fa|fas) ]]; then
            echo "REPETE: Running GLASSgo for the most similar sRNA: " $f
            fullpath=${GLASSgo_OUT_FOLDER}/$f

            output_sRNAs=${f}
            output_json=${output_sRNAs%.*}.json

            #select the most similar sequence
            python ./select_closest_seq.py $fullpath GLASS_go_query_repete.fasta

            if  [ $GI_MODE == 'Y' ]; then   
                gi='-g '$GI_FILE
            else
                gi=''
            fi

            ./GLASSgo.py -d $nt_DATABASE_PATH -l $LONDEN $gi -p $MIN_PIDENT -e $E_value_GG -i GLASS_go_query_repete.fasta -t $CPU_CORE_GG -o rep1_$output_sRNAs
            rm GLASS_go_query_repete.fasta
            python ./combine_seq.py $fullpath rep1_$output_sRNAs $output_sRNAs
            ./createJSON-1_0.py -g $output_sRNAs -o $output_json
            ./generate_GLASSgo_html.py $output_json
            rm rep1_$output_sRNAs
        fi
    done
    mkdir ${GLASSgo_OUT_FOLDER}/ALL_SEQ_NO_REPETE
    mv ${GLASSgo_OUT_FOLDER}/* ${GLASSgo_OUT_FOLDER}/ALL_SEQ_NO_REPETE/
    mv GLASSgo_output_* $GLASSgo_OUT_FOLDER

    conda deactivate
fi


# ------------
# run CopraRNA XXXXXXXXXXXXXXXXXXXXXX
# ------------
if  [ $DO_CopraRNA == 'Y' ]; then
    mkdir -p $CopraRNA_OUT_FOLDER


    # parallel run
    source ./activate JensCopraRNA2
    #./parallel_CopraRNA_prediction.py $GLASSgo2CopraRNA --suffix .fasta --cores $CPU_CORE_CP --batch 20 --out_folder $CopraRNA_OUT_FOLDER \
    #              --ntup $NTUP --ntdown $NTDOWN --region 5utr --enrich 200 --topcount 200 --websrv --noclean 2>&1 | tee ./parallel_CopraRNA_log_file.txt

    #### run each sequence individually ####
    for f in `ls $GLASSgo2CopraRNA`; do
        for k in `ls $GLASSgo2CopraRNA/$f`; do
            if [[ $k == *.fasta ]]; then 
                echo "Running CopraRNA ... >> " $k
                echo "Path input >>" $GLASSgo2CopraRNA/$f/$k

                mkdir $CopraRNA_OUT_FOLDER/$f
                cd $CopraRNA_OUT_FOLDER/$f
                CopraRNA2.pl -srnaseq /working_directory/$GLASSgo2CopraRNA/$f/$k -ntup $NTUP -ntdown $NTDOWN -region 5utr -enrich 200 -topcount 200 -websrv -noclean -cores $CPU_CORE_CP
                cd /working_directory
            fi
        done
    done
    #### End run each sequence individually #### 

    conda deactivate
fi

if  [ $DO_CopraRNAsimplify == 'Y' ]; then

    # ----
    # simplify CopraRNA_results
    # ----

    # prepare RefSeq files
    if  [ $PREPARE_REFSEQ == 'True' ]; then

        if  [ ${ORG_REF} != ${ORG_RUN} ]; then
            python pair_identifiers_by_coordinates_and_identities.py ./option_files/${ORG_REF}.faa ./option_files/${ORG_RUN}.faa -p RefSeq2org
        else
            grep '^>' option_files/${ORG_REF}.faa |cut -d' ' -f1|sed s/'>'/''/g>artifact.temp
            paste -d '\t' artifact.temp artifact.temp > RefSeq2org.tsv
            rm artifact.temp
        fi        
        

        if  [ $DOWNLOAD_INTERPRO == 'True' ]; then
        
            #download interproscan

            if [ -d "interproscan-${INTERPRO_VERSION}" ]; then
                echo InterPro folder exists!
                PATH_INTERPRO=./interproscan-${INTERPRO_VERSION}
            else
                wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${INTERPRO_VERSION}/interproscan-${INTERPRO_VERSION}-64-bit.tar.gz
                tar xvzf interproscan-${INTERPRO_VERSION}-64-bit.tar.gz --directory=./
                PATH_INTERPRO=./interproscan-${INTERPRO_VERSION}           
            fi
        fi

        #set path to interproscan in parameters.conf
        ${PATH_INTERPRO}/interproscan.sh -i /working_directory/option_files/${ORG_REF}.faa -d /working_directory -goterms -pa

        python ipr2goatools.py ${ORG_REF}.faa.tsv
        python split_GO.py ${ORG_REF}_ipr2goatools.tsv ${ORG_REF}_ipr2go.tsv
    fi

    python ./simplify_CopraRNA_results.py $CopraRNA_OUT_FOLDER ./RefSeq2org.tsv
    source ./activate r342
    simplified_folder=03_CopraRNA_Results_simplified
    cd $simplified_folder
    for f in `ls *_locustag.txt`; do
        # GO enrichment
	source ../activate GOstat ##need to extract the conda file
        R --slave -f ../GOstats_analysis.R --args workdir=. organism=${ORG_NAME} target_gene_file=$f ipr2go_file=../${ORG_REF}_ipr2go.tsv 

        source ../activate r342
        # fetch revigo csv
        filename=`basename $f`
        filestem=${filename%.txt}
        ../fetch_revigo_csv.py ${filestem}_GO_enrichment_BP.tsv -o ./
        ../fetch_revigo_csv.py ${filestem}_GO_enrichment_CC.tsv -o ./
        ../fetch_revigo_csv.py ${filestem}_GO_enrichment_MF.tsv -o ./
        # revigo visualization
        R --slave -f ../REVIGO_plotter.R --args workdir=. REVIGO_BP=${filestem}_GO_enrichment_BP_revigo.csv REVIGO_CC=${filestem}_GO_enrichment_CC_revigo.csv REVIGO_MF=${filestem}_GO_enrichment_MF_revigo.csv out_pdf=${filestem}_GO_enrichment.pdf
    done
    conda deactivate
    cd ..
fi
#R --slave -f REVIGO_plotter.R --args workdir=. REVIGO_BP= REVIGO_CC= REVIGO_MF= out_pdf=test_GO_enrichment.pdf

# -----------
# run mafft
# -----------

if  [ $DO_MAFFT == 'Y' ]; then
    mkdir -p $MAFFT_OUT_FOLDER

    source ./activate aligners
    unset MAFFT_BINARIES

    for ID_Folder in `ls $GLASSgo2CopraRNA`; do
        for f in `ls $GLASSgo2CopraRNA/$ID_Folder`; do
            if [[ $f =~ _selection.*\.(fasta|fa|fas) ]]; then 
                echo "Running mafft-linsi for sRNA: " $f
                fullpath=$GLASSgo2CopraRNA/$ID_Folder/$f
                filestem=${f%.*}
                ID=${filestem#GLASSgo_output_}
                # run mafft-linsi
                mafft --localpair --maxiterate 1000 --reorder --thread 20 $fullpath > ${ID}_mafft.fasta
                # run mview
                mview -in fasta -html head -css on -coloring id -colormap clustal ${ID}_mafft.fasta > ${ID}_mafft.html
                # run fasttree
                cat ${ID}_mafft.fasta | tr : _ > ${ID}_modified.fa       
                fasttree -quiet -gtr -gamma -nt ${ID}_modified.fa > ${ID}_fasttree.nwk
                rm ${ID}_modified.fa
                #figtree  -graphic  PDF ${ID}_fasttree.nwk ${ID}_fasttree.pdf

                # move results to MAFFT folder
                mv ${ID}_mafft.* ${ID}_fasttree.* $MAFFT_OUT_FOLDER
            fi
        done
    done
    conda deactivate
fi


# ---------------
# run RNAalifold
# ---------------

if  [ $DO_RNAalifold == 'Y' ]; then
    mkdir -p $RNAalifold_OUT_FOLDER
    source ./activate py27
    for f in `ls $MAFFT_OUT_FOLDER`; do
        echo $f
        if [[ $f =~ .*\.(fasta) ]]; then 
            echo "Running RNAalifold for sRNA: " $f
            fullpath=$MAFFT_OUT_FOLDER/$f
            filestem=${f%.*}
            ID=${filestem%_mafft.fasta}
            # shorten header
            ./shorten_GLASSgo_fasta.py --force $fullpath -p $ID 
            # convert to sto and aln format
            seqmagick convert ${ID}_shorten.fa ${ID}_shorten.sto
            seqmagick convert ${ID}_shorten.fa ${ID}_shorten.aln
            # run RNAalifold and tidy output files
            RNAalifold --color --aln --sci ${ID}_shorten.aln > ${ID}_RNAalifold_output.txt
            mv aln.ps ${ID}_aln.ps
            mv alirna.ps ${ID}_alirna.ps
            ps2pdf -dEPSCrop ${ID}_aln.ps ${ID}_aln.pdf 
            ps2pdf -dEPSCrop ${ID}_alirna.ps ${ID}_alirna.pdf
            # prepare Rscape input 
            ./modify_sto_file_for_Rscape.py ${ID}_shorten.sto ${ID}_RNAalifold_output.txt
            mv ${ID}_* $RNAalifold_OUT_FOLDER
        fi
    done
    conda deactivate
fi

# ------------
# run R-scape 
# ------------

if  [ $DO_RSCAPE == 'Y' ]; then
    source ./activate RNAtools

    mkdir -p $Rscape_OUT_FOLDER
    for f in `ls $RNAalifold_OUT_FOLDER`; do
        if [[ $f =~ .*\.sto ]]; then 
            echo "Running R-scape for sRNA: " $f
            fullpath=$RNAalifold_OUT_FOLDER/$f
            filestem=${f%.*}
            ID=${filestem%_mafft_shorten_for_Rscape.sto}
            # run Rscape
            R-scape -E 0.1 --outmsa $Rscape_OUT_FOLDER/${ID}_Rscape_msa.fa -o $Rscape_OUT_FOLDER/${ID}_Rscape_output.txt --r2rall --outdir $Rscape_OUT_FOLDER $fullpath
        fi
    done
fi

# ------------
# run RNAcode 
# ------------

if  [ $DO_RNACODE == 'Y' ]; then
    source ./activate RNAtools
    mkdir -p $RNAcode_OUT_FOLDER

    for f in `ls $RNAalifold_OUT_FOLDER`; do
        if [[ $f =~ .*\.aln ]]; then 
            echo "Running RNAcode for sRNA: " $f
            fullpath=$RNAalifold_OUT_FOLDER/$f
            filestem=${f%.*}
            ID=${filestem%_sRNA_CopraRNA_input_balanced*.aln}
            # run RNAcode
            mkdir -p $ID
            RNAcode --gtf --tabular --eps --eps-dir $ID --eps-cutoff 0.1 --cutoff 0.1 --outfile ${ID}_RNAcode.tab $fullpath

            if [ -z "$(ls -A $ID)" ]; then
                echo "Protein coding regions were not found!"
                rm ${ID}_*
                rm -rf $ID
            else
                echo "Coding regions detected!"
                cd $ID
                for f in `ls .`; do
                    if [[ $f =~ .*\.(eps|ps) ]]; then
                        ps2pdf -dEPSCrop $f
                    fi
                done
                cd ../
                mv ${ID}_* $ID
                mv $ID $RNAcode_OUT_FOLDER
            fi

        fi
    done
    conda deactivate
fi


# ------------
# run synteny 
# ------------

if  [ $DO_SYNTENY == 'Y' ]; then
    source ./activate r342
    mkdir -p $Synteny_OUT_FOLDER

    # copy sRNA to Synteny_OUT_FOLDER
    echo "Copying GLASSgo output sRNAs to $Synteny_OUT_FOLDER ..."
    cp $GLASSgo_OUT_FOLDER/*.{fa,fas,fasta} $Synteny_OUT_FOLDER
    #rm $Synteny_OUT_FOLDER/*.fa_old
    cd $Synteny_OUT_FOLDER
    cp ../synteny_pdf.r .
    for f in `ls .`; do
        if [[ $f =~ .*\.(fasta|fa|fas) ]]; then 
            echo "Running synteny for sRNA: " $f
            filestem=${f%.*}
            ID=${filestem#GLASSgo_output_}
            if [ ! -f ${ID}_synteny.pdf ]; then
                R --slave -f synteny_pdf.r --args input_sRNA=$f output_prefix=${ID} 2>&1 | tee ./${ID}_synteny_log.txt
            fi
        fi
    done
    conda deactivate
    rm synteny_pdf.r
    cd ../
fi


# -------------------------------------
# run motif_finding in promoter regions 
# -------------------------------------
if  [ $DO_PROMOTER == 'Y' ]; then

    mkdir -p $Promoter_OUT_FOLDER
    source ./activate r342 
    for f in `ls $GLASSgo_OUT_FOLDER`; do
        

        if [[ $f =~ .*\.(fasta|fa|fas) ]]; then 
            echo "Fetching promoter sequences for sRNA $f ... "
            fullpath=$GLASSgo_OUT_FOLDER/$f
            filestem=${f%.*}
            ID=${filestem#GLASSgo_output_}
            # fasta generation
            R --slave -f  GLASSgo2CopraRNA_fasta_generation_by_accession.r --args input_file=$fullpath refpath=taxid_to_refseq cop_path=CopraRNA_available_organisms.txt output_file=${ID}_coprarna_candidates.txt
            # exclusion 
            R --slave -f  GLASSgo2CopraRNA_exclusion_script.r --args datapath=full_GLASSgo_table.Rdata wildcard=$WILDCARD ooi=$ORG_REF
            # fetch promoters
            R --slave -f promoter_sequence_fetcing_from_coor.r --args datapath=refined_GLASSgo_table.Rdata output_prefix=${ID}_promoter
            # move results to its folder
            mkdir -p $Promoter_OUT_FOLDER/$ID
            mv ${ID}_* $Promoter_OUT_FOLDER/$ID

     
        fi
    done
    conda deactivate

    source ./activate RNAtools
    for ID_FOLDER in `ls $Promoter_OUT_FOLDER`; do
        meme $Promoter_OUT_FOLDER/$ID_FOLDER/*_promoter.fasta -oc $Promoter_OUT_FOLDER/$ID_FOLDER/MEME_results -dna -mod anr -nmotifs 50 -revcomp -evt 0.05 -minw 4 -maxw 200 -maxsize 10000000000000
    done
    conda deactivate

fi


if  [ $DO_SUM_HTML == 'Y' ]; then
    #move output to mounted folder
    python summarize_to_html.py ${INPUT}
    mv summarized_information.html ./OUTPUT
fi

if  [ $DO_MOVEOUTPUT == 'Y' ]; then
    #move output to mounted folder
    mv 0*_* ./OUTPUT
    mv PARAMETER_FILE/parameters.conf ./OUTPUT
fi

