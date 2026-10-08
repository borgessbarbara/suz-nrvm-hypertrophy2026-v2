reference_fasta=/media/storage2/jpablo/Rattus_Norvegicus/rat_primary_assembly_orderr8.fa
reference_annotation=/media/storage2/surquiza_docs/GTFs_files/rat/Rattus_norvegicus.GRCr8.114.gtf
file=("H1" "H2" "H3" "H4" "SQ1" "SQ2" "SQ3" "SQ4" "SQ5" "SQ6" "SQ7" "SQ8" "Q13" "Q14" "Q15" "Q16") 
#file=("Q13" "Q14" "Q15" "Q16") 
for sample in ${file[@]} 
      do 
       #java -jar /media/storage2/jpablo/CIRI2/CIRI-vis_v1.4.jar -i /media/storage1/jpablo/CIRI2/CIRI_v2.0.6/AS_results/${sample}.ASout_jav.list -l /media/storage2/jpablo/CIRI2/CIRI_v2.0.6/AS_results/${sample}.ASout_library_length.list -r $reference_fasta -list /media/storage2/jpablo/CIRI2/CIRI_v2.0.6/cardio_bsj.txt -d ${sample}_vis -o ${sample}_vis -min 1 -type pdf
        java -jar /media/storage1/jpablo/circRNAs/CIRI2/CIRI-vis_v1.4.jar -i /media/storage2/jpablo/CIRI2/CIRI_v2.0.6/AS_results/${sample}.ASout_jav.list -l /media/storage2/jpablo/CIRI2/CIRI_v2.0.6/AS_results/${sample}.ASout_library_length.list -r $reference_fasta -list /media/storage2/jpablo/CIRI2/CIRI_v2.0.6/all_circRNA_DE.list -d ${sample}_vis -o ${sample}_vis -min 1 -type pdf
      rm ${sample}_vis/*.pdf
      done
     