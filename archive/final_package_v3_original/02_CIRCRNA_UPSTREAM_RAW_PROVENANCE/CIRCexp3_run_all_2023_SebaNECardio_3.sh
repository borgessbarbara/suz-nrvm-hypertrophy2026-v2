reference_fasta=/media/storage2/jpablo/Rattus_Norvegicus/rat_primary_assembly_orderr8.fa
reference_annotation=/media/storage2/jpablo/Rattus_Norvegicus/rn8.refflat.nochr.txt
STAR_index=/media/storage2/surquiza_docs/index/STAR/Rat/STAR_rn8/
hisat_index=/media/storage2/surquiza_docs/genomes/rat/Ensembl/Rn8
export PATH=$PATH:/media/storage2/software/hisat2-2.2.1

for sample in {13..16}
do
    reads="/media/storage2/surquiza_docs/datasets_rat_HC/dataset_NE_cardios/lncRNA/raw_data_fastq/SQ${sample}_ln*.fq.gz"

    zip_reads=$reads
    echo "Seq Files..."
    echo $zip_reads
    #Decompress the two files a la carpeta donde se va a ejecutar el programa
    echo "unzpping files"
    for f in $zip_reads
    do
      STEM=$(basename "${f}" .gz)
      gunzip -k -c "${f}" > /media/storage2/jpablo/Circexp3/"${STEM}"
    done

    echo "Unziped files"

    #Create a variable to include the two unzip reads
    reads=/media/storage2/jpablo/Circexp3/*.fq
    echo $reads

     for read in $reads;
     do 
     if [[ "$read" =~ .*.ln_1.fq ]];
         then
             echo "Recognized read 1"
             read1="$read"
             echo $read1
     elif [[ "$read" =~ .*.ln_2.fq ]];
         then
             echo "Recognized read 2" 
             read2="$read"
             echo $read2
     fi
     done

     # rescue the the sample's name or number
     unset prefix_name_sample
     prefix_name_sample=$(basename "${read1: -10:-8}")
     echo $prefix_name_sample
    
     unset prefix_name_sample_sam
     prefix_name_sample_sam=$(basename "${read1: -10:-8}".sam)
    
     unset prefix_name_sample_ciri
     prefix_name_sample_ciri=$(basename "${read1: -10:-8}".ciri)
     echo $prefix_name_sample_sam
    
     echo "Sample"
     echo $prefix_name_sample

     echo "cut read length to the same length"
     #trimmomatic PE -phred33 -threads 4 $read1 $read2 $read1.trimmed.fq $read1.un.trimmed.fq  $read2.trimmed.fq $read2.un.trimmed.fq ILLUMINACLIP:/media/storage1/jpablo/circRNAs/CIRI2/z_adapters_sequencing.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 MINLEN:32 #CROP:140
     ../CIRIquant/fastp -i $read1 -I $read2 -o $read1.trimmed.fq -O $read2.trimmed.fq --detect_adapter_for_pe -f 10 -F 10 -3 --cut_tail_mean_quality 28 -q 30 -w 8
     echo "TRIMMO END"

        unset $read1
        unset $read2

        reads=/media/storage2/jpablo/Circexp3/*.fq
        echo $reads

        for read in $reads;
        do 
        if [[ "$read" =~ .*.ln_1.fq.trimmed.fq ]];
            then
                echo "Recognized read 1"
                read1="$read"
                echo $read1
        elif [[ "$read" =~ .*.ln_2.fq.trimmed.fq ]];
            then
                echo "Recognized read 2" 
                read2="$read"
                echo $read2
        fi
        done
        # rescue the the sample's name or number
        unset prefix_name_sample
        prefix_name_sample=$(basename "${read1: -21:-19}")
        echo $prefix_name_sample
        
        unset prefix_name_sample_sam
        prefix_name_sample_sam=$(basename "${read1: -21:-19}".sam)
        
        unset prefix_name_sample_ciri
        prefix_name_sample_ciri=$(basename "${read1: -21:-19}".ciri)
        echo $prefix_name_sample_sam
        echo "Sample"
        echo $prefix_name_sample

        unset CIRCRNA
        CIRCRNA="_CIRCRNA_"
        unset chimeric
        chimeric="Chimeric.out.junction"
        unset parse
        parse="_CIRCRNA_CIRCexplorer2_parse.log"
        unset circRNAknown
        circRNAknown="circularRNA_known.txt"
        unset annotate
        annotate="CIRCexplorer2_annotate.log"
        unset back_spliced_junction
        back_spliced_junction="back_spliced_junction.bed"

        unset unmapped_reads
        unmapped_reads=/media/storage2/jpablo/Circexp3/circexplorer_results/hisat/unmapped.fq
        unset align_bam
        align_bam=/media/storage2/jpablo/Circexp3/circexplorer_results/hisat/align.bam
        unset quant
        quant="quant.txt"

        hisat_out_dir=/media/storage2/jpablo/Circexp3/circexplorer_results

        unset prefix_name
        prefix_name=$prefix_name_sample$CIRCRNA
        unset prefix_CIRCexplorer2_log2
        prefix_CIRCexlorer2_log2=$prefix_name_sample$parse
        unset prefix_Chimeric_out_junction
        prefix_Chimeric_out_junction=$prefix_name$chimeric
        unset prefix_circularRNA_known_txt
        prefix_circularRNA_known_txt=$prefix_name$circRNAknown
        unset prefix_CIRCexplorer2_annotate_log
        prefix_CIRCexplorer2_annotate_log=$prefix_name$annotate
        unset prefix_back_spliced_junctions
        prefix_back_spliced_junction=$prefix_name$back_spliced_junction
        unset prefix_quant
        prefix_quant=$prefix_name$quant


        echo "starting hisat2...."
        # Hisat mapping
            python /media/storage1/jpablo/circRNAs/CIRCexp3/CLEAR/src/run3.py -1 $read1 -2 $read2 -g $reference_fasta -i $hisat_index -G $reference_annotation -o $hisat_out_dir -p 30
        echo "hisat2 finish"

        echo "starting STAR...."
        # Reads map to index
            STAR --chimSegmentMin 10 --runThreadN 20 --genomeDir $STAR_index --readFilesIn $unmapped_reads --outFileNamePrefix $prefix_name
        echo "Read mapped to index finish"
        #Find Chimeric out juctions
            CIRCexplorer2 parse -t STAR  $prefix_Chimeric_out_junction > $prefix_CIRCexlorer2_log2
        echo "Chimeric out junctions finish"
        #Anotacion
            CIRCexplorer2 annotate -r $reference_annotation -g $reference_fasta -b back_spliced_junction.bed -o $prefix_circularRNA_known_txt > $prefix_CIRCexplorer2_annotate_log
        echo "Annotation finish"

        ## quantification circexplorer3
        echo "quant circ3 start"
            python //media/storage1/jpablo/circRNAs/CIRCexp3/CLEAR/src/circ_quant.py -c /media/storage2/jpablo/Circexp3/$prefix_circularRNA_known_txt -b $align_bam -t -r $reference_annotation -o $quant
        echo "quant circ3 finish"

        #elimina los reads descomprimidos
            rm $read1 $read2
            rm *.fq
        #echo  "reads removed finish"   
            rm -rf $hisat_out_dir
        # remover el archivo de alineamientos muy pesado (32Gb)
        unset Aligned
            Aligned="Aligned.out.sam"
        unset prefix_Aligned
            prefix_Aligned=$prefix_name$Aligned
            rm $prefix_Aligned

        echo "Aligned.out.sam removed..."

        #circexplorer_results=/media/storage2/jpablo/Circexp3/$prefix_name*
        echo $circexplorer_results

        echo "creating circexplorer_results directory...."
        unset results
        results="_circexplorer_results"
        mkdir $prefix_name_sample$results
        mv $circexplorer_results /media/storage2/jpablo/Circexp3/$prefix_name_sample$results
        mv $back_spliced_junction /media/storage2/jpablo/Circexp3/$prefix_name_sample$results
        mv quant* /media/storage2/jpablo/Circexp3/$prefix_name_sample$results
        rm -r circexplorer_results
        #mv result_circ /media/storage2/jpablo/Circexp3/$prefix_name_sample$results
        echo "FINISH"

done


