#Create a var8iable to include the two unzip reads
for m in {13..16}
do
    reads="/media/storage2/surquiza_docs/datasets_rat_HC/dataset_NE_cardios/lncRNA/raw_data_fastq/SQ${m}_ln_*.fq.gz"
    reference_fasta=/media/storage2/jpablo/Rattus_Norvegicus/rat_primary_assembly_orderr8.fa
    reference_annotation=/media/storage2/surquiza_docs/GTFs_files/rat/Rattus_norvegicus.GRCr8.114.gtf
    mkdir AS_results
    mkdir CIRI2_results
    mkdir sam_bwa_alignment
    echo $reads
    read1=0
    read2=0

        zip_reads=$reads
        echo "Seq Files..."
        echo $zip_reads
        #Decompress the two files a la carpeta donde se va a ejecutar el programa
        echo "unzpping files"
        for f in $zip_reads
        do
        STEM=$(basename "${f}" .gz)
        gunzip -k -c "${f}" > /media/storage2/jpablo/CIRI2/CIRI_v2.0.6/"${STEM}"
        done

        echo "Unziped files"

        #Create a variable to include the two unzip reads
        reads=/media/storage2/jpablo/CIRI2/CIRI_v2.0.6/*.fq
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
        prefix_name_sample=$(basename "${read1: -11:-8}")
        echo $prefix_name_sample
        
        unset prefix_name_sample_sam
        prefix_name_sample_sam=$(basename "${read1: -11:-8}".sam)
        
        unset prefix_name_sample_ciri
        prefix_name_sample_ciri=$(basename "${read1: -11:-8}".ciri)
        echo $prefix_name_sample_sam
        
        echo "Sample"
        echo $prefix_name_sample

        echo "cut read length to the same length"
        #trimmomatic PE -phred33 -threads 4 $read1 $read2 $read1.trimmed.fq $read1.un.trimmed.fq  $read2.trimmed.fq $read2.un.trimmed.fq ILLUMINACLIP:/media/storage1/jpablo/circRNAs/CIRI2/z_adapters_sequencing.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 MINLEN:140 CROP:140
        ./fastp -i $read1 -I $read2 -o $read1.trimmed.fq -O $read2.trimmed.fq --detect_adapter_for_pe -f 10 -F 10 -3 --cut_tail_mean_quality 28 -q 30 -w 8 -b 140 -B 140 -l 140
    
        echo "end cutadapt"
        reads=/media/storage2/jpablo/CIRI2/CIRI_v2.0.6/*.fq
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
        prefix_name_sample=$(basename "${read1: -22:-19}")
        echo $prefix_name_sample
        
        unset prefix_name_sample_sam
        prefix_name_sample_sam=$(basename "${read1: -22:-19}".sam)
        
        unset prefix_name_sample_ciri
        prefix_name_sample_ciri=$(basename "${read1: -22:-19}".ciri)
        echo $prefix_name_sample_sam
        echo "Sample"
        echo $prefix_name_sample

        echo "START bwa" 
        bwa mem -t 20 -T 19 $reference_fasta $read1 $read2 > $prefix_name_sample_sam
        echo "End bwa"

        echo "START CIRI2"
        perl CIRI2.pl -I $prefix_name_sample_sam -O $prefix_name_sample_ciri -F $reference_fasta -A $reference_annotation -T 20
        echo "END CIRI2"
        
        echo "START CIRI-AS"
        perl CIRI_AS_v1.2.pl -S $prefix_name_sample_sam -C $prefix_name_sample_ciri -F $reference_fasta -A $reference_annotation -T 20 -O  $prefix_name_sample.ASout -D yes
        echo "end CIRI-AS"

        echo "START CIRI-FULL"
        #java -jar /media/storage/jpablo/circRNAs/CIRI2/CIRI-full_v2.0/CIRI-full.jar Pipeline -1 $read1 -2 $read2 -d . -o out -r $reference_fasta -a $reference_annotation -t 20   
        echo "End CIRI-FULL"

        #java -jar /media/storage/jpablo/circRNAs/CIRI2/CIRI-full_v2.0/CIRI-full.jar RO1 -1 01_1_val_1.fq.trimmed.fq -2 01_2_val_2.fq.trimmed.fq 
        #java -jar /media/storage/jpablo/circRNAs/CIRI2/CIRI-full_v2.0/CIRI-full.jar RO2 -r $reference_fasta -s /media/storage/jpablo/circRNAs/CIRI2/CIRI_v2.0.6/sam_bwa_alignment/01.sam -l 140 -o 01
        #java -jar /media/storage/jpablo/circRNAs/CIRI2/CIRI-full_v2.0/CIRI-full.jar Merge -a $reference_annotation -c /media/storage/jpablo/circRNAs/CIRI2/CIRI_v2.0.6/CIRI2_results/01.ciri -as /media/storage/jpablo/circRNAs/CIRI2/CIRI_v2.0.6/AS_results/01.ASout_jav.list -ro 01_ro2_info.list -o out_merge_ -r $reference_fasta

        rm *.fq
        mv *ASout* AS_results/
        mv *.ciri CIRI2_results/
        mv *.sam* sam_bwa_alignment/

        echo "START cirivis"
       #for l in {01..08} 
       #do 
      #  java -jar /media/storage1/jpablo/circRNAs/CIRI2/CIRI-vis_v1.4.jar -i /media/storage1/jpablo/circRNAs/CIRI2/CIRI_v2.0.6/AS_results/${l}.ASout_jav.list -l /media/storage1/jpablo/circRNAs/CIRI2/CIRI_v2.0.6/AS_results/${l}.ASout_library_length.list -r /media/storage/jpablo/Rattus_Norvegicus/BWA_index_rn6/rn6.fa -list select_bsj.list -d ${l}_vis -o ${l}_vis
       #done
        #echo "END ciri vis"

        echo "FINISH"
done
