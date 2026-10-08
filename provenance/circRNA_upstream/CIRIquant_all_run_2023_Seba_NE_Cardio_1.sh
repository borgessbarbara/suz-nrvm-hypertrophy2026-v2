export PATH=$PATH:/media/storage2/software/hisat2-2.2.1
for sample in {01..04}
do
    reads="/media/storage2/surquiza_docs/datasets_rat_HC/dataset_NE_cardios/lncRNA/raw_data_fastq/H${sample: -1}_ln*.fq.gz"

    zip_reads=$reads
    echo "Seq Files..."
    echo $zip_reads
    #Decompress the two files a la carpeta donde se va a ejecutar el programa
    echo "unzpping files"
    for f in $zip_reads
    do
      STEM=$(basename "${f}" .gz)
      gunzip -k -c "${f}" > /media/storage2/jpablo/CIRIquant/"${STEM}"
    done

    echo "Unziped files"

    #Create a variable to include the two unzip reads
    reads=/media/storage2/jpablo/CIRIquant/*.fq
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
    # trimmomatic PE -phred33 -threads 4 $read1 $read2 $read1.trimmed.fq $read1.un.trimmed.fq  $read2.trimmed.fq $read2.un.trimmed.fq ILLUMINACLIP:/media/storage/jpablo/circRNAs/CIRI2/z_adapters_sequencing.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 MINLEN:32 #CROP:140
    ./fastp -i $read1 -I $read2 -o $read1.trimmed.fq -O $read2.trimmed.fq --detect_adapter_for_pe -f 10 -F 10 -3 --cut_tail_mean_quality 28 -q 30 -w 8
    echo "TRIMMO END"

    reads=/media/storage2/jpablo/CIRIquant/*.fq
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
    
    echo "START CIRIquant" # Hisat mapping
     CIRIquant -t 30 -1 $read1 -2 $read2 --config /media/storage2/jpablo/CIRIquant/rn8.yml -o $prefix_name_sample -p $prefix_name_sample
     #CIRIquant -t 30 -1 $read1 -2 $read2 --config /media/storage2/jpablo/CIRIquant/rn8.yml -o $prefix_name_sample -p $prefix_name_sample --circ /media/storage2/jpablo/Circexp3/${prefix_name_sample}_circexplorer_results/${prefix_name_sample}_CIRCRNA_circularRNA_known.txt --tool CIRCexplorer2
    echo "FINISH"
    unset sorted
    sorted=".sorted.bam"
    denovosorted="_denovo.sorted.bam"
     rm *.fq
     rm $read1 $read2
     rm $prefix_name_sample/align/$prefix_name_sample$sorted
     rm $prefix_name_sample/circ/$prefix_name_sample$denovosorted

done
