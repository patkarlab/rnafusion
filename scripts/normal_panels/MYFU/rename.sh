while IFS=, read -r sample normal
do
    mv "${sample}.sorted.bam" "${normal}.bam"
    mv "${sample}.sorted.bam.bai" "${normal}.bam.bai"
done < samplesheet.csv
