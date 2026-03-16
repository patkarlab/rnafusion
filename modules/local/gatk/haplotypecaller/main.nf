process HAPLOTYPECALLER {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val(Sample), file(final_bam), file(final_bam_bai), path(bedfile), path(squid_bam)
		path (GenFile)
		path (GenInd)
		path (GenDict)
	output:
		tuple val(Sample), file("${Sample}.vcf.gz"), file("${Sample}.vcf.gz.tbi")
	script:
	"""
	gatk HaplotypeCaller -I ${final_bam} -R ${GenFile} -L ${bedfile} -O ${Sample}.vcf.gz \
		-dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20 --max-reads-per-alignment-start 0 \
		-A BaseQuality -A BaseQualityRankSumTest -A FragmentLength -A MappingQualityRankSumTest \
		-A MappingQuality -A AlleleFraction -A CountNs -A Coverage -A ReadPosition -A ReadPosRankSumTest
	"""
}
