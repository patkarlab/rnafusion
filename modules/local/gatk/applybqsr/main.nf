process APPLY_BQSR {
	tag "${Sample}"
	label 'process_inter'
	input:
		tuple val(Sample), file(splitncigar_bam), file(recal_table)
		path (GenFile)
		path (GenInd)
		path (GenDict)
	output:
		tuple val(Sample), file("${Sample}_final.bam"), file("${Sample}_final.bam.bai")
	script:
	"""
	gatk ApplyBQSR \
		-R ${GenFile} \
		-I ${splitncigar_bam} \
		--bqsr-recal-file ${recal_table} \
		-O ${Sample}_final.bam

	mv ${Sample}_final.bai ${Sample}_final.bam.bai
	"""
}