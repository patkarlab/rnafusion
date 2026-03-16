process BQSR {
	tag "${Sample}"
	label 'process_inter'
	input:
		tuple val(Sample), file(splitncigar_bam)
		path (GenFile)
		path (GenInd)
		path (GenDict)
		path (knownSNPS)
		path (knownSNPS_index)
	output:
		tuple val(Sample), file("${Sample}_recal.table")
	script:
	"""
	gatk BaseRecalibrator \
		-I ${splitncigar_bam} \
		-R ${GenFile} \
		--known-sites ${knownSNPS} \
		-O ${Sample}_recal.table
	"""
}