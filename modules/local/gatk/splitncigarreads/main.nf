process SPLIT_N_CIGARREADS {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), file(markdups_bam), file(markdups_metrics)
		path (GenFile)
		path (GenInd)
		path (GenDict)
	output:
		tuple val(Sample), file("${Sample}.splitncigar.bam")
	script:
	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}g" SplitNCigarReads -R ${GenFile} -I ${markdups_bam} -O ${Sample}.splitncigar.bam
	"""
}