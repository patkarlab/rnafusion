process MARK_DUPLICATES {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), file(rg_bam)
	output:
		tuple val(Sample), file("${Sample}.marked_duplicates.bam"), file("${Sample}.metrics.txt")
	script:
	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}g" MarkDuplicates -I ${rg_bam} -O ${Sample}.marked_duplicates.bam -M ${Sample}.metrics.txt --CREATE_INDEX true
	"""
}