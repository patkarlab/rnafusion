process FEATURECOUNTS {
    label 'process_inter'
	tag "${Sample}"
	input:
		tuple val (Sample), file(bam), file(bamBai)
		file (ensembl_gtf)
	output:
		tuple val (Sample), file("${Sample}.featureCounts.out"), file("${Sample}.featureCounts.out.summary")
	script:
	"""
	featureCounts -p -a ${ensembl_gtf} -o ${Sample}.featureCounts.out -T ${task.cpus} ${bam}
	"""	
}
