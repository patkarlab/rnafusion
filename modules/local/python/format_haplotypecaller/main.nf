process FORMAT_HAPLOTYPECALLER {
	tag "${Sample}"
	label 'process_inter'
	input:
		tuple val (Sample), path(multianno_csv)
	output:
		tuple val (Sample), file("${Sample}_varRNA.csv")
	script:
	"""
	${params.format_haplotypecaller} ${multianno_csv} ${Sample}_varRNA.csv
	"""
}