process FORMAT_VARDICT {
	tag "${Sample}"
	label 'process_inter'
	input:
		tuple val (Sample), path(multianno_csv)
	output:
		tuple val (Sample), file("${Sample}_vardict.csv")
	script:
	"""
	${params.format_vardict} ${multianno_csv} ${Sample}_vardict.csv
	"""
}