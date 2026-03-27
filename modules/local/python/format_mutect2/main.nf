process FORMAT_MUTECT2 {
	tag "${Sample}"
	label 'process_inter'
	input:
		tuple val (Sample), path(multianno_csv)
	output:
		tuple val (Sample), file("${Sample}_mutect2.csv")
	script:
	"""
	format_mutect2.py ${multianno_csv} ${Sample}_mutect2.csv
	"""
}
