process MERGE_FEATURECOUNTS {
	label 'process_inter'
	tag "${Sample}"
	publishDir "${PWD}/Final_Output/${Sample}/", mode: 'copy', pattern: '*merged.featureCounts.out'
	input:
		tuple val(Sample), path(bedfile), path(squid_bam)
		tuple val(Sample), file(featurecounts_out), file(featurecounts_summary)
		file(pon_myfu)
		file(pon_ball)
		file(pon_tall)
	output:
		tuple val(Sample), path(bedfile), file("${Sample}.merged.featureCounts.out")
	script:
	def bed_name = bedfile.getName()
	def pon_file = 
		bed_name == "myeloid_fusion02062022_hg38.bed" ? pon_myfu :
		bed_name == "RADICALv3_hg38_sortd.bed"       ? pon_ball :
		bed_name == "T-ALL02062022_hg38.bed"         ? pon_tall :
		null

	if (pon_file == null) {
		error "Unknown bedfile: ${bed_name}"
	}

	"""
	merge_featurecounts.py ${pon_file} ${featurecounts_out} ${Sample}.merged.featureCounts.out
	"""
}