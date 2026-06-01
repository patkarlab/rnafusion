process HEATMAP_GEN {
	label 'process_inter'
	tag "${Sample}"
	publishDir "${PWD}/Final_Output/${Sample}/", mode: 'copy'
	input:
		tuple val(Sample), path(bedfile), file(merged_counts)
		file (gene_map)
		file (myfu_genes)
		file (ball_genes)
		file (tall_genes)        
	output:
		tuple val (Sample), file("${Sample}.merged.featureCounts_filtered_TPM_matrix.tsv"), file("${Sample}.merged.featureCounts_TPM_heatmap.png")
	script:
	def bed_name = bedfile.getName()
	def gene_file =
		bed_name == "myeloid_fusion02062022_hg38.bed" ? myfu_genes :
		bed_name == "RADICALv3_hg38_sortd.bed"        ? ball_genes :
		bed_name == "T-ALL02062022_hg38.bed"          ? tall_genes :
		null

	if (gene_file == null) {
		error "Unknown bedfile: ${bed_name}"
	}

	"""
	generate_heatmap.py ${merged_counts} ${gene_map} ${gene_file}
	"""
}