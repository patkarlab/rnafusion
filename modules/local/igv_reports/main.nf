process IGV_REPORTS {
	tag "${Sample}"
	label 'process_inter'
	publishDir "${PWD}/Final_Output/${Sample}/", mode: 'copy', pattern: '*.html'
	input:
		tuple val(Sample), file (bam), file (bamBai), file (MultiannoVcf)
		path (GenFile)
		path (GenInd)
		val(variant_caller)
	output:
		tuple val(Sample), file("${Sample}_${variant_caller}_igv.html")
	script:
	"""

	if [ \$(grep -vcE '^#|^[[:space:]]*\$' ${MultiannoVcf}) -eq 0 ]; then
		echo "No variants found" > ${Sample}_${variant_caller}_igv.html
	else
		create_report ${MultiannoVcf} --fasta ${GenFile} --standalone --tracks ${bam} --output ${Sample}_${variant_caller}_igv.html
	fi
	"""		
}