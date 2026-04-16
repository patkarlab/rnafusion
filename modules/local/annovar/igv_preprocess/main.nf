process IGV_PREPROCESS {
	tag "${Sample}"
	label 'process_inter'
	input:
		tuple val (Sample), path(Vcf)
		val(variant_caller)
	output:
		tuple val (Sample), path("${Sample}_${variant_caller}.annovar.hg38_multianno.vcf")
	script:
	"""
	table_annovar.pl ${Vcf} --out ${Sample}_${variant_caller}.annovar --remove --protocol refGene,cytoBand,avsnp151,intervar_20180118,1000g2015aug_all,cosmic70,clinvar_20250721,gnomad211_exome \
	--operation g,r,f,f,f,f,f,f --buildver hg38 --nastring . --otherinfo --thread ${task.cpus} /databases/humandb -xreffile /databases/gene_fullxref.txt -vcfinput 
	"""
}
