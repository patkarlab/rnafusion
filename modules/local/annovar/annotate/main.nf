process ANNOVAR {
	tag "${Sample}"
	label 'process_inter'
	input:
		tuple val (Sample), path(Vcf), path(VcfTbi)
		val(variant_caller)
	output:
		tuple val (Sample), path("${Sample}_${variant_caller}.out.hg38_multianno.csv")
	script:
	"""
	convert2annovar.pl -format vcf4 ${Vcf} --outfile ${Sample}_${variant_caller}.avinput --withzyg --includeinfo
	table_annovar.pl ${Sample}_${variant_caller}.avinput --out ${Sample}_${variant_caller}.out --remove --protocol refGene,cytoBand,avsnp151,intervar_20180118,1000g2015aug_all,cosmic70,clinvar_20250721,gnomad211_exome --operation g,r,f,f,f,f,f,f --buildver hg38 --nastring '-1' --otherinfo --csvout \
	--thread ${task.cpus} /databases/humandb -xreffile /databases/gene_fullxref.txt

	touch ${Sample}_${variant_caller}.out.hg38_multianno.csv
	"""
}
