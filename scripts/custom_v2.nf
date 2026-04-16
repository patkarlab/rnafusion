#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
"""

vardict = params.vardict
deepvariant = params.deepvariant
mutect2 = params.mutect2
genome_loc = file("${params.genome}", checkIfExists: true)
index_file = file("${params.genome}.fai", checkIfExists: true)
dict_file = file("${params.genome_dict}", checkIfExists: true)
known_SNPs = file("${params.dbsnp}", checkIfExists: true)
known_SNPs_index = file("${params.dbsnp_index}", checkIfExists: true)
sv_lib = file("${params.sv_lib}", checkIfExists: true)
sv_anno = file("${params.sv_anno}", checkIfExists: true)

include { VAR_RNA } from '../workflows/var_rna.nf'
include { ANNOVAR as ANNOVAR_VARDICT ; ANNOVAR as ANNOVAR_DEEPVARIANT ; ANNOVAR as ANNOVAR_MUTECT2 } from '../modules/local/annovar/annotate/main'
include { IGV_PREPROCESS as IGV_PREPROCESS_VARDICT ; IGV_PREPROCESS as IGV_PREPROCESS_MUTECT2 } from '../modules/local/annovar/igv_preprocess/main'
include { IGV_REPORTS as IGV_REPORTS_VARDICT ; IGV_REPORTS as IGV_REPORTS_MUTECT2 } from '../modules/local/igv_reports/main'
include { FORMAT_VARDICT } from '../modules/local/python/format_vardict/main'
include { FORMAT_MUTECT2 } from '../modules/local/python/format_mutect2/main'


process COUNTS {
	tag "${sampleId}"
	publishDir "${PWD}/Final_Output/${sampleId}/", mode: 'copy'
	input:
		tuple val(sampleId), path(bedfile), path(squid_bam)
	output:
		tuple val (sampleId), file("${sampleId}.counts_squid.bed")
	script:
	"""
	bedtools coverage -counts -a ${bedfile} -b ${squid_bam} > ${sampleId}.counts_squid.bed
	"""
}

process BAM {
	tag "${sampleId}"
	label 'process_inter'
	publishDir "${PWD}/Final_Output/${sampleId}/", mode: 'copy'
	input:
		tuple val(sampleId), path(bedfile), path(squid_bam)
	output:
		tuple val (sampleId), file ("${sampleId}.sorted.bam"), file ("${sampleId}.sorted.bam.bai")
	script:
	"""
	samtools sort -@ ${task.cpus} ${squid_bam} -o ${sampleId}.sorted.bam
	samtools index ${sampleId}.sorted.bam
	"""
}

process FILE_COPY {
	tag "${sampleId}"
	input:
		tuple val (sampleId), file(counts_squid), file(vardict_csv), file(haplotypecaller_csv), file(mutect2_csv)
	output:
		val (sampleId), emit:sample_id
	script:
	"""
	if [ -f ${PWD}/arriba/${sampleId}.arriba.fusions.tsv ]; then
		cp ${PWD}/arriba/${sampleId}.arriba.fusions.tsv ${PWD}/Final_Output/${sampleId}/
	fi

	if [ -f ${PWD}/arriba_visualisation/${sampleId}.pdf ]; then
		cp ${PWD}/arriba_visualisation/${sampleId}.pdf ${PWD}/Final_Output/${sampleId}/
	fi

	if [ -f ${PWD}/squid/${sampleId}.squid.fusions.annotated.txt ]; then
		cp ${PWD}/squid/${sampleId}.squid.fusions.annotated.txt ${PWD}/Final_Output/${sampleId}/
	fi

	if [ -f ${PWD}/pizzly/${sampleId}.pizzly.txt ]; then
		cp ${PWD}/pizzly/${sampleId}.pizzly.txt ${PWD}/Final_Output/${sampleId}/
	fi

	if [ -f ${PWD}/fusioncatcher/${sampleId}.fusioncatcher.fusion-genes.txt ]; then
		cp ${PWD}/fusioncatcher/${sampleId}.fusioncatcher.fusion-genes.txt ${PWD}/Final_Output/${sampleId}/
	fi

	if [ -f ${PWD}/fusioncatcher/${sampleId}.fusioncatcher.summary.txt ]; then
		sed.sh ${PWD}/fusioncatcher/${sampleId}.fusioncatcher.summary.txt
		cp ${PWD}/fusioncatcher/${sampleId}.fusioncatcher.summary.txt ${PWD}/Final_Output/${sampleId}/
	fi

	if [ -f ${PWD}/starfusion/${sampleId}.starfusion.fusion_predictions.tsv ]; then
		cp ${PWD}/starfusion/${sampleId}.starfusion.fusion_predictions.tsv ${PWD}/Final_Output/${sampleId}/
	fi

	if [ -d ${PWD}/fusionreport/${sampleId} ]; then
		cp -r ${PWD}/fusionreport/${sampleId} ${PWD}/Final_Output/${sampleId}/${sampleId}_fusionreport
	fi

	if [ -f ${PWD}/fusioninspector/${sampleId}.fusion_inspector_web.html ]; then 
		cp -r ${PWD}/fusioninspector/${sampleId}.fusion_inspector_web.html ${PWD}/Final_Output/${sampleId}/
	fi

	merge_variants.py ${sampleId}_variants.csv ${mutect2_csv} ${vardict_csv} ${haplotypecaller_csv}
	cp  ${sampleId}_variants.csv ${PWD}/Final_Output/${sampleId}/

	merge-csv_v3.py ${sampleId} ${PWD}/Final_Output/${sampleId}/${sampleId}.xlsx \
		${PWD}/Final_Output/${sampleId}/${sampleId}.counts_squid.bed \
		${PWD}/Final_Output/${sampleId}/${sampleId}.arriba.fusions.tsv \
		${PWD}/Final_Output/${sampleId}/${sampleId}.squid.fusions.annotated.txt \
		${PWD}/Final_Output/${sampleId}/${sampleId}.pizzly.txt \
		${PWD}/Final_Output/${sampleId}/${sampleId}.fusioncatcher.fusion-genes.txt \
		${PWD}/Final_Output/${sampleId}/${sampleId}.fusioncatcher.summary.txt \
		${PWD}/Final_Output/${sampleId}/${sampleId}.starfusion.fusion_predictions.tsv \
		${PWD}/Final_Output/${sampleId}/${sampleId}_variants.csv \
	"""
}

process FUSVIZ_ALL_SAMPLES {
	publishDir "${PWD}/Final_Output/", mode: 'copy', pattern: 'sv_output'
	input:
		val all_samples 
		path (sv_lib)
		path (sv_anno)
	output:
		path "sv_output"

	script:
	"""
	mkdir -p sv_input

	sample_ids="${all_samples.join(' ')}"

	for sampleId in \$sample_ids; do
		mkdir -p sv_input/\$sampleId

		if [ -f ${PWD}/arriba/\${sampleId}.arriba.fusions.tsv ]; then
			cp ${PWD}/arriba/\$sampleId.arriba.fusions.tsv sv_input/\${sampleId}/Arriba.tsv
		fi

		if [ -f ${PWD}/fusioncatcher/\${sampleId}.fusioncatcher.fusion-genes.txt ]; then
			cp ${PWD}/fusioncatcher/\$sampleId.fusioncatcher.fusion-genes.txt sv_input/\${sampleId}/Fusioncatcher.txt
		fi

		if [ -f ${PWD}/starfusion/\${sampleId}.starfusion.fusion_predictions.tsv ]; then
			cp ${PWD}/starfusion/\$sampleId.starfusion.fusion_predictions.tsv sv_input/\${sampleId}/STAR-fusion.tsv
		fi
	done

	export PERL5LIB="${sv_lib}"

	SV_standard.pl --genome hg38 --type RNA --anno ${sv_anno} --input sv_input --output sv_output
	"""
}

process CFF_FILEGEN {
	tag "${sampleId}"
	input:
		val(sampleId)
	output:
		tuple val (sampleId), file ("*.cff")	
	script:
	"""
	liftover.py ${params.chain_hg38_to_hg19} ${sampleId} ${PWD}/Final_Output/${sampleId}/${sampleId}.starfusion.fusion_predictions.tsv ${PWD}/Final_Output/${sampleId}/${sampleId}.fusioncatcher.fusion-genes.txt ${PWD}/Final_Output/${sampleId}/${sampleId}.squid.fusions.annotated.txt ${PWD}/Final_Output/${sampleId}/${sampleId}.arriba.fusions.tsv
	"""
}

process METAFUSION {
	tag "${sampleId}"
	input:
		tuple val(sampleId), file(cff_file)
	output:
		tuple val (sampleId), file("${sampleId}_final.n2.cluster.xlsx")
	script:
	"""
	if [ -s ${cff_file} ];then 
		mkdir ${sampleId}
		cp ${cff_file} ${sampleId}
		cd ${sampleId}

		metafus_gen.sh ${sampleId}.cff ${sampleId} > temp.sh

		# tool cutoff for docker
		tool_cutoff=\$(grep -i 'num_tools' temp.sh | sed 's:[^0-9]::g')
		num_tools=\$(awk 'BEGIN{FS="\\t"}{print \$11}' ${sampleId}.cff | uniq | sort | wc -l)

		# No. of tools in the input cff file
		if [ \${num_tools} -ge \${tool_cutoff} ]; then
			bash temp.sh || true
		fi
		
		if [ -f ${sampleId}/final.n2.cluster.xlsx ]; then
			cp ${sampleId}/final.n2.cluster.xlsx ../${sampleId}_final.n2.cluster.xlsx
		else
			touch ../${sampleId}_final.n2.cluster.xlsx
		fi
	else
		touch ${sampleId}_final.n2.cluster.xlsx
	fi
	"""
}

process FILTER_METAFUSION {
	tag "${sampleId}"
	input:
		tuple val(sampleId), file(final_cluster)
	output:
		tuple val (sampleId), file("${sampleId}_metafuse.xlsx")
	script:
	"""
	if [ -s ${final_cluster} ];then 
		filter_metafuse.py ${final_cluster} ${sampleId}_metafuse.xlsx
	else
		touch ${sampleId}_metafuse.xlsx
	fi
	"""
}

process UPDATE_METAFUSION_DB {
	tag "${sampleId}"
	input:
		tuple val(sampleId), file(metafusion_excel)
	output:
		tuple val (sampleId), file("${sampleId}_metafuse.xlsx")
	script:
	"""
	if [ -s ${metafusion_excel} ];then 
		metafus_table_append.sh ${metafusion_excel} > append_table.sh
		bash append_table.sh
	else
		echo "No metafusion output to append to database for ${sampleId}"
	fi
	"""
}

process LIFTOVER_METAFUSION {
	tag "${sampleId}"
	publishDir "${PWD}/Final_Output/${sampleId}/", mode: 'copy'
	input:
		tuple val(sampleId), file(metafusion_excel)
	output:
		tuple val (sampleId), file("${sampleId}_metafuse_hg38.xlsx")
	script:
	"""
	if [ -s ${metafusion_excel} ];then 
		convert_metafuse.py ${params.chain_hg19_to_hg38} ${sampleId}_metafuse.xlsx ${sampleId}_metafuse_hg38.xlsx
	else
		touch ${sampleId}_metafuse_hg38.xlsx
	fi
	"""	
}

process DASHBOARD {
	tag "${sampleId}"
	publishDir "${PWD}/Final_Output/${sampleId}/", mode: 'copy'
	input:
		tuple val(sampleId)
	output:
		tuple val (sampleId), file ("${sampleId}_dashboard.html")
	script:
	"""
	generate_Illumina_dashboard_v31.py --fusions ${PWD}/Final_Output/${sampleId}/${sampleId}.xlsx --cytoband ${params.cytoBand} --output ${sampleId}_dashboard.html
	"""
}

process VARDICT {
	tag "${sampleId}"
	label 'process_inter'
	input:
		tuple val(sampleId), path(bedfile), path(squid_bam), file (gatk_bam), file (gatk_bam_bai)
		path (GenFile)
		path (GenInd)
	output:
		tuple val(sampleId), file ("${sampleId}.vardict.vcf.gz")
	script:
	"""
	java -Xmx${task.memory.toGiga()}g -jar /usr/local/share/vardict-java-1.8.3-0/lib/VarDict-1.8.3.jar -G ${GenFile} -th ${task.cpus} -f 0.05 -N ${sampleId} -b ${gatk_bam} -c 1 -S 2 -E 3 -g 4 -L 10000000 --verbose ${bedfile} | teststrandbias.R | var2vcf_valid.pl | gzip > ${sampleId}.vardict.vcf.gz
	"""
}

process DEEPVARIANT {
	tag "${sampleId}"
	label 'process_inter'
	input:
		tuple val(sampleId), path(bedfile), path(squid_bam), file (gatk_bam), file (gatk_bam_bai)
		path (GenFile)
		path (GenInd)
	output:
		tuple val(sampleId), file ("${sampleId}.deepvar.vcf")
	script:
	"""
	/opt/deepvariant/bin/run_deepvariant --model_type=WES --ref=${GenFile} --regions=${bedfile} --reads=${gatk_bam} --output_vcf=${sampleId}.deepvar.vcf --num_shards=${task.cpus}
	"""
}

process MUTECT2 {
	tag "${sampleId}"
	label 'process_inter'
	input:
		tuple val(sampleId), path(bedfile), path(squid_bam), file (gatk_bam), file (gatk_bam_bai)
		path (GenFile)
		path (GenInd)
		path (GenDict)
		path (known_SNPs)
		path (known_SNPs_index)		
	output:
		tuple val (sampleId), file ("${sampleId}_mutect_passed.vcf")
	script:
	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}g" Mutect2 -R ${GenFile} -I:tumor ${gatk_bam} -O ${sampleId}_mutect.vcf --germline-resource ${known_SNPs} -L ${bedfile} --native-pair-hmm-threads ${task.cpus}
	gatk --java-options "-Xmx${task.memory.toGiga()}g" FilterMutectCalls -R ${GenFile} -V ${sampleId}_mutect.vcf -O ${sampleId}_mutect_filtered.vcf --stats ${sampleId}_mutect.vcf.stats
	gatk --java-options "-Xmx${task.memory.toGiga()}g" SelectVariants -R ${GenFile} -V ${sampleId}_mutect_filtered.vcf --exclude-filtered -O ${sampleId}_mutect_passed.vcf
	"""
}

process VEP {		
	tag "${sampleId}"
	label 'process_inter'
	publishDir "${PWD}/Final_Output/${sampleId}/", mode: 'copy'
	input:
		tuple val (sampleId), file(vardict_vcf), file(vardict_vcf_index)
	output:
		tuple val (sampleId), file ("${sampleId}_vardict_vep.txt")
	script:
	"""
	vep -i ${vardict_vcf} --fasta ${params.genome} --cache_version 113 -o ${sampleId}_vep.txt --offline --fork ${task.cpus} --tab --force_overwrite --symbol --protein --af --max_af  --no_check_alleles \
	--sift b --variant_class --canonical --allele_number --hgvs --shift_hgvs 1 --af_1kg --af_gnomadg --pubmed

	filter_vep -i ${sampleId}_vep.txt -o ${sampleId}_filtered.txt --filter "(CANONICAL is YES) and (AF < 0.01 or not AF)" --force_overwrite
	grep -v "##" ${sampleId}_filtered.txt > ${sampleId}_vep_delheaders.txt
	extract_vepdata.py ${sampleId}_vep_delheaders.txt ${sampleId}.extractedvepdelheaders.csv
	extract_vaf.py ${vardict_vcf} ${sampleId}.extracted.csv
	mergeDeepVariantVep.py ${sampleId}.extracted.csv ${sampleId}.extractedvepdelheaders.csv ${sampleId}_vardict_vep.txt
	"""
}


workflow COVERAGE {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.set { samples_ch }

	main:
	COUNTS(samples_ch)
	BAM(samples_ch)
	haplotypecaller_csv = VAR_RNA(samples_ch)
	VARDICT(samples_ch.join(haplotypecaller_csv.gatk_bam), genome_loc, index_file)
	//VEP(VARDICT.out)
	ANNOVAR_VARDICT(VARDICT.out, vardict)
	FORMAT_VARDICT(ANNOVAR_VARDICT.out)
	DEEPVARIANT(samples_ch.join(haplotypecaller_csv.gatk_bam), genome_loc, index_file)
	ANNOVAR_DEEPVARIANT(DEEPVARIANT.out, deepvariant)
	MUTECT2(samples_ch.join(haplotypecaller_csv.gatk_bam), genome_loc, index_file, dict_file, known_SNPs, known_SNPs_index)
	ANNOVAR_MUTECT2(MUTECT2.out, mutect2)
	FORMAT_MUTECT2(ANNOVAR_MUTECT2.out)
	FILE_COPY(COUNTS.out.join(FORMAT_VARDICT.out.join(haplotypecaller_csv.variants.join(FORMAT_MUTECT2.out))))
	FUSVIZ_ALL_SAMPLES(FILE_COPY.out.sample_id.collect(), sv_lib, sv_anno)
	CFF_FILEGEN(FILE_COPY.out.sample_id)
	METAFUSION(CFF_FILEGEN.out)
	FILTER_METAFUSION(METAFUSION.out)
	UPDATE_METAFUSION_DB(FILTER_METAFUSION.out)
	LIFTOVER_METAFUSION(FILTER_METAFUSION.out)
	DASHBOARD(FILE_COPY.out.sample_id)
	IGV_PREPROCESS_VARDICT(VARDICT.out, vardict)
	IGV_PREPROCESS_MUTECT2(MUTECT2.out, mutect2)
	IGV_REPORTS_VARDICT(haplotypecaller_csv.gatk_bam.join(IGV_PREPROCESS_VARDICT.out), genome_loc, index_file, vardict)
	IGV_REPORTS_MUTECT2(haplotypecaller_csv.gatk_bam.join(IGV_PREPROCESS_MUTECT2.out), genome_loc, index_file, mutect2)
}
workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}
