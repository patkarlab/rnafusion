
genome_loc = file("${params.genome}", checkIfExists: true)
index_file = file("${params.genome}.fai", checkIfExists: true)
dict_file = file("${params.genome_dict}", checkIfExists: true)
known_SNPs = file("${params.dbsnp}", checkIfExists: true)
known_SNPs_index = file("${params.dbsnp_index}", checkIfExists: true)
haplotypecaller = params.haplotypecaller


include { ADD_READGROUPS } from '../modules/local/gatk/addreadgroups/main'
include { MARK_DUPLICATES } from '../modules/local/gatk/markduplicates/main'
include { SPLIT_N_CIGARREADS } from '../modules/local/gatk/splitncigarreads/main'
include { BQSR } from '../modules/local/gatk/baserecalibrator/main'
include { APPLY_BQSR } from '../modules/local/gatk/applybqsr/main'
include { HAPLOTYPECALLER } from '../modules/local/gatk/haplotypecaller/main'
include { ANNOVAR as ANNOVAR_HAPLOTYPECALLER } from '../modules/local/annovar/annotate/main'
include { FORMAT_HAPLOTYPECALLER } from '../modules/local/python/format_haplotypecaller/main'

workflow VAR_RNA {
	take:
		samples_ch
	main:
	ADD_READGROUPS(samples_ch)
	MARK_DUPLICATES(ADD_READGROUPS.out)
	SPLIT_N_CIGARREADS(MARK_DUPLICATES.out, genome_loc, index_file, dict_file)
	BQSR(SPLIT_N_CIGARREADS.out, genome_loc, index_file, dict_file, known_SNPs, known_SNPs_index)
	APPLY_BQSR(SPLIT_N_CIGARREADS.out.join(BQSR.out), genome_loc, index_file, dict_file)
	HAPLOTYPECALLER(APPLY_BQSR.out.join(samples_ch), genome_loc, index_file, dict_file)
	ANNOVAR_HAPLOTYPECALLER(HAPLOTYPECALLER.out, haplotypecaller)
	FORMAT_HAPLOTYPECALLER(ANNOVAR_HAPLOTYPECALLER.out)
	emit:
		FORMAT_HAPLOTYPECALLER.out
}

