process ADD_READGROUPS {
	tag "${Sample}"
	label 'process_inter'
	input:
		tuple val(Sample), path(bedfile), path(squid_bam)
	output:
		tuple val(Sample), file("${Sample}.rg.sortd.bam")
	script:
	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}g" AddOrReplaceReadGroups --I ${squid_bam} --O ${Sample}.rg.sortd.bam --RGID 1 --RGLB lib1 --RGPL ILLUMINA --RGPU unit1 --RGSM ${Sample} --SORT_ORDER coordinate
	"""
}