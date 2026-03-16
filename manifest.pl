#! /usr/bin/perl -w
# This script will output the pipeline manifest for npm1-mrd pipeline on the Amazon server no. 2

$sample_file = $ARGV[0];
chomp ($sample_file);

opendir(DIR, "./samples");
@files = grep(/\.fastq.gz$/,readdir(DIR));
closedir(DIR);

print "run_number,barcode_id,project,assay,sample_type,forward_read,reverse_read\n";
open (INFILE, $sample_file) || die ("cannot open $sample_file\n");
while ($samp_id=<INFILE>) {
	chomp ($samp_id);
	@file_name = grep (/$samp_id/,@files);
	@R1 = grep (/R1/,@file_name);
	@S = split (/_/, $R1[0]);
	@R2 = grep (/R2/,@file_name);	
	#print "220419,$samp_id,220419,NPM1_MRD,$S[1],/data/home2/hematopath/hemoseq_v2/npm1_mrd/samples/@R1,/data/home2/hematopath/hemoseq_v2/npm1_mrd/samples/@R2\n";
	print "$samp_id,/home/diagnostics/pipelines/nf-core/rnafusion/samples/@R1,/home/diagnostics/pipelines/nf-core/rnafusion/samples/@R2,forward\n";
	}
close (INFILE);
