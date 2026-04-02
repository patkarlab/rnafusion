#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX;
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use File::Spec;
use List::Util qw< min max sum >;
use Anno;
use RNAcaller;
use DNAcaller;

	#-----------------------------------------
	# Set options in the script command line
	#-----------------------------------------
	my @usage;
	push @usage, "Usage: ".basename($0)." [options]\n";
	push @usage, "  --help\n";
	push @usage, "  --type {DNA, RNA}		Set input SV calls from DNA-seq or RNA-seq\n";
	push @usage, "  --genome {hg19, hg38, GRCm39}	Set genome version used for SV calls\n";
	push @usage, "  --filter {PASS, ALL}		Set an option for filtering raw SV calls before merging (default: PASS), only available for DNA SVs\n";
	push @usage, "  --support {min, max, median}	Set a method to obtain split and spanning read support if SVs from multiple callers are available (default: median)\n";
	push @usage, "  --offset			Set an offset value for extending a gene interval, e.g. [start-offset, end+offset], default:1000\n";
	push @usage, "  --anno			Set annotation file directory\n";
	push @usage, "  --input			Set input directory\n";
	push @usage, "  --output			Set output directory\n";

	my $help;
	my $type;
	my $genome;
	my $filter;
	my $support;
	my $offset;
	my $anno;
	my $input;
	my $output;

	GetOptions
	(
		'help'		=> \$help,
		'type=s'	=> \$type,
		'genome=s'	=> \$genome,
		'filter=s'	=> \$filter,
		'support=s'	=> \$support,
		'offset=i'	=> \$offset,
		'anno=s'	=> \$anno,
		'input=s'	=> \$input,
		'output=s'	=> \$output
	);
	not defined $help or die @usage;

	defined $input or die @usage;
	defined $output or die @usage;
	defined $type or die strftime("%Y-%m-%d %H:%M:%S", localtime), " [option control] --type is a required option.\n";
	defined $genome or die strftime("%Y-%m-%d %H:%M:%S", localtime), " [option control] --genome is a required opton.\n";
	defined $anno or die " [option control] --anno is a required opton.\n";
	if (! defined($support) ) { $support = 'median'; }
	if (! defined($offset) ) { $offset = 1000; }
	if (! -e "$anno/gene_interval_hg19.bed.gz") { die strftime("%Y-%m-%d %H:%M:%S", localtime), " [option control] Gene interval annotation for hg19 not exists, quit!\n"; }
	if (! -e "$anno/gene_interval_hg38.bed.gz") { die strftime("%Y-%m-%d %H:%M:%S", localtime), " [option control] Gene interval annotation for hg38 not exists, quit!\n"; }
	if (! -e "$anno/gene_interval_GRCm39.bed.gz") { die strftime("%Y-%m-%d %H:%M:%S", localtime), " [option control] Gene interval annotation for GRCm39 not exists, quit!\n"; }
	if (! -e "$anno/exon_interval_hg19.bed.gz") { die strftime("%Y-%m-%d %H:%M:%S", localtime), " [option control] Exon interval annotation for hg19 not exists, quit!\n"; }
	if (! -e "$anno/exon_interval_hg38.bed.gz") { die strftime("%Y-%m-%d %H:%M:%S", localtime), " [option control] Exon interval annotation for hg38 not exists, quit!\n"; }
	if (! -e "$anno/exon_interval_GRCm39.bed.gz") { die strftime("%Y-%m-%d %H:%M:%S", localtime), " [option control] Exon interval annotation for GRCm39 not exists, quit!\n"; }
	if (! -e "$anno/human_ensembl_symbol_1_1.txt") { die strftime("%Y-%m-%d %H:%M:%S", localtime), " [option control] Human Ensembl_id and Gene_symbol info not exists, quit!\n"; }
	if (! -e "$anno/mouse_ensembl_symbol_1_1.txt") { die strftime("%Y-%m-%d %H:%M:%S", localtime), " [option control] Mouse Ensembl_id and Gene_symbol info not exists, quit!\n"; }

	if (! defined($filter) ) { $filter = "PASS"; } else {
		if ( $filter eq 'PASS'  || $filter eq 'ALL' ) {} else { die strftime("%Y-%m-%d %H:%M:%S", localtime), " [option control] --filter MUST be 'PASS' or 'ALL'.\n"; } }
		
	# check output folder status
	if ( -e $output ) { # if output folder present, please clean up firstly.
		`rm -rf $output/*`;
	} else {
		`mkdir $output`; 
	}

	#------------------------------------------------------
	# load annotation files (either $type 'DNA' or 'RNA')
	#------------------------------------------------------
	my %gene_interval_dna; # - gene interval for DNA SVs
	my %gene_interval_rna; # - gene interval for RNA SVs
	my %exon_interval_rna; # - exon interval per gene for RNA SVs
	my %ens_symbol_rna; # - ensembl id => gene symbol
	if ( $type eq 'DNA' ) {
		if ( $genome eq 'hg19' ) {
			Anno::process_gene_dna($anno, $genome, \%gene_interval_dna, $offset);
		} elsif ( $genome eq 'hg38' ) {
			Anno::process_gene_dna($anno, $genome, \%gene_interval_dna, $offset);
		} elsif ( $genome eq 'GRCm39' ) {
			Anno::process_gene_dna($anno, $genome, \%gene_interval_dna, $offset);
		} else {
			die strftime("%Y-%m-%d %H:%M:%S", localtime), " [load annotation] Genome version not correct for DNA SVs!\n";
		}
	} elsif ( $type eq 'RNA' ) {
		if ( $genome eq 'hg19' ) {
			Anno::process_gene_rna($anno, $genome, \%gene_interval_rna, \%exon_interval_rna, \%ens_symbol_rna, $offset);
		} elsif ( $genome eq 'hg38' ) {
			Anno::process_gene_rna($anno, $genome, \%gene_interval_rna, \%exon_interval_rna, \%ens_symbol_rna, $offset);
		} elsif ( $genome eq 'GRCm39' ) {
			Anno::process_gene_rna($anno, $genome, \%gene_interval_rna, \%exon_interval_rna, \%ens_symbol_rna, $offset);
		} else {
			die strftime("%Y-%m-%d %H:%M:%S", localtime), " [load annotation] Genome version not correct for RNA SVs!\n";
		}
	} else {
		die strftime("%Y-%m-%d %H:%M:%S", localtime), " [load annotation] Sequencing type not correct!\n";
	}

	
	if ( $type eq 'RNA' ) {
		#-------------------------------------------------------------------------------
		# Load and merge RNA SVs from csv file (can be one caller or multiple callers
		#-------------------------------------------------------------------------------
		# // [merge RNA SVs] - step.1: check input file format status, then load them
		my %input_hash_rna;
		my %count_hash_rna; # if both geneA-geneB and geneB-geneA are present, count their frequency, respectively.
		open(DIR, "ls $input |") || die strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge RNA SVs] - step.1: cannot open directory for RNA SVs input:$!\n";
		while ( <DIR> ) {
			chomp $_; my $sample = $_;
			if ( -d "$input/$sample" ) {
				print "#####################\n#$sample\n####################\n";
				open(FL, "ls $input/$sample/* |") || die strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge RNA SVs] - step.1: cannot open input file:$!\n"; # txt or tsv
				while ( <FL>) {
					chomp $_; my $file = $_;
					if ( $file =~/Arriba\.(tsv|txt)$/ ) {
						my $text = `cat $file | cut -s -f3,4,5,6,10,11,12,21,22`; # select relevant columns of Arriba output
						my %sss; $input_hash_rna{$sample}{'Arriba'} = \%sss; # define a reference to a hash structure
						RNAcaller::Arriba_support($sample, $input_hash_rna{$sample}{'Arriba'}, $text, $offset, \%gene_interval_rna, \%exon_interval_rna, $output);
					} elsif ( $file =~/STAR\-fusion\.(tsv|txt)$/ ) {
						my $text = `cat $file | cut -s -f2,3,7,8,9,10`; # select relevant columns of STAR-fusion output
						my %sss; $input_hash_rna{$sample}{'STAR-fusion'} = \%sss; # define a reference to a hash structure
						RNAcaller::STAR_fusion_support($sample, $input_hash_rna{$sample}{'STAR-fusion'}, $text, $offset, \%gene_interval_rna, \%exon_interval_rna, $output);
					} elsif ( $file =~/Fusioncatcher\.(tsv|txt)$/ ) {
						my $text = `cat $file | cut -s -f5,6,9,10,11,12`; # select relevant columns of Fusioncatcher output
						my %sss; $input_hash_rna{$sample}{'Fusioncatcher'} = \%sss; # define a reference to a hash structure
						RNAcaller::Fusioncatcher_support($sample, $input_hash_rna{$sample}{'Fusioncatcher'}, $text, $offset, \%gene_interval_rna, \%exon_interval_rna, $output);
					} elsif ( $file =~/deFuse\.(tsv|txt)$/ ) {
						my $text = `cat $file | cut -s -f3,4,5,6,12,21,22,25,26,38,39,40,41,51,57,65`; # select relevant columns of deFuse output
						my %sss; $input_hash_rna{$sample}{'deFuse'} = \%sss; # define a reference to a hash structure
						RNAcaller::deFuse_support($sample, $input_hash_rna{$sample}{'deFuse'}, $text, $offset, \%gene_interval_rna, \%exon_interval_rna, $output);
					} elsif ( $file =~/Dragen\.(tsv|txt)$/ ) {
						my $text = `cat $file | cut -s -f3,4,9,10,11,13`; # select relevant columns of deFuse output
						my %sss; $input_hash_rna{$sample}{'Dragen'} = \%sss; # define a reference to a hash structure
						RNAcaller::Dragen_support($sample, $input_hash_rna{$sample}{'Dragen'}, $text, $offset, \%gene_interval_rna, \%exon_interval_rna, $output);
					} else {
						print strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge RNA SVs] - step.1: Unknown $file, the caller is not in acceptable list!\n";
					}
				}
				close FL;

				# merge RNA SVs from multiple callers per sample
				my @caller_num = sort {$a cmp $b} keys %{$input_hash_rna{$sample}};
				my %sss; $input_hash_rna{$sample}{'merge'} = \%sss;
				for (my $i=0; $i < scalar(@caller_num); $i++ ) {
					foreach my $break1 ( sort {$a <=> $b} keys %{$input_hash_rna{$sample}{$caller_num[$i]}} ) {
						foreach my $break2 ( sort {$a <=> $b} keys %{$input_hash_rna{$sample}{$caller_num[$i]}{$break1}} ) {
							my $chrA = $input_hash_rna{$sample}{$caller_num[$i]}{$break1}{$break2}[0];
							my $chrB = $input_hash_rna{$sample}{$caller_num[$i]}{$break1}{$break2}[1];
							my $strandA = $input_hash_rna{$sample}{$caller_num[$i]}{$break1}{$break2}[2];
							my $strandB = $input_hash_rna{$sample}{$caller_num[$i]}{$break1}{$break2}[3];
							my $geneA = $input_hash_rna{$sample}{$caller_num[$i]}{$break1}{$break2}[4];
							my $geneB = $input_hash_rna{$sample}{$caller_num[$i]}{$break1}{$break2}[5];
							my $split = $input_hash_rna{$sample}{$caller_num[$i]}{$break1}{$break2}[6];
							my $span = $input_hash_rna{$sample}{$caller_num[$i]}{$break1}{$break2}[7];
							if ( exists($input_hash_rna{$sample}{'merge'}{$geneA}{$geneB}) ) {
								$count_hash_rna{$geneA}{$geneB}++;
								push @{$input_hash_rna{$sample}{'merge'}{$geneA}{$geneB}{$break1}{$break2}}, [$chrA, $chrB, $strandA, $strandB, $split, $span];
							} elsif ( exists($input_hash_rna{$sample}{'merge'}{$geneB}{$geneA}) ) {
								$count_hash_rna{$geneB}{$geneA}++;
								push @{$input_hash_rna{$sample}{'merge'}{$geneB}{$geneA}{$break2}{$break1}}, [$chrB, $chrA, $strandB, $strandA, $split, $span];
							} else {
								$count_hash_rna{$geneA}{$geneB}++;
								push @{$input_hash_rna{$sample}{'merge'}{$geneA}{$geneB}{$break1}{$break2}}, [$chrA, $chrB, $strandA, $strandB, $split, $span];
							}
						}
					}
				}
			} else {
				print strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge RNA SVs] - step.1: $sample is not a available folder\n";
			}
		}
		close DIR;

		# // [merge RNA SVs] - step.2: final print merging results
		open (OUT, ">$output/Final_RNA_SVs.txt") || die "final print merging results of RNA SVs:$!\n";
		print OUT "chrom1\tpos1\tgene1\tchrom2\tpos2\tgene2\tname\tsplit\tspan\tstrand1\tstrand2\tuntemplated_insert\tcomment\n";
		foreach my $sample ( keys %input_hash_rna ) {
			foreach my $geneA ( keys %{$input_hash_rna{$sample}{'merge'}} ) {
				foreach my $geneB ( keys %{$input_hash_rna{$sample}{'merge'}{$geneA}} ) {
					foreach my $breakpoint1 ( keys %{$input_hash_rna{$sample}{'merge'}{$geneA}{$geneB}} ) {
						foreach my $breakpoint2 ( keys %{$input_hash_rna{$sample}{'merge'}{$geneA}{$geneB}{$breakpoint1}} ) {
							my $chrA = $input_hash_rna{$sample}{'merge'}{$geneA}{$geneB}{$breakpoint1}{$breakpoint2}[0][0];
							my $chrB = $input_hash_rna{$sample}{'merge'}{$geneA}{$geneB}{$breakpoint1}{$breakpoint2}[0][1];
							my $strandA = $input_hash_rna{$sample}{'merge'}{$geneA}{$geneB}{$breakpoint1}{$breakpoint2}[0][2];
							my $strandB = $input_hash_rna{$sample}{'merge'}{$geneA}{$geneB}{$breakpoint1}{$breakpoint2}[0][3];

							my @reads = @{$input_hash_rna{$sample}{'merge'}{$geneA}{$geneB}{$breakpoint1}{$breakpoint2}};
							my @split;	my @span;	my $split_value;	my $span_value;
							for ( my $i=0; $i < scalar(@reads); $i++ ) { push @split, $reads[$i]->[4]; push @span, $reads[$i]->[5]; }
							@split = sort @split;	@span = sort @span;
							if ( $support eq 'min' ) {
								$split_value = $split[0];	$span_value = $span[0];
							} elsif ( $support eq 'max' ) {
								$split_value = $split[-1];	$span_value = $span[-1];
							} elsif ( $support eq 'median' ) {
								my $index = int(scalar(@split)/2);
								if ( @split % 2 ) {
									$split_value = $split[$index];	$span_value = $span[$index];
								} else {
									$split_value = int(($split[$index-1] + $split[$index])/2);	$span_value = int(($span[$index-1] + $span[$index])/2);
								}
							} else {
								die "--support option should be one of [min, max, median]\n";
							}

							my $symbolA; if ( exists($ens_symbol_rna{$geneA}) ) { if ( $ens_symbol_rna{$geneA}[0] eq '' ) { $symbolA = $geneA; } else { $symbolA = $ens_symbol_rna{$geneA}[0]; } } else { $symbolA = $geneA; }
							my $symbolB; if ( exists($ens_symbol_rna{$geneB}) ) { if ( $ens_symbol_rna{$geneB}[0] eq '' ) { $symbolB = $geneB; } else { $symbolB = $ens_symbol_rna{$geneB}[0]; } } else { $symbolB = $geneB; }
							if ( exists($count_hash_rna{$geneB}{$geneA}) ) {
								if ( $count_hash_rna{$geneA}{$geneB} > $count_hash_rna{$geneB}{$geneA} ) { # geneA-geneB is more frequent than geneB-geneA
									print OUT "$chrA\t$breakpoint1\t$symbolA\t$chrB\t$breakpoint2\t$symbolB\t$sample\t$split_value\t$span_value\t$strandA\t$strandB\t\t\n";
								} elsif ( $count_hash_rna{$geneA}{$geneB} < $count_hash_rna{$geneB}{$geneA} ) { # geneB-geneA is more frequent than geneA-geneB
									print OUT "$chrB\t$breakpoint2\t$symbolB\t$chrA\t$breakpoint1\t$symbolA\t$sample\t$split_value\t$span_value\t$strandB\t$strandA\t\t\n";
								} else {
									if ( $geneA lt $geneB ) {
										print OUT "$chrA\t$breakpoint1\t$symbolA\t$chrB\t$breakpoint2\t$symbolB\t$sample\t$split_value\t$span_value\t$strandA\t$strandB\t\t\n";
									} else {
										print OUT "$chrB\t$breakpoint2\t$symbolB\t$chrA\t$breakpoint1\t$symbolA\t$sample\t$split_value\t$span_value\t$strandB\t$strandA\t\t\n";
									}
								}
							} else {
								print OUT "$chrA\t$breakpoint1\t$symbolA\t$chrB\t$breakpoint2\t$symbolB\t$sample\t$split_value\t$span_value\t$strandA\t$strandB\t\t\n";
							}
						}
					}
				}
			}
		}
		close OUT;
	} else {
		#--------------------------------------------------------------------------------
		# Load and merge DNA SVs from VCF files (can be one caller or multiple callers)
		#--------------------------------------------------------------------------------
		# // [merge DNA SVs] - step.1: check input file format status (either 'vcf' or 'vcf.gz' is acceptable), then load them
		my %input_hash_dna; 
		open(DIR, "ls $input |") || die strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge DNA SVs] - step.1: cannot open directory for DNA SVs input:$!\n";
		while ( <DIR> ) {
			chomp $_; my $sample = $_;
			if ( -d "$input/$sample" ) {
				open(FL, "ls $input/$sample/*.vcf* |") || die strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge DNA SVs] - step.1: cannot open input file (vcf or vcf.gz):$!\n"; # either vcf or vcf.gz
				while ( <FL> ) {
					chomp $_; my $file = $_;
					my $text; # get vcf content
					if ( $file =~/\.vcf$/ ) {
						if ( $filter eq 'PASS' ) {
 							$text = `grep -P '^#|\tPASS\t' $file`; # only select variants tagged by 'PASS'
						} else {
							$text = `cat $file`; # select all variants
						}
					} elsif ( $file =~/\.vcf\.gz$/ ) {
						if ( $filter eq 'PASS' ) {
 							$text = `zgrep -P '^#|\tPASS\t' $file`; # only select variant tagged by 'PASS'
						} else {
							$text = `zcat $file`; # select all variants
						}
					} else {
						print strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge DNA SVs] - step.1: $file is not a vcf or vcf.gz format! Continue!"; next;
					}
					# assign to different caller category
					if ( $file =~/Manta\.vcf/ ) {
						$input_hash_dna{$sample}{'Manta'}[0] = $text;
						my %sss; $input_hash_dna{$sample}{'Manta'}[1] = \%sss; # define a reference to a hash structure
						DNAcaller::Menta_support($sample, $text, $input_hash_dna{$sample}{'Manta'}[1]);
					} elsif ( $file =~/Svaba\.vcf/ ) {
						$input_hash_dna{$sample}{'Svaba'}[0] = $text;
						my %sss; $input_hash_dna{$sample}{'Svaba'}[1] = \%sss; # define a reference to a hash structure
						DNAcaller::Svaba_support($sample, $text, $input_hash_dna{$sample}{'Svaba'}[1]);
					} elsif ( $file =~/Delly\.vcf/ ) {
						$input_hash_dna{$sample}{'Delly'}[0] = $text;
						my %sss; $input_hash_dna{$sample}{'Delly'}[1] = \%sss; # define a reference to a hash structure
						DNAcaller::Delly_support($sample, $text, $input_hash_dna{$sample}{'Delly'}[1]);
					} elsif ( $file =~/Gridss\.vcf/ ) {
						$input_hash_dna{$sample}{'Gridss'}[0] = $text;
						my %sss; $input_hash_dna{$sample}{'Gridss'}[1] = \%sss; # define a reference to a hash structure
						DNAcaller::Gridss_support($sample, $text, $input_hash_dna{$sample}{'Gridss'}[1]);
					} elsif ( $file =~/Lumpy\.vcf/ ) {
						$input_hash_dna{$sample}{'Lumpy'}[0] = $text;
						my %sss; $input_hash_dna{$sample}{'Lumpy'}[1] = \%sss; # define a reference to a hash structure
						DNAcaller::Lumpy_support($sample, $text, $input_hash_dna{$sample}{'Lumpy'}[1]);
					} else {
						print strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge DNA SVs] - step.1: Unknown $file, the caller is not in acceptable list!\n";
					}
				}
				close FL;
			} else {
				print strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge DNA SVs] - step.1: $sample is not a available folder\n";
			}
		}
		close DIR;

		# // [merge DNA SVs] - step.2-4: run SURVIVOR for merging variant callings from different callers...
		my %collect_support_dna;
		my %remove_duplicate; # remove duplicate in bedpe format
		foreach my $sample ( sort {$a cmp $b} keys %input_hash_dna ) {
			`mkdir $output/$sample`;
			my @caller_per; # caller order 'Delly', 'Gridss', 'Lumpy', 'Manta' and 'Svaba'
			foreach my $caller ( sort {$a cmp $b} keys %{$input_hash_dna{$sample}} ) {
				push @caller_per, $caller;
				open(TMP, ">$output/$sample/${caller}.vcf") || die strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge DNA SVs] - step.2: preparation of $sample for SURVIVOR ERROR:$!\n";
				print TMP $input_hash_dna{$sample}{$caller}[0];
				close TMP;
			
#				foreach my $id ( keys %{$input_hash_dna{$sample}{$caller}[1]} ) { 
#					print "[merge DNA SVs] - step.2: $sample\t$caller\t$id\t$input_hash_dna{$sample}{$caller}[1]{$id}[0]\t$input_hash_dna{$sample}{$caller}[1]{$id}[1]\t$input_hash_dna{$sample}{$caller}[1]{$id}[2]\n"; 
#				}
			}
			`ls $output/$sample/*.vcf > $output/${sample}_list`;
			# open an option for parameter settings - TO DO
			`SURVIVOR merge $output/${sample}_list 300 1 0 0 0 50 $output/${sample}_merge.vcf`;
			`SURVIVOR vcftobed $output/${sample}_merge.vcf 0 -1 $output/${sample}_merge.bedpe`;

			# // [merge DNA SVs] - step.2: process merged vcf file
			if ( ! -z "$output/${sample}_merge.vcf" ) { # if merged vcf file is OK
				open(IN, "grep -v '#' $output/${sample}_merge.vcf | cut -s -f3,9- |") || die strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge DNA SVs] - step.2: merged vcf file for $sample unavailable:$!\n";
				while ( <IN> ) {
					chomp $_;
					my @array = (split /\t/, $_);
					my $id = shift @array; # obtain id column of merged vcf
					my $format = shift @array; # obtain FORMAT column of merged vcf
					if ( scalar(@array) == scalar(@caller_per) ) { # How many callers shows positive in SV entry id
						my @merge_tag = split /\:/, $format;
						for (my $i=0; $i < scalar(@merge_tag); $i++) { # loop tag in FORMAT field of merged vcf file
							if ( $merge_tag[$i] eq 'ID' ) {
								for (my $j=0; $j < scalar(@array); $j++) { # loop caller in sample field 
									my @merge_info = split /\:/, $array[$j];
									if ( $merge_info[$i] eq 'NaN' || $merge_info[$i] eq 'NAN' ) {
										next;
									} else { # matched
										if ( exists($input_hash_dna{$sample}{$caller_per[$j]}) ) { # match to the conresponding caller
											my $name = $merge_info[$i];
											if ( $caller_per[$j] eq 'Manta' || $caller_per[$j] eq 'Svaba' ) { $name =~s/_/\:/g; } # replace '_' wiht ':' 
											if ( exists($input_hash_dna{$sample}{$caller_per[$j]}[1]{$name}) ) {
												push @{$collect_support_dna{$sample}{$id}{'split'}}, $input_hash_dna{$sample}{$caller_per[$j]}[1]{$name}[0];
												push @{$collect_support_dna{$sample}{$id}{'span'}}, $input_hash_dna{$sample}{$caller_per[$j]}[1]{$name}[1];
												push @{$collect_support_dna{$sample}{$id}{'type'}}, $input_hash_dna{$sample}{$caller_per[$j]}[1]{$name}[2];
												push @{$collect_support_dna{$sample}{$id}{'caller'}}, $caller_per[$j];
												$collect_support_dna{$sample}{$id}{'bedpe'} = undef;
#												print " [merge DNA SVs] - step.2: $sample\t$caller_per[$j]\t$input_hash_dna{$sample}{$caller_per[$j]}[1]{$name}[0]\t$input_hash_dna{$sample}{$caller_per[$j]}[1]{$name}[1]\t$input_hash_dna{$sample}{$caller_per[$j]}[1]{$name}[2]\n";
											} else {
												print "[merge DNA SVs] - step.2: $name of $caller_per[$j] is fitered out in initial control\n";
											}
										} else {
											die strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge DNA SVs] - step.2: the caller $caller_per[$j] does not match to merged vcf column, and it format incorrect!\n";
										}
									}
								}
								last;
							}
						}
					} else {
						die strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge DNA SVs] - step.2: the number of sample columns (", scalar(@array), "!=", scalar(@caller_per), ") incorrect in merged VCF!\n";
					}
				}
				close IN;
			} else {
				die strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge DNA SVs] - step.2: ERROR in sample $sample merging process - merged vcf unavailable!\n";
			}

			# // [merge DNA SVs] - step.3: process merged bedpe file
			my %single_break_dna; # - collect unique breakpoint entry 'chr:start-end'
			if ( ! -z "$output/${sample}_merge.bedpe" ) { # if merged bedpe file is OK
				open(IN, "grep -v '#' $output/${sample}_merge.bedpe | cut -s -f1-7,11 |") || die strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge DNA SVs] - step.3: merged bedpe file for $sample unavailable:$!\n";
				while ( <IN> ) {
					chomp $_;
					my ($chr1, $start1, $end1, $chr2, $start2, $end2, $id, $type) = (split /\t/, $_)[0, 1, 2, 3, 4, 5, 6, 7];
					if ( ! exists($remove_duplicate{$sample}{$chr1}{$start1}{$end1}{$chr2}{$start2}{$end2}) ) { # remove replicate entries in one sample
						$remove_duplicate{$sample}{$chr1}{$start1}{$end1}{$chr2}{$start2}{$end2} = 1;
					} else {
						next;
					}

					if ( exists($collect_support_dna{$sample}{$id}) ) { # if entry id in bedpe file exists also in vcf file
						if ( $start1 <= $end1 ) {
							if ( $start2 <= $end2 ) {
								$single_break_dna{$chr1}{$start1}{$end1} = undef;
								$single_break_dna{$chr2}{$start2}{$end2} = undef;
								if ( $chr1 eq $chr2 ) { # intra-chrom
									if ( ($start1 + $end1)/2 <= ($start2 + $end2)/2 ) {
										$collect_support_dna{$sample}{$id}{'bedpe'} = [$chr1, $start1, $end1, $chr2, $start2, $end2, $type];
									} else {
										$collect_support_dna{$sample}{$id}{'bedpe'} = [$chr2, $start2, $end2, $chr1, $start1, $end1, $type];
									}
								} else {
									$collect_support_dna{$sample}{$id}{'bedpe'} = [$chr1, $start1, $end1, $chr2, $start2, $end2, $type];
								}
							} else {
								$single_break_dna{$chr1}{$start1}{$end1} = undef;
								$single_break_dna{$chr2}{$end2}{$start2} = undef;
								if ( $chr1 eq $chr2 ) { # intra-chrom
									if ( ($start1 + $end1)/2 <= ($start2 + $end2)/2 ) {
										$collect_support_dna{$sample}{$id}{'bedpe'} = [$chr1, $start1, $end1, $chr2, $end2, $start2, $type];
									} else {
										$collect_support_dna{$sample}{$id}{'bedpe'} = [$chr2, $end2, $start2, $chr1, $start1, $end1, $type];
									}
								} else {
									$collect_support_dna{$sample}{$id}{'bedpe'} = [$chr1, $start1, $end1, $chr2, $end2, $start2, $type];
								}
							}
						} else {
							if ( $start2 <= $end2 ) {
								$single_break_dna{$chr1}{$end1}{$start1} = undef;
								$single_break_dna{$chr2}{$start2}{$end2} = undef;
								if ( $chr1 eq $chr2 ) { # intra-chrom
									if ( ($start1 + $end1)/2 <= ($start2 + $end2)/2 ) {
										$collect_support_dna{$sample}{$id}{'bedpe'} = [$chr1, $end1, $start1, $chr2, $start2, $end2, $type];
									} else {
										$collect_support_dna{$sample}{$id}{'bedpe'} = [$chr2, $start2, $end2, $chr1, $end1, $start1, $type];
									}
								} else {
									$collect_support_dna{$sample}{$id}{'bedpe'} = [$chr1, $end1, $start1, $chr2, $start2, $end2, $type];
								}
							} else {
								$single_break_dna{$chr1}{$end1}{$start1} = undef;
								$single_break_dna{$chr2}{$end2}{$start2} = undef;
								if ( $chr1 eq $chr2 ) { # intra-chrom
									if ( ($start1 + $end1)/2 <= ($start2 + $end2)/2 ) {
										$collect_support_dna{$sample}{$id}{'bedpe'} = [$chr1, $end1, $start1, $chr2, $end2, $start2, $type];
									} else {
										$collect_support_dna{$sample}{$id}{'bedpe'} = [$chr2, $end2, $start2, $chr1, $end1, $start1, $type];
									}
								} else {
									$collect_support_dna{$sample}{$id}{'bedpe'} = [$chr1, $end1, $start1, $chr2, $end2, $start2, $type];
								}
							}
						}
					} else {
						die strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge DNA SVs] - step.3: ERROR, entry id in bedpe does not exist in vcf file -- NOT make sense!\n";
					}
				}
				close IN;
			} else {
				die strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge DNA SVs] - step.3: ERROR in sample $sample merging process - merged bedpe unavailable!\n";
			}

			# // [merge DNA SVs] - step.4: assign gene annotation to breakpoints
			DNAcaller::dna_break_anno($sample, \%single_break_dna, \%gene_interval_dna, \%collect_support_dna, $output);

		}

		# // [merge DNA SVs] - step.5: final print merging results
		open (OUT, ">$output/Final_DNA_SVs.txt") || die "final print merging results of DNA SVs:$!\n";
		print OUT "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\ttype\tsplit\tspan\tgene1\tgene2\tcomment\n";
		foreach my $sample (keys %collect_support_dna ) {
			foreach my $id ( keys %{$collect_support_dna{$sample}} ) {
				if ( defined($collect_support_dna{$sample}{$id}{'bedpe'}) ) {
					my $chr1 = $collect_support_dna{$sample}{$id}{'bedpe'}[0];	my $chr2 = $collect_support_dna{$sample}{$id}{'bedpe'}[3];
					my $start1 = $collect_support_dna{$sample}{$id}{'bedpe'}[1];	my $start2 = $collect_support_dna{$sample}{$id}{'bedpe'}[4];
					my $end1 = $collect_support_dna{$sample}{$id}{'bedpe'}[2];	my $end2 = $collect_support_dna{$sample}{$id}{'bedpe'}[5];
					my @split = @{$collect_support_dna{$sample}{$id}{'split'}};	@split = sort @split;
					my @span = @{$collect_support_dna{$sample}{$id}{'span'}};	@span = sort @span;
					my $split_value;	my $span_value;
					if ( $support eq 'min' ) {
						$split_value = $split[0];	$span_value = $span[0];
					} elsif ( $support eq 'max' ) {
						$split_value = $split[-1];	$span_value = $span[-1];
					} elsif ( $support eq 'median' ) {
						my $index = int(scalar(@split)/2);
						if ( @split % 2 ) {
							$split_value = $split[$index];	$span_value = $span[$index];
						} else {
							$split_value = int(($split[$index-1] + $split[$index])/2);	$span_value = int(($span[$index-1] + $span[$index])/2);
						}
					} else {
						die "--support option should be one of [min, max, median]\n";
					}
					my $gene1;	my $gene2;
					if ( defined($collect_support_dna{$sample}{$id}{'bedpe'}[7]) ) { 
						foreach my $gene_tmp ( @{$collect_support_dna{$sample}{$id}{'bedpe'}[7]} ) {
							if ( $gene_tmp =~/^ENSG/ || $gene_tmp =~/^ENSMUS/ ) { $gene1 = $gene_tmp; } else { $gene1 = $gene_tmp; last; }
						}
					} else { 
						$gene1 = "*"; 
					}
					if ( defined($collect_support_dna{$sample}{$id}{'bedpe'}[8]) ) {
						foreach my $gene_tmp ( @{$collect_support_dna{$sample}{$id}{'bedpe'}[8]} ) {
							if ( $gene_tmp =~/^ENSG/ || $gene_tmp =~/^ENSMUS/ ) { $gene2 = $gene_tmp; } else { $gene2 = $gene_tmp; last; }
						}
					} else {
						$gene2 = "*";
					}

					if ( $chr1 ne $chr2 ) { # inter-chromosome events
						print OUT "$chr1\t$start1\t$end1\t$chr2\t$start2\t$end2\t$sample\tBND\t$split_value\t$span_value\t$gene1\t$gene2\t\n";
#						print "type:[@{$collect_support_dna{$sample}{$id}{'type'}}]\t";
					} else { # intra-chromosome events -- abount INS insertion????????
						my $type_assign;
						if ( $collect_support_dna{$sample}{$id}{'bedpe'}[6] eq 'BND' ) { # if incorrect type info
							foreach my $type_tmp ( @{$collect_support_dna{$sample}{$id}{'type'}} ) {
								if ( $type_tmp ne 'BND' ) { $type_assign = $type_tmp; last; }
							}
						} else {
							$type_assign = $collect_support_dna{$sample}{$id}{'bedpe'}[6];
						}
						
						if ( defined($type_assign) ) { 	
							print OUT "$chr1\t$start1\t$end1\t$chr2\t$start2\t$end2\t$sample\t$type_assign\t$split_value\t$span_value\t$gene1\t$gene2\t\n";
#							print "type:[@{$collect_support_dna{$sample}{$id}{'type'}}]\t";
						}
					}
				}
			}
		}
		close OUT;
	}

