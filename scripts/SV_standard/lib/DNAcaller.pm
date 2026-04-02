package DNAcaller;
use strict;
use warnings;

	# // [load DNA SVs] from Lumpy caller
	sub Lumpy_support {
		my ($sample_name, $content, $hash_ref) = @_;
		my @all = split /\n/, $content;
		foreach my $line ( @all ) {
			next if ( $line =~/^#/ ); my $type = undef;
			my ($id, $info, $tag, $type1, $type2) = (split /\t/, $line)[2, 7, 8, 9, 10]; # print "[DNAcaller class] - Lumpy output for $sample_name: $id, $tag, $type1\n";
			if ( defined($info) ) {
				my $PR_num; # collect spanning read pairs number
				my $SR_num; # collect split read number
				if ( $info =~/SVTYPE\=([\w]+)/ ) { # match sv type
					$type = $1;
				} else {
					die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Lumpy output for $sample_name: SVTYPE info unavailable!\n";
				}
				if ( $info =~/\;SR\=([\d]+)/ ) { $SR_num = $1; }
				if ( $info =~/\;PE\=([\d]+)/ ) { $PR_num = $1; }
				if ( defined($PR_num) ) {
					if ( defined($SR_num) ) {
						next if ( $PR_num == 0 && $SR_num == 0 );
						$hash_ref->{$id} = [$SR_num, $PR_num, $type];
					} else {
						next if ( $PR_num == 0 );
						$hash_ref->{$id} = [0, $PR_num, $type];
					}
				} else {
					if ( defined($SR_num) ) {
						next if ( $SR_num == 0 );
						$hash_ref->{$id} = [$SR_num, 0, $type];
					} else {
						die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Lumpy output for $sample_name: num_split_read and num_spanning_read unavailable!\n";
					}
				}
			} else {
				die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Lumpy output for $sample_name: INFO format incorrect!\n";
			}
		}
	}

	# // [load DNA SVs] from Gridss caller
	sub Gridss_support {
		my ($sample_name, $content, $hash_ref) = @_;
		my @all = split /\n/, $content;
		foreach my $line ( @all ) {
			next if ( $line =~/^#/ ); my $type = undef;
			my ($id, $info, $tag, $type1, $type2) = (split /\t/, $line)[2, 7, 8, 9, 10]; # print "[DNAcaller class] - Gridss output for $sample_name: $id, $tag, $type1\n";
			if ( defined($info) ) {
				my $PR_num; # collect spanning read pairs number
				my $SR_num; # collect split read number
				my $BSC_num; # collect soft-clip read number
				my $BASSR_num; # collect read number for single breakpoint
				if ( $info =~/SVTYPE\=([\w]+)/ ) { # match sv type
					$type = $1;
				} else {
					die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Gridss output for $sample_name: SVTYPE info unavailable!\n";
				}
				if ( $info =~/\;SR\=([\d]+)/ ) { $SR_num = $1; }
				if ( $info =~/\;RP\=([\d]+)/ ) { $PR_num = $1; }
				if ( $info =~/\;BSC\=([\d]+)/ ) { $BSC_num = $1; }
				if ( $info =~/\;BASSR\=([\d]+)/ ) { $BASSR_num = $1; }
				if ( defined($PR_num) ) {
					if ( defined($SR_num) ) {
						if ( defined($BSC_num) ) {
							$SR_num = $SR_num + $BSC_num;
						}
						next if ( $PR_num == 0 && $SR_num == 0 );
						$hash_ref->{$id} = [$SR_num, $PR_num, $type];
					} else {
						if ( defined($BSC_num) ) {
							$SR_num = $BSC_num;
						} else {
							$SR_num = 0;
						}
						next if ( $PR_num == 0 && $SR_num == 0 );
						$hash_ref->{$id} = [$SR_num, $PR_num, $type];
					}
				} else {
					if ( defined($SR_num) ) {
						if ( defined($BSC_num) ) {
							$SR_num = $SR_num + $BSC_num;
						}
						next if ( $SR_num == 0 );
						$hash_ref->{$id} = [$SR_num, 0, $type];
					} else {
						if ( defined($BSC_num) ) {
							next if ( $BSC_num == 0 );
							$hash_ref->{$id} = [$BSC_num, 0, $type];
						} else {
							if ( defined($BASSR_num) ) {
								next if ( $BASSR_num == 0 );
								$hash_ref->{$id} = [$BASSR_num, 0, $type];
							} else {
								die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Gridss output for $sample_name: num_split_read and num_spanning_read unavailable!\n";
							}
						}
					}
				}
			} else {
				die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Gridss output for $sample_name: INFO format incorrect!\n";
			}
		}
	}

	# // [load DNA SVs] from Delly caller
	sub Delly_support {
		my ($sample_name, $content, $hash_ref) = @_;
		my @all = split /\n/, $content;
		foreach my $line ( @all ) {
			next if ( $line =~/^#/ ); my $type = undef;
			my ($id, $info, $tag, $type1, $type2) = (split /\t/, $line)[2, 7, 8, 9, 10]; # print "[DNAcaller class] - Delly output for $sample_name: $id, $tag, $type1\n";
			if ( $info =~/SVTYPE\=([\w]+)/ ) { # match sv type
				$type = $1;
			} else {
				die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Delly output for $sample_name: SVTYPE info unavailable!\n";
			}
			if ( defined($id) && defined($tag) && defined($type1) ) {
				my @element = split /\:/, $tag;
				my @first = split /\:/, $type1;
				my $PR_num; # collect spanning read pairs number
				my $SR_num; # collect split read number
				for (my $i = 0; $i < scalar(@element); $i++ ) {
					if ( $element[$i] eq 'RV' ) { # for split read
						if ( defined($type2) ) {
							my @second = split /\:/, $type2;
							if ( $first[$i] > $second[$i] ) { # first column -> tumor
								$SR_num = $first[$i];
							} else { # second column -> tumor
								$SR_num = $second[$i];
							}
						} else {
							 $SR_num = $first[$i];
						}
					} elsif ( $element[$i] eq 'DV' ) { # for spanning read
						if ( defined($type2) ) {
							my @second = split /\:/, $type2;
							if ( $first[$i] > $second[$i] ) { # first column -> tumor
								$PR_num = $first[$i];
							} else { # second column -> tumor
								$PR_num = $second[$i];
							}
						} else {
							$PR_num = $first[$i];
						}
					}
				}
				if ( defined($PR_num) && defined($SR_num) ) {
					next if ( $PR_num == 0 && $SR_num == 0 );
					$hash_ref->{$id} = [$SR_num, $PR_num, $type];
				} else {
					die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Delly output for $sample_name: num_split_read and num_spanning_read unavailable!\n";
				}
			} else {
				die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Delly output for $sample_name: INFO format incorrect!\n";
			}
		}
	}

	# // [load DNA SVs] from Svaba caller
	sub Svaba_support {
		my ($sample_name, $content, $hash_ref) = @_;
		my @all = split /\n/, $content;
		foreach my $line ( @all ) {
			next if ( $line =~/^#/ ); my $type = undef;
			my ($id, $info, $tag, $type1, $type2) = (split /\t/, $line)[2, 7, 8, 9, 10]; # print "[DNAcaller class] - Svaba output for $sample_name: $id, $tag, $type1\n";
			if ( $info =~/SVTYPE\=([\w]+)/ ) { # match sv type
				$type = $1;
			} else {
				die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Svaba output for $sample_name: SVTYPE info unavailable!\n";
			}
			if ( defined($id) && defined($tag) && defined($type1) ) {
				my @element = split /\:/, $tag;
				my @first = split /\:/, $type1;
				my $PR_num; # collect spanning read pairs number
				my $SR_num; # collect split read number
				for (my $i = 0; $i < scalar(@element); $i++ ) {
					if ( $element[$i] eq 'SR' ) { # for split read
						if ( defined($type2) ) {
							my @second = split /\:/, $type2;
							if ( $first[$i] > $second[$i] ) { # first column -> tumor
								$SR_num = $first[$i];
							} else { # second column -> tumor
								$SR_num = $second[$i];
							}
						} else {
							$SR_num = $first[$i];
						}
					} elsif ( $element[$i] eq 'DR' ) { # for spanning read
						if ( defined($type2) ) {
							my @second = split /\:/, $type2;
							if ( $first[$i] > $second[$i] ) { # first column -> tumor
								$PR_num = $first[$i];
							} else { # second column -> tumor
								$PR_num = $second[$i];
							}
						} else {
							$PR_num = $first[$i];
						}
					}
				}
				if ( defined($PR_num) && defined($SR_num) ) {
					next if ( $PR_num == 0 && $SR_num == 0 );
					$hash_ref->{$id} = [$SR_num, $PR_num, $type]; # print "[DNAcaller class] - Svaba output for $sample_name: $id $SR_num, $PR_num, $type\n";
				} else {
					die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Svaba output for $sample_name: num_split_read and num_spanning_read unavailable!\n";
				}
			} else {
				die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Svaba output for $sample_name: INFO format incorrect!\n";
			}
		}
	}
	
	# // [load DNA SVs] from Manta caller
	sub Menta_support {
		my ($sample_name, $content, $hash_ref) = @_;
		my @all = split /\n/, $content;
		foreach my $line ( @all ) {
			next if ( $line =~/^#/ ); my $type = undef;
			my ($id, $info, $tag, $type1, $type2) = (split /\t/, $line)[2, 7, 8, 9, 10]; # print "[DNAcaller class] - Manta output for $sample_name: $id, $tag, $type1\n";
			if ( $info =~/SVTYPE\=([\w]+)/ ) { # match sv type
				$type = $1;
			} else {
				die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Manta output for $sample_name: SVTYPE info unavailable!\n";
			}
			if ( defined($id) && defined($tag) && defined($type1) ) {
				if ( $tag =~/PR|SR/ ) {
					if ( $tag eq 'PR' ) { # only spanning read available
						if ( defined($type2) ) {
							my $num1 = (split /\,/, $type1)[1];
							my $num2 = (split /\,/, $type2)[1];
							if ( defined($num1) && defined($num2) ) {
								if ( $num1 > $num2 ) { $hash_ref->{$id} = [0, $num1, $type]; } else { $hash_ref->{$id} = [0, $num2, $type];  }
							} else {
								die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Manta output for $sample_name: num_spanning_read unavailable in normal and tumor!\n";
							}
						} else {
							my $num1 = (split /\,/, $type1)[1];
							if ( defined($num1) ) { $hash_ref->{$id} = [0, $num1, $type]; } else { die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Manta output for $sample_name: num_spanning_read unavailable in tumor!\n"; }
						}
					} elsif ( $tag eq 'PR:SR' ) {
						if ( defined($type2) ) { # both split and spanning read available
							my ($PR1, $SR1) = (split /\:/, $type1)[0, 1];
							my ($PR2, $SR2) = (split /\:/, $type2)[0, 1];
							my $PR_num1 = (split /\,/, $PR1)[1];
							my $SR_num1 = (split /\,/, $SR1)[1];
							my $PR_num2 = (split /\,/, $PR2)[1];
							my $SR_num2 = (split /\,/, $SR2)[1];

							if ( defined($PR_num1) && defined($PR_num2) && defined($SR_num1) && defined($SR_num2) ) {
								if ( $PR_num1 > $PR_num2 ) {
									if ( $SR_num1 > $SR_num2 ) {
										$hash_ref->{$id} = [$SR_num1, $PR_num1, $type];
									} else {
										$hash_ref->{$id} = [$SR_num2, $PR_num1, $type];
									}
								} else {
									if ( $SR_num1 > $SR_num2 ) {
										$hash_ref->{$id} = [$SR_num1, $PR_num2, $type];
									} else {
										$hash_ref->{$id} = [$SR_num2, $PR_num2, $type];
									}
								}
							} else {
								die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Manta output for $sample_name: num_split_read and num_spanning_read unavailable in normal and tumor!\n";
							}
						} else {
							my ($PR1, $SR1) = (split /\:/, $type1)[0, 1];
							my $PR_num1 = (split /\,/, $PR1)[1];
							my $SR_num1 = (split /\,/, $SR1)[1];
							if ( defined($PR_num1) && defined($SR_num1) ) { 
								$hash_ref->{$id} = [$SR_num1, $PR_num1, $type]; 
							} else { 
								die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Manta output for $sample_name: num_split_read and num_spanning_read unavailable in tumor!\n"; 
							}
						}
					} elsif ( $tag eq 'SR' ) { # only split read available
						if ( defined($type2) ) {
							my $num1 = (split /\,/, $type1)[1];
							my $num2 = (split /\,/, $type2)[1];
							if ( defined($num1) && defined($num2) ) {
								if ( $num1 > $num2 ) { $hash_ref->{$id} = [$num1, 0, $type]; } else { $hash_ref->{$id} = [$num2, 0, $type];  }
							} else {
								die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Manta output for $sample_name: num_split_read unavailable in normal and tumor!\n";
							}
						} else {
							my $num1 = (split /\,/, $type1)[1];
							if ( defined($num1) ) { $hash_ref->{$id} = [$num1, 0, $type]; } else { die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Manta output for $sample_name: num_split_read unavailable in tumor!\n"; }
						}
					} else {
						die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Manta output for $sample_name: split and spanning read info incorrect!\n";
					}
				} else {
					die strftime("%Y-%m-%d %H:%M:%S", localtime), " [DNAcaller class] - Manta output for $sample_name: INFO format incorrect!\n";
				}
			}
		}
	}

	#// assign gene annotation to the breakpoints
	sub dna_break_anno {
		my ($sample_name, $ref_single_break_dna, $ref_gene_interval_dna, $ref_collect_support_dna, $file) = @_;
		
		open (LOG, ">>$file/log") || die "finally print log history error:$!\n";
		foreach my $chr ( sort {$a cmp $b} keys %{$ref_single_break_dna} ) {
			my $i = 0;
			foreach my $start ( sort {$a <=> $b} keys %{$ref_single_break_dna->{$chr}} ) {
				foreach my $end ( sort {$a <=> $b} keys %{$ref_single_break_dna->{$chr}{$start}} ) {
					while ( $i < scalar(@{$ref_gene_interval_dna->{$chr}}) ) {
						# annotate breakpoint with gene region
						if ( ($start >= $ref_gene_interval_dna->{$chr}[$i][0] and $start <= $ref_gene_interval_dna->{$chr}[$i][1]) or ($end >= $ref_gene_interval_dna->{$chr}[$i][0] and $end <= $ref_gene_interval_dna->{$chr}[$i][1]) ) {
							my $j = $i;
							while ( $j < scalar(@{$ref_gene_interval_dna->{$chr}}) ) {
								if ( ($start >= $ref_gene_interval_dna->{$chr}[$j][0] and $start <= $ref_gene_interval_dna->{$chr}[$j][1]) or ($end >= $ref_gene_interval_dna->{$chr}[$j][0] and $end <= $ref_gene_interval_dna->{$chr}[$j][1]) ) {
									$ref_single_break_dna->{$chr}{$start}{$end}{$ref_gene_interval_dna->{$chr}[$j][2]} = $ref_gene_interval_dna->{$chr}[$j][3];
								} else {
#									print "[merge DNA SVs] - step.4: i: $i; j: $j; start: $start; [$gene_interval_dna{$chr}[$j][0], $gene_interval_dna{$chr}[$j][1] $gene_interval_dna{$chr}[$j][2]]\n";
									last;
								}
								$j++;
							}
							last;
						} else {
							if ( $end < $ref_gene_interval_dna->{$chr}[$i][0] ) {
								last;
							}
							if ( $start > $ref_gene_interval_dna->{$chr}[$i][1] ) {
								$i++;
							}
						}
					}

					# process breakpoint with >2 gene annotations
					if ( defined($ref_single_break_dna->{$chr}{$start}{$end}) ) {
						my @array = keys %{$ref_single_break_dna->{$chr}{$start}{$end}}; # get gene_name as array
						my $ref_hash = $ref_single_break_dna->{$chr}{$start}{$end}; # define a reference
						
						if ( scalar(@array) > 1 ) { # if a breakpoint with >2 gene annotations
#							foreach my $senz ( keys %{$ref_hash} ) { print "[merge DNA SVs] - step.4: unfiltered $sample_name, $chr, $start, $end, $senz, $ref_hash->{$senz}\n"; }
							
							# remove gene name with 'N-N': antise-RNA or presumed fusions
							delete @{$ref_hash}{ grep { $_ =~/[A-Za-z0-9]+\-[A-Za-z]/ } keys %{$ref_hash} };
							# remove common RNA genes: scRNA, vault_RNA, Mt_rDNA, sRNA, ribozyme, Mt_tRNA, rRNA, scaRNA, snRNA, misc_RNA
							delete @{$ref_hash}{ grep { $ref_hash->{$_} =~/scRNA|vault_RNA|Mt_rDNA|sRNA|ribozyme|Mt_tRNA|rRNA|scaRNA|snRNA|misc_RNA/ } keys %{$ref_hash} };
							# if $ref_hash is empty, set as 'undef'
							if (! %{$ref_hash} ) { $ref_single_break_dna->{$chr}{$start}{$end} = undef; next; }
							
							# remove snoRNA, miRNA, lncRNA and pseudogene if coding_protein annotations are present
							my %ref_hash1 = %{$ref_hash};
							delete @ref_hash1{ grep { $ref_hash1{$_} =~/snoRNA|miRNA|lncRNA|pseudogene/ } keys %ref_hash1 };
							my @array1 = keys %ref_hash1;
							if ( scalar(@array1) == 0 ) {
								print LOG "[merge DNA SVs] - step.4: kept snoRNA-miRNA-lncRNA-pseudogene $sample_name, $chr, $start, $end, ", keys %{$ref_hash}, "\n";
							} else {
								# remove 'orf' gene if other coding protein annotations are present
								my %ref_hash2 = %ref_hash1;
								delete @ref_hash2{ grep { $_ =~/orf/ } keys %ref_hash2 };
								my @array2 = keys %ref_hash2;
								if ( scalar(@array2) == 0 ) {
									$ref_single_break_dna->{$chr}{$start}{$end} = \%ref_hash1; # only 'orf' gene present, keep it
								} else {
									$ref_single_break_dna->{$chr}{$start}{$end} = \%ref_hash2; # both 'orf' and non-'orf' genes present, keep non-'orf' gene
								}
							}
						}
					}
				}
			}
		}

		# [merge DNA SVs] - step.4: link gene annotation to %collect_support_dna
		foreach my $id ( keys %{$ref_collect_support_dna->{$sample_name}} ) {
			if ( defined($ref_collect_support_dna->{$sample_name}{$id}{'bedpe'}) ) {
				my $chr1 = $ref_collect_support_dna->{$sample_name}{$id}{'bedpe'}[0];	my $chr2 = $ref_collect_support_dna->{$sample_name}{$id}{'bedpe'}[3];
				my $start1 = $ref_collect_support_dna->{$sample_name}{$id}{'bedpe'}[1];	my $start2 = $ref_collect_support_dna->{$sample_name}{$id}{'bedpe'}[4];
				my $end1 = $ref_collect_support_dna->{$sample_name}{$id}{'bedpe'}[2];	my $end2 = $ref_collect_support_dna->{$sample_name}{$id}{'bedpe'}[5];
				if ( exists($ref_single_break_dna->{$chr1}{$start1}{$end1}) and exists($ref_single_break_dna->{$chr2}{$start2}{$end2}) ) {
					if ( defined($ref_single_break_dna->{$chr1}{$start1}{$end1}) ) {
						my @partner1 = keys %{$ref_single_break_dna->{$chr1}{$start1}{$end1}};
						if ( defined($ref_single_break_dna->{$chr2}{$start2}{$end2}) ) {
							my @partner2 = keys %{$ref_single_break_dna->{$chr2}{$start2}{$end2}};
							$ref_collect_support_dna->{$sample_name}{$id}{'bedpe'}[7] = \@partner1;
							$ref_collect_support_dna->{$sample_name}{$id}{'bedpe'}[8] = \@partner2;
						} else {
							$ref_collect_support_dna->{$sample_name}{$id}{'bedpe'}[7] = \@partner1;
							$ref_collect_support_dna->{$sample_name}{$id}{'bedpe'}[8] = undef;
						}
					} else {
						if ( defined($ref_single_break_dna->{$chr2}{$start2}{$end2}) ) {
							my @partner2 = keys %{$ref_single_break_dna->{$chr2}{$start2}{$end2}};
							$ref_collect_support_dna->{$sample_name}{$id}{'bedpe'}[7] = undef;
							$ref_collect_support_dna->{$sample_name}{$id}{'bedpe'}[8] = \@partner2;
						} else {
							$ref_collect_support_dna->{$sample_name}{$id}{'bedpe'}[7] = undef;
							$ref_collect_support_dna->{$sample_name}{$id}{'bedpe'}[8] = undef;
						}
					}
				} else {
					die strftime("%Y-%m-%d %H:%M:%S", localtime), " [merge DNA SVs] - step.4: Unexpected $id in bedpe file fails to get gene annotation for sample $sample_name!\n";
				}
			}
		}
		close LOG;
	}
1;
