package RNAcaller;
use strict;
use warnings;

	# // parse RNA SVs from Dragen caller
	sub Dragen_support {
		my ($sample_name, $hash_ref, $content, $offset, $ref_gene_interval_rna, $ref_exon_interval_rna, $file) = @_;
		my @all = split /\n/, $content;
		for (my $i=0; $i < scalar(@all); $i++) {
			my ($tmp1, $tmp2, $gene1, $gene2, $split, $discordant) = (split /\t/, $all[$i])[0, 1, 2, 3, 4, 5];
			if ( $i == 0 ) { # judge header format
				if ( $tmp1 eq 'LeftBreakpoint' && $tmp2 eq 'RightBreakpoint' && $gene1 eq 'Gene1Id' && $gene2 eq 'Gene2Id' && $split eq 'NumSplitReads' && $discordant eq 'NumPairedReads' ) {
					next;
				} else {
					print "Dragen output header: {$tmp1, $tmp2, $gene1, $gene2, $split, $discordant}\n";
					die "[RNAcaller class] Header format of Dragen caller ouput not correct for $sample_name, and make sure the right version in usage!\n";
				}
			}
			next if ( $split < 2 ); # remove split number < 2
			my ($chr1, $breakpoint1, $strand1) = (split /\:/, $tmp1)[0, 1, 2];
			my ($chr2, $breakpoint2, $strand2) = (split /\:/, $tmp2)[0, 1, 2];
			# print "Dragen: $strand1, $strand2, $breakpoint1, $breakpoint2, $split, $discordant, $gene1, $gene2\n";
			my ($type1, $gene1_new, $break1_new) = &align_breakpoint($strand1, $chr1, $breakpoint1, $gene1, $ref_gene_interval_rna, $ref_exon_interval_rna, $offset, "upstream", $file);
			my ($type2, $gene2_new, $break2_new) = &align_breakpoint($strand2, $chr2, $breakpoint2, $gene2, $ref_gene_interval_rna, $ref_exon_interval_rna, $offset, "downstream", $file);
			
			if ( $type1 =~/boundary|exon|intron/ and $type2 =~/boundary|exon|intron/ ) {
				$gene1 = $gene1_new;	$gene2 = $gene2_new;
				next if ( $gene1 eq $gene2);
				if ( $type1 eq 'boundary' ) { $breakpoint1 = $break1_new; }
				if ( $type2 eq 'boundary' ) { $breakpoint2 = $break2_new; }
				$hash_ref->{$breakpoint1}{$breakpoint2} = [$chr1, $chr2, $strand1, $strand2, $gene1, $gene2, $split, $discordant];
			}
		}
	}

	# // parse RNA SVs from Arriba caller
	sub Arriba_support {
		my ($sample_name, $hash_ref, $content, $offset, $ref_gene_interval_rna, $ref_exon_interval_rna, $file) = @_;
		my @all = split /\n/, $content;
		for (my $i=0; $i < scalar(@all); $i++) {
			my ($strand1, $strand2, $coordinate1, $coordinate2, $split1, $split2, $discordant, $gene1, $gene2) = (split /\t/, $all[$i])[0, 1, 2, 3, 4, 5, 6, 7, 8];
			if ( $i == 0 ) { # judge header format
				if ( $strand1 eq 'strand1(gene/fusion)' && $strand2 eq 'strand2(gene/fusion)' && $coordinate1 eq 'breakpoint1' && $coordinate2 eq 'breakpoint2' && $split1 eq 'split_reads1' && $split2 eq 'split_reads2' && $discordant eq 'discordant_mates' && $gene1 eq 'gene_id1' && $gene2 eq 'gene_id2' ) {
					next;
				} else {
					print "Arriba output header: {$strand1, $strand2, $coordinate1, $coordinate2, $split1, $split2, $discordant, $gene1, $gene2}\n";
					die "[RNAcaller class] Header format of Arriba caller ouput not correct for $sample_name, and make sure the right version in usage!\n";
				}
			}
			my $split = $split1 + $split2;
			next if ( $split < 2 ); # remove split number < 2
			$gene1 =~s/\.[\d]+//g;  $gene2 =~s/\.[\d]+//g; # remove suffix of ensembl id
			$strand1 = (split /\//, $strand1)[1];   $strand2 = (split /\//, $strand2)[1];
			my ($chr1, $breakpoint1) = (split /\:/, $coordinate1)[0, 1];	$chr1 =~s/chr//g;	$chr1 = 'chr'.$chr1;
			my ($chr2, $breakpoint2) = (split /\:/, $coordinate2)[0, 1];	$chr2 =~s/chr//g;	$chr2 = 'chr'.$chr2;
			# print "Arriba: $strand1, $strand2, $coordinate1, $coordinate2, $split1, $split2, $discordant, $gene1, $gene2\n";
			my ($type1, $gene1_new, $break1_new) = &align_breakpoint($strand1, $chr1, $breakpoint1, $gene1, $ref_gene_interval_rna, $ref_exon_interval_rna, $offset, "upstream", $file);
			my ($type2, $gene2_new, $break2_new) = &align_breakpoint($strand2, $chr2, $breakpoint2, $gene2, $ref_gene_interval_rna, $ref_exon_interval_rna, $offset, "downstream", $file);

			if ( $type1 =~/boundary|exon|intron/ and $type2 =~/boundary|exon|intron/ ) {
				$gene1 = $gene1_new;    $gene2 = $gene2_new;
				next if ( $gene1 eq $gene2);
				if ( $type1 eq 'boundary' ) { $breakpoint1 = $break1_new; }
				if ( $type2 eq 'boundary' ) { $breakpoint2 = $break2_new; }
				$hash_ref->{$breakpoint1}{$breakpoint2} = [$chr1, $chr2, $strand1, $strand2, $gene1, $gene2, $split, $discordant];
			}
		}
	}

	# // parse RNA SVs from STAR-fusion caller
	sub STAR_fusion_support {
		my ($sample_name, $hash_ref, $content, $offset, $ref_gene_interval_rna, $ref_exon_interval_rna, $file) = @_;
		my @all = split /\n/, $content;
		for (my $i=0; $i < scalar(@all); $i++) {
			my ($split, $discordant, $gene1, $tmp1, $gene2, $tmp2) = (split /\t/, $all[$i])[0, 1, 2, 3, 4, 5];
			if ( $i == 0 ) { # judge header format
				if ( $split eq 'JunctionReadCount' && $discordant eq 'SpanningFragCount' && $gene1 eq 'LeftGene' && $gene2 eq 'RightGene' && $tmp1 eq 'LeftBreakpoint' && $tmp2 eq 'RightBreakpoint' ) {
					next;
				} else {
					print "STAR-fusion output header: {$split, $discordant, $gene1, $tmp1, $gene2, $tmp2}\n";
					die "[RNAcaller class] Header format of STAR-fusion caller ouput not correct for $sample_name, and make sure the right version in usage!\n";
				}
			}
			next if ( $split < 2 ); # remove split number < 2
			$gene1 = (split /\^/, $gene1)[1];	$gene1 =~s/\.[\d]+//g;
			$gene2 = (split /\^/, $gene2)[1];	$gene2 =~s/\.[\d]+//g;
			my ($chr1, $breakpoint1, $strand1) = (split /\:/, $tmp1)[0, 1, 2];	$chr1 =~s/chr//g;	$chr1 = 'chr'.$chr1;
			my ($chr2, $breakpoint2, $strand2) = (split /\:/, $tmp2)[0, 1, 2];	$chr2 =~s/chr//g;	$chr2 = 'chr'.$chr2;
			# print "STAR_fusion: $split, $discordant, $gene1, $tmp1, $gene2, $tmp2\n";
			my ($type1, $gene1_new, $break1_new) = &align_breakpoint($strand1, $chr1, $breakpoint1, $gene1, $ref_gene_interval_rna, $ref_exon_interval_rna, $offset, "upstream", $file);
			my ($type2, $gene2_new, $break2_new) = &align_breakpoint($strand2, $chr2, $breakpoint2, $gene2, $ref_gene_interval_rna, $ref_exon_interval_rna, $offset, "downstream", $file);

			if ( $type1 =~/boundary|exon|intron/ and $type2 =~/boundary|exon|intron/ ) {
				$gene1 = $gene1_new;    $gene2 = $gene2_new;
				next if ( $gene1 eq $gene2);
				if ( $type1 eq 'boundary' ) { $breakpoint1 = $break1_new; }
				if ( $type2 eq 'boundary' ) { $breakpoint2 = $break2_new; }
				$hash_ref->{$breakpoint1}{$breakpoint2} = [$chr1, $chr2, $strand1, $strand2, $gene1, $gene2, $split, $discordant];
			}
		}
	}

	# // parse RNA SVs from Fusioncatcher caller
	sub Fusioncatcher_support {
		my ($sample_name, $hash_ref, $content, $offset, $ref_gene_interval_rna, $ref_exon_interval_rna, $file) = @_;
		my @all = split /\n/, $content;
		for (my $i=0; $i < scalar(@all); $i++) {
			my ($discordant, $split, $tmp1, $tmp2, $gene1, $gene2) = (split /\t/, $all[$i])[0, 1, 2, 3, 4, 5];
			if ( $i == 0 ) { # judge header format
				if ( $discordant eq 'Spanning_pairs' && $split eq 'Spanning_unique_reads' && $tmp1 eq 'Fusion_point_for_gene_1(5end_fusion_partner)' && $tmp2 eq 'Fusion_point_for_gene_2(3end_fusion_partner)' && $gene1 eq 'Gene_1_id(5end_fusion_partner)' && $gene2 eq 'Gene_2_id(3end_fusion_partner)' ) {
					next;
				} else {
					print "Fusioncatcher output header: {$discordant, $split, $tmp1, $tmp2, $gene1, $gene2}\n";
					die "[RNAcaller class] Header format of Fusioncatcher caller ouput not correct for $sample_name, and make sure the right version in usage!\n";
				}
			}
			next if ( $split < 2 ); # remove split number < 2
			$gene1 =~s/\.[\d]+//g;	$gene2 =~s/\.[\d]+//g;
			my ($chr1, $breakpoint1, $strand1) = (split /\:/, $tmp1)[0, 1, 2];	$chr1 =~s/chr//g;	$chr1 = 'chr'.$chr1;
			my ($chr2, $breakpoint2, $strand2) = (split /\:/, $tmp2)[0, 1, 2];	$chr2 =~s/chr//g;	$chr2 = 'chr'.$chr2;
			# print "Fusioncatcher: $discordant, $split, $tmp1, $tmp2, $gene1, $gene2\n";
			my ($type1, $gene1_new, $break1_new) = &align_breakpoint($strand1, $chr1, $breakpoint1, $gene1, $ref_gene_interval_rna, $ref_exon_interval_rna, $offset, "upstream", $file);
			my ($type2, $gene2_new, $break2_new) = &align_breakpoint($strand2, $chr2, $breakpoint2, $gene2, $ref_gene_interval_rna, $ref_exon_interval_rna, $offset, "downstream", $file);

			if ( $type1 =~/boundary|exon|intron/ and $type2 =~/boundary|exon|intron/ ) {
				$gene1 = $gene1_new;	$gene2 = $gene2_new;
				next if ( $gene1 eq $gene2);
				if ( $type1 eq 'boundary' ) { $breakpoint1 = $break1_new; }
				if ( $type2 eq 'boundary' ) { $breakpoint2 = $break2_new; }
				$hash_ref->{$breakpoint1}{$breakpoint2} = [$chr1, $chr2, $strand1, $strand2, $gene1, $gene2, $split, $discordant];
			}
		}
	}

	# // arse RNA SVs from deFuse caller
	sub deFuse_support {
		my ($sample_name, $hash_ref, $content, $offset, $ref_gene_interval_rna, $ref_exon_interval_rna, $file) = @_;
                my @all = split /\n/, $content;
                for (my $i=0; $i < scalar(@all); $i++) {
                        my ($split, $split_span_p, $split_pos_p, $split_min_p, $homo, $gene1, $gene2, $chr1, $chr2, $breakpoint1, $breakpoint2, $strand1, $strand2, $span_mul, $span, $prob) = (split /\t/, $all[$i])[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
			if ( $i == 0 ) { # judge header format
				if ( $split eq 'splitr_count' && $split_span_p eq 'splitr_span_pvalue' && $split_pos_p eq 'splitr_pos_pvalue' && $split_min_p eq 'splitr_min_pvalue' && $homo eq 'breakpoint_homology' && $gene1 eq 'gene1' && $gene2 eq 'gene2' && $chr1 eq 'gene_chromosome1' && $chr2 eq 'gene_chromosome2' && $breakpoint1 eq 'genomic_break_pos1' && $breakpoint2 eq 'genomic_break_pos2' && $strand1 eq 'genomic_strand1' && $strand2 eq 'genomic_strand2' && $span_mul eq 'num_multi_map' && $span eq 'span_count' && $prob eq 'probability') {
					next;
				} else {
					print "deFuse output header: {$split, $split_span_p, $split_pos_p, $split_min_p, $homo, $gene1, $gene2, $chr1, $chr2, $breakpoint1, $breakpoint2, $strand1, $strand2, $span_mul, $span, $prob}\n";
					die "[RNAcaller class] Header format of deFuse caller ouput not correct for $sample_name, and make sure the right version in usage!\n";
				}
			}
			next if ( $split < 2 ); # remove split number < 2
			next if ( $split_span_p < 0.05 || $split_pos_p < 0.05 || $split_min_p < 0.05 || $prob < 0.05 );
			if ( $span ) { next if ( $span_mul/$span > 0.2 ); }
			next if ( $homo > 5 );
			$gene1 =~s/\.[\d]+//g;	$gene2 =~s/\.[\d]+//g;
			$chr1 =~s/chr//g;	$chr1 = 'chr'.$chr1;
			$chr2 =~s/chr//g;	$chr2 = 'chr'.$chr2;
			# print "deFuse: $split, $split_span_p, $split_pos_p, $split_min_p, $homo, $gene1, $gene2, $chr1, $chr2, $breakpoint1, $breakpoint2, $strand1, $strand2, $span_mul, $span\n";
			my ($type1, $gene1_new, $break1_new) = &align_breakpoint($strand1, $chr1, $breakpoint1, $gene1, $ref_gene_interval_rna, $ref_exon_interval_rna, $offset, "upstream", $file);
			my ($type2, $gene2_new, $break2_new) = &align_breakpoint($strand2, $chr2, $breakpoint2, $gene2, $ref_gene_interval_rna, $ref_exon_interval_rna, $offset, "downstream", $file);

			if ( $type1 =~/boundary|exon|intron/ and $type2 =~/boundary|exon|intron/ ) {
				$gene1 = $gene1_new;	$gene2 = $gene2_new;
				next if ( $gene1 eq $gene2);
				if ( $type1 eq 'boundary' ) { $breakpoint1 = $break1_new; }
				if ( $type2 eq 'boundary' ) { $breakpoint2 = $break2_new; }
				$hash_ref->{$breakpoint1}{$breakpoint2} = [$chr1, $chr2, $strand1, $strand2, $gene1, $gene2, $split, $span];
			}
		}
	}

	
	sub align_breakpoint {
		my ($s, $c, $bbb, $g, $gene_ref, $exon_ref, $offset, $comment, $file) = @_;
		my %b_pos_new; # - adjusted breakpoints and re-annotated partner genes
		open (LOG, ">>$file/log") || die "finally print log history error:$!\n";
		if ( $g =~/^ENSG/ || $g =~/^ENSMUS/ ) { # ens_id available
			if ( exists($exon_ref->{$g}) ) { # ens_id in caller output matches to that in annotation database
				my $start_new = $exon_ref->{$g}[0][1] - $offset;
				my $end_new = $exon_ref->{$g}[-1][2] + $offset;
				if ( $bbb >= $start_new and $bbb <= $end_new ) { # breakpoint wihin the given gene interval
					&collect_gene($bbb, $s, $gene_ref->{$c}, $exon_ref, \%b_pos_new, $comment, $g);
				} else { # if breakpoint not within the given gene interval -- need re-annotation
					if ( exists($gene_ref->{$c}) ) {
						&collect_gene($bbb, $s, $gene_ref->{$c}, $exon_ref, \%b_pos_new, $comment, "NOT");
					} else {
						print LOG "[RNAcaller class] given annotation has $g but breakpoint out of $g interval. Re-annotation does not work because chromosome $c not available in given annotation\n";
					}
				}
			} else { # if breakpoint has a gene_id which is not avaialable in our database -- need re-annotation
				if ( exists($gene_ref->{$c}) ) {
					&collect_gene($bbb, $s, $gene_ref->{$c}, $exon_ref, \%b_pos_new, $comment, "NOT");
				} else {
					print LOG "[RNAcaller class] given annotation does not have $g. Re-annotation does not work because chromosome $c not available in given annotation\n";
				}
			}
		} else { # ens_id not available -- need re-annotation
			if ( exists($gene_ref->{$c}) ) {
				&collect_gene($bbb, $s, $gene_ref->{$c}, $exon_ref, \%b_pos_new, $comment, "NOT");
			} else {
				print LOG "[RNAcaller class] breakpoint is annotated at intergenic region in fusion caller output. Re-annotation does not work because chromosome $c not available in given annotation\n";
			}
		}

		my $ens_return; my $bbb_return;
		# // prioritisation order: boundary -> exon -> intron -> intergenic
		if ( %b_pos_new ) {
			if ( exists($b_pos_new{'boundary'}) ) {
				print LOG "[RNAcaller class] boundary{$g-$bbb}: "; 
				foreach my $ens_id ( keys %{$b_pos_new{'boundary'}} ) { 
					print LOG "$ens_id=>";	$ens_return = $ens_id;
					foreach my $site ( reverse sort {$b_pos_new{'boundary'}{$ens_id}->{$a} <=> $b_pos_new{'boundary'}{$ens_id}->{$b}} keys %{$b_pos_new{'boundary'}{$ens_id}} ) {
						print LOG "($site - $b_pos_new{'boundary'}{$ens_id}{$site}) ";	$bbb_return = $site;
					}
				} 
				print LOG "\n";
				return("boundary", $ens_return, $bbb_return);
			} else {
				if ( exists($b_pos_new{'exon'}) ) {
					print LOG "[RNAcaller class] exon{$g-$bbb}: ";
					foreach my $ens_id ( keys %{$b_pos_new{'exon'}} ) { 
						print LOG "$ens_id=>";	$ens_return = $ens_id;
						foreach my $site ( sort {$b_pos_new{'exon'}{$ens_id}->{$a} <=> $b_pos_new{'exon'}{$ens_id}->{$b}} keys %{$b_pos_new{'exon'}{$ens_id}} ) {
							print LOG "($site - $b_pos_new{'exon'}{$ens_id}{$site}) ";	$bbb_return = $site;
						}
					}
					print LOG "\n";
					return("exon", $ens_return, $bbb_return);
				} else {
					if ( exists($b_pos_new{'intron'}) ) {
						print LOG "[RNAcaller class] intron{$g-$bbb}: ";
						foreach my $ens_id ( keys %{$b_pos_new{'intron'}} ) { 
							print LOG "($ens_id=>undef) ";	$ens_return = $ens_id;
						}
						print LOG "\n";
						return("intron", $ens_return, undef);
					} else {
						if ( exists($b_pos_new{'intergenic'}) ) {
							print LOG "[RNAcaller class] intergenic{$g-$bbb}: \n";
							return("intergenic", undef, undef);
						} else {
							die "Error: unexpeced annotation type for [$g-$bbb]\n";
						}
					}
				}
			}
		} else {
			return("unknown", undef, undef);
		}
		close LOG;
	}

	sub collect_gene { # get a set of genes harbor the given breakpoint
		my ($b_pos, $strand, $ref, $ref_exon, $ref_anno, $comment, $gene_id) = @_;
		my $i = 0;
		my %collect_match; # gene strand matches to that of breakpoint
		my %collect_no; # gene strand does not matches to that of breakpoint

		if ( $gene_id eq "NOT" ) { # start to re-annotate: assign new ens_id
			while ( $i < scalar(@{$ref}) ) {
				if ( $b_pos >= $ref->[$i][0] and $b_pos <= $ref->[$i][1] ) {
					my $j = $i;
					while ( $j < scalar(@{$ref}) ) {
						if ( $b_pos >= $ref->[$j][0] and $b_pos <= $ref->[$j][1] ) {
							if ( $strand eq $ref->[$j][3] ) { # strand of breakpoint == gene strand
								$collect_match{$ref->[$j][2]} = $ref->[$j][3]; # ens_id => strand
							} else {
								$collect_no{$ref->[$j][2]} = $ref->[$j][3]; # ens_id => strand
							}	
						} else {
							last;
						}
						$j++;
					}
					last;
				} else {
					if ( $b_pos < $ref->[$i][0] ) {
						last;
					}
					if ( $b_pos > $ref->[$i][1] ) {
						$i++;
					}
				}
			}
		} else {
			if ( $strand eq $ref_exon->{$gene_id}[0][3] ) { # strand of breakpoint == gene strand
				$collect_match{$gene_id} = $ref_exon->{$gene_id}[0][3];
			} else {
				$collect_no{$gene_id} = $ref_exon->{$gene_id}[0][3];
			}
		}
			
		if ( %collect_match ) {
			foreach my $id ( keys %collect_match ) {
				if ( exists($ref_exon->{$id}) ) {
					foreach my $sss ( @{$ref_exon->{$id}} ) { # each exon interval
						if ( ($sss->[2] - $sss->[1] + 1) > 5 ) { # if exon width > 5
							if ( $comment eq 'upstream' ) {
								if ( $sss->[3] eq '+' ) { # '+' strand
									if ( $b_pos >= ($sss->[2] - 5) and $b_pos <= ($sss->[2] + 5) ) { # breakpoint adjustment to exon boundary
										$ref_anno->{'boundary'}{$id}->{$sss->[2]} = abs($sss->[2] - $b_pos);
									} elsif ( $b_pos >= $sss->[1] and $b_pos < ($sss->[2] - 5) ) {
										$ref_anno->{'exon'}{$id}->{$b_pos}++;
									}
								} else {
									if ( $b_pos >= ($sss->[1] - 5) and $b_pos <= ($sss->[1] + 5) ) { # breakpoint adjustment to exon boundary
										$ref_anno->{'boundary'}{$id}->{$sss->[1]} = abs($sss->[1] - $b_pos);
									} elsif ( $b_pos > ($sss->[1] + 5) and $b_pos <= $sss->[2] ) { 
										$ref_anno->{'exon'}{$id}->{$b_pos}++;
									}
								}
							} else { # for downstream partner
								if ( $sss->[3] eq '+' ) { # '+' strand
									if ( $b_pos >= ($sss->[1] - 5) and $b_pos <= ($sss->[1] + 5) ) { # breakpoint adjustment to exon boundary
										$ref_anno->{'boundary'}{$id}->{$sss->[1]} = abs($sss->[1] - $b_pos);
									} elsif ( $b_pos > ($sss->[1] + 5) and $b_pos <= $sss->[2] ) {
										$ref_anno->{'exon'}{$id}->{$b_pos}++;
									}
								} else {
									if ( $b_pos >= ($sss->[2] - 5) and $b_pos <= ($sss->[2] + 5) ) { # breakpoint adjustment to exon boundary
										$ref_anno->{'boundary'}{$id}->{$sss->[2]} = abs($sss->[2] - $b_pos);
									} elsif ( $b_pos >= $sss->[1] and $b_pos < ($sss->[2] - 5) ) {
										$ref_anno->{'exon'}{$id}->{$b_pos}++;
									}
								}
							}
						} else {
							if ( $comment eq 'upstream' ) {
								if ( $sss->[3] eq '+' ) { # '+' strand
									if ( $b_pos >= $sss->[1] and $b_pos <= ($sss->[2] + 5) ) { # breakpoint adjustment to exon boundary
										$ref_anno->{'boundary'}{$id}->{$sss->[2]} = abs($sss->[2] - $b_pos);
									}
								} else {
									if ( $b_pos >= ($sss->[1] - 5) and $b_pos <= $sss->[2] ) { # breakpoint adjustment to exon boundary
										$ref_anno->{'boundary'}{$id}->{$sss->[1]} = abs($sss->[1] - $b_pos);
									}
								}
							} else { # for downstream partner
								if ( $sss->[3] eq '+' ) { # '+' strand
									if ( $b_pos >= ($sss->[1] - 5) and $b_pos <= $sss->[2] ) { # breakpoint adjustment to exon boundary
										$ref_anno->{'boundary'}{$id}->{$sss->[1]} = abs($sss->[1] - $b_pos);
									}
								} else {
									if ( $b_pos >= $sss->[1] and $b_pos <= ($sss->[2] + 5) ) { # breakpoint adjustment to exon boundary
										$ref_anno->{'boundary'}{$id}->{$sss->[2]} = abs($sss->[2] - $b_pos);
									}
								}
							}
						}
					}

					if ( exists($ref_anno->{'boundary'}) or exists($ref_anno->{'exon'}) ) {
					} else {
						$ref_anno->{'intron'}{$id} = undef;
					}
				}
			}
		} else {
			if ( %collect_no ) {
				foreach my $id ( keys %collect_no ) {
					if ( exists($ref_exon->{$id}) ) {
						foreach my $sss ( @{$ref_exon->{$id}} ) {
							if ( $b_pos >= $sss->[1] and $b_pos <= $sss->[2] ) {
								$ref_anno->{'exon'}{$id}->{$b_pos}++;
							}
						}
						if (! exists($ref_anno->{'exon'}) ) {
							$ref_anno->{'intron'}{$id} = undef;
						}
					}
				}
			} else {
				$ref_anno->{'intergenic'} = undef;
			}
		}
	}
1;
