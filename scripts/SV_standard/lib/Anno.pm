package Anno;
use strict;
use warnings;

	#----------------------------------------------------------------------------
	# load annotation file for RNA SV (i.e. gene and exon interval bed file)
	#----------------------------------------------------------------------------
	sub process_gene_rna {
		my ($path, $version, $file_ref, $exon_ref, $ens_ref, $offset) = @_;

		# load gene interval info
		open(RNA, "zcat $path/gene_interval_${version}.bed.gz |") || die strftime("%Y-%m-%d %H:%M:%S", localtime), " [Anno Class] - RNA: cannot open gene interval in this $path:$!\n";
		while ( <RNA> ) { 
			chomp $_; next if ( $_ =~/hromo/ );
			my ($chr, $start, $end, $ens_id, $strand) = (split /\t/, $_)[0, 1, 2, 3, 4];
			next if (! defined($ens_id) ); next if ( $ens_id eq ''); # if ensembl id is empty, just jump out
			# re-define gene interval [start-offset, end+offset] - default offset value: 1000 bp
			$start = $start - $offset;
			$end = $end + $offset;
			push @{$file_ref->{$chr}}, [$start, $end, $ens_id, $strand];
		}
                close RNA;

                # load exon interval info
                open(RNA, "zcat $path/exon_interval_${version}.bed.gz |") || die strftime("%Y-%m-%d %H:%M:%S", localtime), " [Anno Class] - RNA: cannot open exon interval in this $path:$!\n";
                while ( <RNA> ) { 
                        chomp $_; next if ( $_ =~/hromo/ );
                        my ($chr, $start, $end, $ens_id, $strand) = (split /\t/, $_)[0, 1, 2, 3, 4];
			$start = $start + 1; # convert bed start coordinate to normal coordinate
                        next if (! defined($ens_id) ); next if ( $ens_id eq ''); # if ensembl id is empty, just jump out
                        push @{$exon_ref->{$ens_id}}, [$chr, $start, $end, $strand];
                }
                close RNA;

                # load ensembl_id and gene symbol in a pair
		if ( $version eq 'hg19' || $version eq 'hg38' ) {
                	open(INFO, "$path/human_ensembl_symbol_1_1.txt") || die strftime("%Y-%m-%d %H:%M:%S", localtime), " [Anno Class] - RNA: cannot open human ensembl_id and gene_symbol in this $path:$!\n";
		} elsif ( $version eq 'GRCm39' ) {
			open(INFO, "$path/mouse_ensembl_symbol_1_1.txt") || die strftime("%Y-%m-%d %H:%M:%S", localtime), " [Anno Class] - RNA: cannot open mouse ensembl_id and gene_symbol in this $path:$!\n";
		}
		while ( <INFO> ) {
                        chomp $_; next if ( $_ =~/Gene_ID/ );
                        my ($ens_id, $gene_symbol, $type) = (split /\t/, $_)[0, 1, 2];
                        $ens_ref->{$ens_id} = [$gene_symbol, $type];
                }
                close INFO;
        }

	#------------------------------------------------------------------
	# load annotation file for SV DNA SV (i.e. gene interval bed file)
	#------------------------------------------------------------------
	sub process_gene_dna {
		my ($path, $version, $file_ref, $offset) = @_;

		# load ensembl_id and gene symbol in a pair
		my %ens_hash; #
		if ( $version eq 'hg19' || $version eq 'hg38' ) {
			open(INFO, "$path/human_ensembl_symbol_1_1.txt") || die strftime("%Y-%m-%d %H:%M:%S", localtime), " [Anno Class] - DNA: cannot open human ensembl_id and gene_symbol in this $path:$!\n";
		} elsif ( $version eq 'GRCm39' ) {
			open(INFO, "$path/mouse_ensembl_symbol_1_1.txt") || die strftime("%Y-%m-%d %H:%M:%S", localtime), " [Anno Class] - RNA: cannot open mouse ensembl_id and gene_symbol in this $path:$!\n";
		}
		while ( <INFO> ) {
			chomp $_; next if ( $_ =~/Gene_ID/ );
			my ($ens_id, $gene_symbol, $type) = (split /\t/, $_)[0, 1, 2];
			$ens_hash{$ens_id} = [$gene_symbol, $type];
		}
		close INFO;

		# load gene interval info
		open(DNA, "zcat $path/gene_interval_${version}.bed.gz |") || die strftime("%Y-%m-%d %H:%M:%S", localtime), " [Anno Class] - DNA: cannot open gene interval in this $path:$!\n";
		while ( <DNA> ) {
			chomp $_; next if ( $_ =~/hromo/ );
			my ($chr, $start, $end, $ens_id, $strand) = (split /\t/, $_)[0, 1, 2, 3, 4];
			next if (! defined($ens_id) ); next if ( $ens_id eq ''); # if ensembl id is empty, just jump out
			# re-define gene interval [if '+' strand -> start-offset; if '-' -> end+offset] - default offset value: 1000 bp
			if ( $strand eq '+' ) {
				$start = $start - $offset;
			} else {
				$end = $end + $offset;
			}
			if ( exists($ens_hash{$ens_id}) ) {
				if ( $ens_hash{$ens_id}[0] eq '' ) {
					push @{$file_ref->{$chr}}, [$start, $end, $ens_id, $ens_hash{$ens_id}[1]]; # print "ensembl_id: $start, $end, $ens_id, $ens_hash{$ens_id}[1]\n";
				} else {
					push @{$file_ref->{$chr}}, [$start, $end, $ens_hash{$ens_id}[0], $ens_hash{$ens_id}[1]]; # print "gene_symbol: $start, $end, $ens_hash{$ens_id}[0], $ens_hash{$ens_id}[1]\n";
				}
			}
		}
		close DNA;
	}
1;
