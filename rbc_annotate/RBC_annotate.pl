#!/usr/bin/perl -w

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;

use Getopt::Long;
use FindBin;
use lib $FindBin::Bin;

#use SNPlist;
#use IndelList;
use GeneSNPlist;
use Varlist;
use List::Util qw[min max];

### Variables
my $refseq_file = "";
my $refseq_error_file = "";
my $gff = "";
my $snp_file = "";
my $trans_file = "";
my $chromosome;
my $start;
my $end;
my $target_region = 0;
my $outfile_prefix = "";
my $all_isoforms = 0;
my $cds_alignments = 0;
my $vcf = 1;
my $mut = 0;
my $majmin_allele = "all";
my @filter = ();

########### Configs #######################
my $config_file = "";
# Rice Genbank defaults:
my $gff_sorted = 1;

	my %gene_hash = ( "CDS"=>"CDS", );
my %gene_features = ( "gene"=>\%gene_hash, );
my %coding_feature = ( "gene"=>"CDS", );
my %other_features = ();

# chromosome name adjustments
my $crop_gff_start = "";
my $crop_gff_end = "";

#determination of gene name
my $gene_tag = "ID";
my $gene_splitstring = "";

my $coding_tag = "Parent";
my $coding_splitstring = "";

my $other_tag = "Parent";
my $splitstring = "";

# determination of isoform
my $isoform_feature = "mRNA";
my $isoform_tag = "ID";
my $isoform_splitstring = "t2";

my %coding_subfeatures = ();
###########################################


my @args = @ARGV;

### Get command line options
my %CMD;
GetCom();

my $filepos_var = 0;
my $filepos_gff = 0;
my $filepos_genome = 0;


my $file_ending = $mut ? "mut" : "rbca";
open SNPOUT, ">>$outfile_prefix.$file_ending" or die "Cannot open output file $outfile_prefix.$file_ending\n";
open GENEOUT, ">$outfile_prefix.genes" or die "Cannot open output file $outfile_prefix.genes\n" if ($cds_alignments);

## Store header if output is rbca
if (!$mut) {
	system("grep \"^##\" $snp_file > $outfile_prefix.rbca");
	print SNPOUT "##$0 ".join(" ",@args)."\n";
}


### Parse chromosome translation file
my %trans = ();
my %trans_rev = ();
if ($trans_file ne "") {
	open TRANS, $trans_file or die "Cannot open chromosome translation file: $trans_file\n";
	while (<TRANS>) {
		my ($chr1, $chr2) = split/\t/;
		$trans{$chr2} = $chr1;
		$trans_rev{$chr1} = $chr2;
	}

	### Translate chrom names of snp file
	print STDERR "Translating chromosome names of input files...\n";
	if ($snp_file ne "") {
		open F,"<$snp_file";
		open OUT, ">$snp_file.TMP.transl";
		while (<F>) {
			my @a = split/\t/;
			if (exists $trans{$a[0]}) {
				$a[0] = $trans{$a[0]};
			}
			else { die "Could not find chromosome '".$a[0]."' of snp file '$snp_file' in chromosome translation file '$trans_file'! Aborted\n"; }

			unshift(@a, "S");
			print OUT join("\t", @a);
		}
		close F; close OUT;
		$snp_file = "$snp_file.TMP.transl";
	}
}


### Sort snp file
print STDERR "Sorting input files...\n";
system("sort -g -k1 -k2 $snp_file | awk '{ if (substr(\$1,1,1) == \"#\") { next; } str=\"S\t\"; if (length(\$4)>1 && index(\$4,\",\")==0) { str=\"D\t\" } else if (length(\$5)>1 && index(\$5,\",\")==0) { str=\"I\t\" } print str \$0}' > VAR.TMP");

my $var_file = "VAR.TMP";

### clean up
system("rm $snp_file.TMP*") if ($trans_file ne "");


my $peak_pos = 1;

if ($gff ne "") {

	### Read in config file:
	if ($config_file ne "") {

		# resetting hashes:
		%gene_features = ();
		%coding_feature = ();
		%other_features = ();

		open CONFIG, "<$config_file";
		while (<CONFIG>) {

			#comment and empty lines
			next if ($_ =~ m/^#/ || $_ =~ m/^\s+$/);

			chomp;
			my @a = split/\t/;
			if ($a[0] eq "GENE_FEATURE") {
				if ($a[1] =~ m/#/) { print STDERR "Alias is not allowed for gene features!\n"; exit; }

				for (my $i=2; $i!=@a; ++$i) {
					if ($a[$i] =~ m/#/) {
						my ($feature, $name) = split"#",$a[$i];
						$gene_features{$a[1]}{$feature} = $name;
					}
					else {
						$gene_features{$a[1]}{$a[$i]} = $a[$i];
					}
				}
			}
			elsif ($a[0] eq "CODING_FEATURE") {
				if ($a[1] =~ m/#/) {
					my ($feature, $name) = split"#", $a[1];
					$coding_feature{$feature} = $name;
				}
				else {
					$coding_feature{$a[1]} = $a[2];
				}
			}
			elsif ($a[0] eq "OTHER_FEATURES") {
				for (my $i=1; $i!=@a; ++$i) {
					if ($a[$i] =~ m/#/) {
						my ($feature, $name) = split"#", $a[$i];
						$other_features{$feature} = $name;
					}
					else {
						$other_features{$a[$i]} = $a[$i];
					}
				}
			}
			elsif ($a[0] eq "CROP_CHR_GFF_START") {
				$crop_gff_start = defined $a[1] ? $a[1] : "";
			}
			elsif ($a[0] eq "CROP_CHR_GFF_END") {
				 $crop_gff_end = defined $a[1] ? $a[1] : "";
			}
			elsif ($a[0] eq "ISOFORM_FEATURE") {
                                $isoform_feature = defined $a[1] ? $a[1] : "";
                        }
			elsif ($a[0] eq "GENE_SPLITSTRING") {
				$gene_splitstring = defined $a[1] ? $a[1] : "";
			}
			elsif ($a[0] eq "ISOFORM_SPLITSTRING") {
                                $isoform_splitstring = defined $a[1] ? $a[1] : "";
                        }
			elsif ($a[0] eq "CODING_SPLITSTRING") {
				$coding_splitstring = defined $a[1] ? $a[1] : "";
			}
			elsif ($a[0] eq "SPLITSTRING") {
				$splitstring = defined $a[1] ? $a[1] : "";
			}
                        elsif ($a[0] eq "GENE_TAG") {
                                $gene_tag = defined $a[1] ? $a[1] : "";
                        }
                        elsif ($a[0] eq "ISOFORM_TAG") {
                                $isoform_tag = defined $a[1] ? $a[1] : "";
                        }
                        elsif ($a[0] eq "CODING_TAG") {
                                $coding_tag = defined $a[1] ? $a[1] : "";
                        }
                        elsif ($a[0] eq "TAG") {
                                $other_tag = defined $a[1] ? $a[1] : "";
                        }
                        elsif ($a[0] eq "GFF_SORTED") {
				$gff_sorted = $a[1];
			}
			else {
				print STDERR "Do not recognize config variable $a[0]\n";
				exit;
			}			
		}
		close CONFIG;

	}

}

foreach (keys %coding_feature) {
	$coding_subfeatures{$coding_feature{$_}} = 1;
}



### annotate for each chromosome unless target region was specified
my $next_chr = 1;
while ($next_chr != 0) {

	### Get Variation lists
	my $vars = new Varlist();
	if ($target_region) {
		$filepos_var = $vars->get_and_prioritize($chromosome, $start, $end, $var_file, $peak_pos, $filepos_var, $vcf, $majmin_allele, \@filter);

		if (scalar(@{$vars->{vars}}) == 0) {
			$next_chr = 0;
			last;
		}
	}
	else {
		$filepos_var = $vars->get_and_prioritize("-1", -1, -1, $var_file, $peak_pos, $filepos_var, $vcf, $majmin_allele, \@filter);
		if (scalar(@{$vars->{vars}}) == 0) {
			$next_chr = 0;
			last;
		}
		$chromosome = $vars->{vars}[0]->{chromosome};

		$start = $vars->{vars}[0]->{position};
		$end = $vars->{maxend};
	}


	my $var_array_pos = 0;
	my $var_maxnr = scalar(@{$vars->{vars}});

	if ($trans_file ne "") {
		print "chr ".$trans_rev{$chromosome}." #Variations: $var_maxnr\n";
	}
	else {
		print "chr $chromosome #Variations: $var_maxnr\n";
	}

	### Functional analysis of SNPs and indels
	my %coding_ann = ();
	my %gene_ann = ();
	my %seq_type = ();
	my %gene_type = ();
	my %feat_by_pos = ();
	my %genes = ();
	my $chr = "";

	if( ($gff ne "") && ($refseq_file ne "") ) {

		open GFF, "<$gff";

		### Get gene annotation from gff for $chromosome
		my $region_visited = 0;
		my $mainfeature = "";
		my $gene_name = "";
		my $isoform = "";

		my $gff_chr   = "";
		my $gff_start = 0;
		my $gff_end   = 0;
		my $gff_type  = "";
		my $gff_ori   = "";
		my $gff_tags  = "";

		while( <GFF> ) {

			chomp;
			next if ($_ =~ m/^#/);

			# GFF format: chr, source, seq_type, start, end, score, orientation, frame, description
			($gff_chr, $gff_type, $gff_type, $gff_start, $gff_end, $gff_ori, $gff_ori, $gff_tags, $gff_tags) = split("\t", $_);
	
			# Adapt equal chromosome names
			if ($crop_gff_start ne "" && $gff_chr =~ m/^$crop_gff_start/) { $gff_chr = substr($gff_chr, length($crop_gff_start)); }
			if ($crop_gff_end ne "" && $gff_chr =~ m/$crop_gff_end$/) { $gff_chr = substr($gff_chr, 0, length($gff_chr)-length($crop_gff_end)); }
			if ($trans_file ne "") {
				if (exists $trans{$gff_chr}) { $gff_chr = $trans{$gff_chr}; }
				else { die "Could not find chromosome '$gff_chr' of GFF file '$gff' in chromosome translation file '$trans_file'! Aborted\n"; }
			}

			last if ($region_visited && $gff_chr ne $chromosome);

			# Read annotation
			if( ($gff_chr eq $chromosome) && (!$gff_sorted || ($gff_end >= $start && $gff_start <= $end)) ) {
				$start = min($gff_start, $start);
				$end   = max($gff_end, $end);

				$region_visited = 1;

				# Determine gene name or isoform:
				my @tags = split/;/, $gff_tags;
				foreach my $tag_str (@tags) {
					if ($tag_str =~ m/=/) {
						my @t = split/=/, $tag_str;
						if ( $gff_type eq $isoform_feature && $t[0] eq $isoform_tag ) {
							for (my $i=0; $i<length($isoform_splitstring); $i+=2) {
								my $splitchar = substr($isoform_splitstring, $i, 1);
								my $splitcol  = substr($isoform_splitstring, $i+1, 1)-1;
								$t[1] = int((split/$splitchar/, $t[1])[$splitcol]);
							}
							if ($t[1] !~ m/^fgene/ && $t[1] =~ m/\./) {
								$isoform = (split/\./, $t[1])[-1];
							}
							else {
								$isoform = 1;
							}
						}

						$isoform = 1 if ($isoform eq "");

						if ( exists $coding_feature{$mainfeature} && $coding_feature{$mainfeature} eq $gff_type && $t[0] eq $coding_tag ) {
							for (my $i=0; $i!=length($coding_splitstring); $i+=2) {
								my $splitchar = substr($coding_splitstring, $i, 1);
								my $splitcol  = substr($coding_splitstring, $i+1, 1)-1;
								$t[1] = (split/$splitchar/, $t[1])[$splitcol];
							}

							$t[1] = substr( $t[1], 0, rindex($t[1], "\." )) if ($t[1] =~ m/\./);

							if (defined $gene_name && $t[1] ne $gene_name) { $mainfeature = ""; $isoform = ""; }
							$gene_name = $t[1];
						}
						elsif ( exists $gene_features{$gff_type} && $t[0] eq $gene_tag ) {
							for (my $i=0; $i!=length($gene_splitstring); $i+=2) {
								my $splitchar = substr($gene_splitstring, $i, 1);
								my $splitcol  = substr($gene_splitstring, $i+1, 1)-1;

								$t[1] = (split/$splitchar/, $t[1])[$splitcol];
							}

							$t[1] = substr( $t[1], 0, rindex($t[1], "\." )) if ($t[1] =~ m/\./);

							if (defined $gene_name && $t[1] ne $gene_name) { $mainfeature = ""; $isoform = ""; }
							$gene_name = $t[1];
						}
						elsif ( (exists $other_features{$gff_type} || (exists $gene_features{$mainfeature} && exists $gene_features{$mainfeature}{$gff_type})) && $t[0] eq $other_tag ) {
							for (my $i=0; $i!=length($splitstring); $i+=2) {
								my $splitchar = substr($splitstring, $i, 1);
								my $splitcol  = substr($splitstring, $i+1, 1)-1;
								$t[1] = (split/$splitchar/, $t[1])[$splitcol];
							}

							$t[1] = substr( $t[1], 0, rindex($t[1], "\." )) if ($t[1] =~ m/\./);

							if (defined $gene_name && $t[1] ne $gene_name) { $mainfeature = ""; $isoform = ""; }
							$gene_name = $t[1];
						}
					}
				}

				next if ($isoform ne "" && !$all_isoforms && $isoform ne "1");

				# Gene locus (= mainfeature)
				if (exists $gene_features{$gff_type} && !exists $coding_ann{$gene_name}) {
					my @gene_locus = ($gff_start, $gff_end, $gff_ori, "");
					$gene_ann{$gene_name} = \@gene_locus;

					$gene_type{$gene_name} = $gff_type;
					$feat_by_pos{$gff_start} = $gene_name."#".$gff_end;
					$mainfeature = $gff_type;
				}
				# Coding sequence
				elsif (exists $coding_feature{$mainfeature} && $gff_type eq $coding_feature{$mainfeature}) {
					if(!exists $coding_ann{$gene_name}) {
						my %gene = ();
						$coding_ann{$gene_name} = \%gene;
					}

					my @cds = ($gff_start, $gff_end, $gff_ori, "");
					$coding_ann{$gene_name}{$gff_start} = \@cds;
					$seq_type{$gene_name}{$gff_start}{"type"} = $coding_feature{$mainfeature};
					$seq_type{$gene_name}{$gff_start}{"end"} = $gff_end;
				}
				# Other sequence types
				elsif (exists $other_features{$gff_type} || (exists $gene_features{$mainfeature} && exists $gene_features{$mainfeature}{$gff_type})) {
					my $seq_type = "";
					$seq_type = $other_features{$gff_type} if (exists $other_features{$gff_type});
					$seq_type = $gene_features{$mainfeature}{$gff_type} if (exists $gene_features{$mainfeature} && exists $gene_features{$mainfeature}{$gff_type});

					$seq_type{$gene_name}{$gff_start}{"type"} = $seq_type;
					$seq_type{$gene_name}{$gff_start}{"end"} = $gff_end;
					$feat_by_pos{$gff_start} = $gene_name."#".$gff_end if ($mainfeature eq "");
				}
			} # end of 'gene is in range'
		} # end of while <GFF>

		close GFF;


		if (scalar(keys %feat_by_pos) == 0) {
			print STDERR "\n***WARNING***\nNo features for chromosome '".($trans_file ne ""?$trans_rev{$chromosome}:$chromosome)."' found!\nWill skip this chromosome!\nPossible reasons: Discordant chromosome name or \nno specified features detected in gff file\n***********\n";
                        next;
                }


		## Load chromosome
		open GENOME, $refseq_file or die "Cannot open reference sequence file\n";
		my $chr_seq = "";
		while( <GENOME> ) {
			chomp;
			if(substr($_, 0, 1)  eq ">") {
				my $current_chr = substr($_, 1);
				$current_chr =~ s/\s.*//g;
				if ($trans_file ne "") {
					if (exists $trans{$current_chr}) { $current_chr = $trans{$current_chr}; }
					else { die "Could not find chromosome '$current_chr' of genome file '$refseq_file' in chromosome translation file '$trans_file'! Aborted\n"; }
				}

				if($current_chr eq $chromosome) {
					while(<GENOME>) {
						chomp;
						if(substr($_, 0, 1)  eq ">") { 
							last; 
						}
						else {
							$chr_seq .= $_;
						}
					}
				}
			}

			if($chr_seq ne "" ) { last; }

		}
		close GENOME;

		if ($chr_seq eq "") {
			print STDERR "\n***ERROR***\nSequence for chromosome '$chromosome' not found in genome fasta file.\nCheck the conformity of the chromosome names of the different files!\n***********\n";
			exit;
		}

		### Start iterating over whole seq features and concomitantly over all SNPs and Indels (of current chr)
		foreach my $startpos (sort{$a<=>$b} keys %feat_by_pos) {
			my ($gene_name, $endpos) = split "#", $feat_by_pos{$startpos};

			### Create GeneSNPlist for this gene
			my $gene = new GeneSNPlist($gene_name, $startpos, $endpos, \%coding_subfeatures, \%{$seq_type{$gene_name}});

			### Get SNPs for this gene and store them in gene's SNPlist
			my $last_i = $var_array_pos;
			for (my $i = $var_array_pos; $i != $var_maxnr; $i++) {

				if ($vars->{vars}[$i]->{position} < $startpos && $vars->{vars}[$i]->{end} < $startpos) {
					$var_array_pos++ if ($last_i == $i-1);
					next;
				}

				$last_i = $i;

				# insertion is annotated as the preceding base in genome, if it is endpos -> it's outside gene
				last if ( ($vars->{vars}[$i]->{vartype} eq "I" && $vars->{vars}[$i]->{position} >= $endpos) 
					|| $vars->{vars}[$i]->{position} > $endpos );

				my $var;
				if ($vars->{vars}[$i]->{stype} ne "") {
					# duplicate SNP (more than one feature cover this SNP

					if ($vars->{vars}[$i]->{vartype} eq "S") {
						$var = new SNP();
						$var->init( $vars->{vars}[$i]->{ecotype}, $vars->{vars}[$i]->{chromosome},
							$vars->{vars}[$i]->{position}, $vars->{vars}[$i]->{ref_base},
							$vars->{vars}[$i]->{new_base}, 
							"S", $peak_pos, \@{$vars->{vars}[$i]->{fields}} );
						push @{$vars->{vars}}, $var;
					}
					else {
						$var = new Indel();
						$var->init( $vars->{vars}[$i]->{ecotype}, $vars->{vars}[$i]->{chromosome},
							$vars->{vars}[$i]->{global_start}, $vars->{vars}[$i]->{global_end},
							max($vars->{vars}[$i]->{position}, $startpos), min($vars->{vars}[$i]->{end}, $endpos),
							$vars->{vars}[$i]->{ref_base}, $vars->{vars}[$i]->{new_base}, $vars->{vars}[$i]->{seq}, 
							$vars->{vars}[$i]->{vartype}, $peak_pos, \@{$vars->{vars}[$i]->{fields}} );
						push @{$vars->{vars}}, $var;
					}
				}
				else {
					$var = $vars->{vars}[$i];
					if ($var->{vartype} eq "D") {
						$var->{position} = max($var->{position}, $startpos);
						$var->{end} = min($var->{end}, $endpos);
					}
				}

				$var->{gene_id} = $gene_name;
				if (!exists $coding_ann{$gene_name}) {
					if (!exists $seq_type{$gene_name}{$startpos}{"type"}) {
						$var->{stype} = $gene_type{$gene_name};
					}
					else {
						$var->{stype} = $seq_type{$gene_name}{$startpos}{"type"};
					}
				}

				push @{$gene->{Varlist}}, $var;
			}

			next if (scalar(@{$gene->{Varlist}}) == 0);

			$genes{$gene_name} = $gene;


			# for each coding feature:
			if (exists $coding_ann{$gene_name}) {

				### Extract gene seq and CDS seqs
				$gene_ann{$gene_name}[3] = substr($chr_seq,
								  $gene_ann{$gene_name}[0] - 1,
								  $gene_ann{$gene_name}[1] - $gene_ann{$gene_name}[0] + 1);

				foreach my $cds_start ( sort{$a<=>$b} keys %{$coding_ann{$gene_name}} ) {
					$coding_ann{$gene_name}{$cds_start}[3] = 
						substr(	$chr_seq,
							$cds_start - 1,
							$coding_ann{$gene_name}{$cds_start}[1] - $cds_start + 1);
				}

				if (exists $coding_ann{$gene_name}) {
					$gene->get_cds_seq( $chromosome, $gene_ann{$gene_name}[2], 1,
								$gene_ann{$gene_name}[3], %{$coding_ann{$gene_name}});

					### Annotate SNPs and fill SNPs in GeneSNPlist:
					$gene->get_protein_changes($gene_features{ $gene_type{$gene_name} }{$coding_feature{ $gene_type{$gene_name} }});
				}

			}

		} # for each feature

	} # if gff and genome file are specified



	### Print SNP results

	my $output_string="";
	my $oldpos = 0;

	if (-e $snp_file && -z $snp_file) { print STDERR "$snp_file is empty!\n"; }
	if ($mut) { print SNPOUT "chr\tstart\tend\tsample\ttype\tBPchange/sequence/length\tgene\tgene_pos\tcds_pos\tAAchange\tdegeneracy\n"; }

	foreach my $var (sort {$a->{peak_distance} <=> $b->{peak_distance}} @{$vars->{vars}} ) {

		$var->{stype} = $var->{stype} eq ""? "intergenic" : $var->{stype};

		$output_string = "";
		my $gene_pos = 0;
		my $cds_pos = 0;
		if ($var->{gene_id} ne "NA") {
			if ($var->{vartype} eq "S") {
				$gene_pos = ($genes{$var->{gene_id}}->{orientation} eq "+") ? $var->{position}-$genes{$var->{gene_id}}->{start}+1 : $genes{$var->{gene_id}}->{end}-$var->{position}+1;
				$cds_pos = ($genes{$var->{gene_id}}->{orientation} eq "+") ? $var->{cds_pos} : $genes{$var->{gene_id}}->{cds_length}-$var->{cds_pos}+1 if ($var->{cds_pos} != 0);
			}
			elsif ($var->{vartype} eq "I") {
				$gene_pos = ($genes{$var->{gene_id}}->{orientation} eq "+") ? $var->{global_start}-$genes{$var->{gene_id}}->{start}+1 : $genes{$var->{gene_id}}->{end}-$var->{global_end};
				$cds_pos = ($genes{$var->{gene_id}}->{orientation} eq "+") ? $var->{cds_pos} : $genes{$var->{gene_id}}->{cds_length}-$var->{cds_pos} if ($var->{cds_pos} != 0);
			}
			else {
				$gene_pos = ($genes{$var->{gene_id}}->{orientation} eq "+") ? max($var->{global_start}, $genes{$var->{gene_id}}->{start})-$genes{$var->{gene_id}}->{start}+1 : $genes{$var->{gene_id}}->{end} - min($var->{global_end}, $genes{$var->{gene_id}}->{end}) + 1;
				if ($var->{cds_pos_left} > 0) {
					my $sep = $mut ? "-" : "\t";
					$cds_pos = $var->{cds_pos_left} . $sep . $var->{cds_pos_right} if ($genes{$var->{gene_id}}->{orientation} eq "+");
					$cds_pos = ($genes{$var->{gene_id}}->{cds_length} - $var->{cds_pos_right} + 1) . $sep . ($genes{$var->{gene_id}}->{cds_length} - $var->{cds_pos_left} + 1) if ($genes{$var->{gene_id}}->{orientation} eq "-");
				}
			}
		}


		if (!$mut) {
			##### Regular output: #####

                        $output_string .= ($trans_file ne "" ? $trans_rev{$var->{chromosome}} : $var->{chromosome}) . "\t";
                        $output_string .= $var->{position} . "\t";
			$output_string .= $var->{ecotype} . "\t";
			$output_string .= $var->{ref_base} . "\t" . $var->{new_base} . "\t";

			if ($var->{vartype} eq "S") {
				#$output_string .= $var->{ecotype} . "\t";
				#$output_string .= ($trans_file ne "" ? $trans_rev{$var->{chromosome}} : $var->{chromosome}) . "\t";
				#$output_string .= $var->{position} . "\t" . $var->{ref_base} . "\t" . $var->{new_base} . "\t";
		
				foreach (@{$var->{fields}}) {
					$output_string .= $_ . "\t";
				}

				if( ($gff ne "") && ($refseq_file ne "") ) {
					$output_string .= $var->{stype};

					if($var->{gene_id} ne "NA") {
						$output_string .= "\t" . $var->{gene_id} . "\t1\t";
						$output_string .= $genes{$var->{gene_id}}->{start} . "\t" . $genes{$var->{gene_id}}->{end} . "\t" . $genes{$var->{gene_id}}->{orientation};
						$gene_pos = ($genes{$var->{gene_id}}->{orientation} eq "+") ? $var->{position}-$genes{$var->{gene_id}}->{start}+1 : $genes{$var->{gene_id}}->{end}-$var->{position}+1;
						$output_string .= "\t" . $gene_pos;
					}
		
					if($var->{cds_pos} != 0) {
						$cds_pos = ($genes{$var->{gene_id}}->{orientation} eq "+") ? $var->{cds_pos} : $genes{$var->{gene_id}}->{cds_length}-$var->{cds_pos}+1;
						my $syn_nonsyn = "Syn";
						if($var->{ref_aa} ne $var->{new_aa}) { $syn_nonsyn = "Nonsyn"; }
						$output_string .= "\t" . $genes{$var->{gene_id}}->{cds_length} . "\t" .
								$cds_pos . "\t" . 
								$var->{codon_pos} . "\t$syn_nonsyn\t" .
								$var->{ref_aa} . "\t" . $var->{new_aa} . "\t" .
								$var->{codon} . "\t" . $var->{degeneracy};
					}
					elsif ($var->{splicechange}) {
						$output_string .= "\t" . $var->{codon_pos} . "\t" . $var->{codon};
					}
				}
			}
			elsif ($var->{vartype} eq "I") {
				#$output_string .= $var->{ecotype} . "\t";

				#$output_string .= ($trans_file ne "" ? $trans_rev{$var->{chromosome}} : $var->{chromosome}) . "\t";
				#$output_string .= $var->{position} . "\t" . $var->{end} . "\t" . $var->{seq} . "\t";
				
		
                                foreach (@{$var->{fields}}) {
                                        $output_string .= $_ . "\t";
                                }

				if ( ($gff ne "") && ($refseq_file ne "") ) {
					$output_string .= $var->{stype};
					if ($var->{gene_id} ne "NA") {
						$output_string .= "\t" . $var->{gene_id} . "\t1\t";
						$output_string .= $genes{$var->{gene_id}}->{start} ."\t" . $genes{$var->{gene_id}}->{end} . "\t" . $genes{$var->{gene_id}}->{orientation} . "\t";
						$gene_pos = ($genes{$var->{gene_id}}->{orientation} eq "+") ? $var->{global_start}-$genes{$var->{gene_id}}->{start}+1 : $genes{$var->{gene_id}}->{end}-$var->{global_end};
						$output_string .= $gene_pos . "\t";
					}
					if ($var->{cds_pos} != 0) {
						$cds_pos = ($genes{$var->{gene_id}}->{orientation} eq "+") ? $var->{cds_pos} : $genes{$var->{gene_id}}->{cds_length}-$var->{cds_pos};
						$output_string .= $genes{$var->{gene_id}}->{cds_length} . "\t" . $cds_pos;
					}
				}
			}
			elsif ($var->{vartype} eq "D") {
				#$output_string .= $var->{ecotype} . "\t";

				#$output_string .= ($trans_file ne "" ? $trans_rev{$var->{chromosome}} : $var->{chromosome}) . "\t";
				#$output_string .= $var->{global_start} . "\t" . $var->{global_end} . "\t" . ($var->{global_end}-$var->{global_start}+1) . "\t";

                                foreach (@{$var->{fields}}) {
                                        $output_string .= $_ . "\t";
                                }

				if ( ($gff ne "") && ($refseq_file ne "") ) {
					$output_string .= $var->{stype};

					if ($var->{gene_id} ne "NA") {
						$output_string .= "\t" . $var->{gene_id} . "\t1\t";
						$output_string .= $genes{$var->{gene_id}}->{start} . "\t" . $genes{$var->{gene_id}}->{end} . "\t" . $genes{$var->{gene_id}}->{orientation} . "\t";
						$gene_pos = ($genes{$var->{gene_id}}->{orientation} eq "+") ? max($var->{global_start}, $genes{$var->{gene_id}}->{start})-$genes{$var->{gene_id}}->{start}+1 : $genes{$var->{gene_id}}->{end} - min($var->{global_end}, $genes{$var->{gene_id}}->{end}) + 1;
						$output_string .= $gene_pos;
					}

					if ($var->{cds_pos_left} > 0) {
						$output_string .= "\t" . $genes{$var->{gene_id}}->{cds_length};
						$output_string .= "\t" . $var->{cds_pos_left} . "\t" . $var->{cds_pos_right} if ($genes{$var->{gene_id}}->{orientation} eq "+");
						$output_string .= "\t" . ($genes{$var->{gene_id}}->{cds_length} - $var->{cds_pos_right} + 1) . "\t" . ($genes{$var->{gene_id}}->{cds_length} - $var->{cds_pos_left} + 1) if ($genes{$var->{gene_id}}->{orientation} eq "-");

						my $cds_len = $var->{cds_pos_right}-$var->{cds_pos_left}+1;
						$output_string .= "\t" . $cds_len;
					}
				}
			}
		}
		else {
			##### .mut output #####

			$output_string .= ($trans_file ne "" ? $trans_rev{$var->{chromosome}} : $var->{chromosome}) . "\t";
			if ($var->{vartype} eq "D") { $output_string .= $var->{global_start} . "\t" . $var->{global_end} . "\t" . $var->{ecotype} . "\t"; }
			if ($var->{vartype} eq "I") { $output_string .= $var->{position} . "\t" . $var->{end} . "\t" . $var->{ecotype} . "\t"; }
			if ($var->{vartype} eq "S") { $output_string .= $var->{position} . "\t" . $var->{position} . "\t" . $var->{ecotype} . "\t"; }
			
			if ($var->{new_stop}>0)		{ $output_string .= "premature_stop (".$var->{new_stop}.")"; }
			elsif ($var->{lost_stop})	{ $output_string .= "lost_stop"; }
			elsif ($var->{new_start})	{ $output_string .= "new_start"; }
			elsif ($var->{lost_start})	{ $output_string .= "lost_start"; }
			elsif ($var->{splicechange})	{ $output_string .= "splicechange"; }
			elsif ($var->{vartype} eq "S" && $var->{ref_aa} ne $var->{new_aa})	{ $output_string .= "synonymous"; }
			elsif ($var->{vartype} eq "S" && $var->{ref_aa} ne $var->{new_aa})	{ $output_string .= "nonsynonymous"; }
			else				{ $output_string .= $var->{stype}; }

			if ($var->{vartype} eq "I") {
				$output_string .= "\t+" . $var->{seq};
			}
			elsif ($var->{vartype} eq "D") {
				$output_string .= "\t-" . $var->{seq} if ($var->{seq} ne "");	# seq is available in vcf files
				$output_string .= "\t" . ($var->{cds_pos_right}-$var->{cds_pos_left}+1) if ($var->{seq} eq "");	# seq is NOT available in shore files
			}
			else {
				$output_string .= "\t" . $var->{ref_base}.">".$var->{new_base};
			}

			if ($var->{gene_id} ne "NA") {
				$output_string .=  "\t" . $var->{gene_id} . "\t" . $gene_pos;
				if ($cds_pos ne "0") {
					$output_string .= "\t" . $cds_pos;
					if ($var->{vartype} eq "S") {
						$output_string .= "\t" . $var->{codon}.">".$var->{codon_alt} . "\t" . $var->{ref_aa}.">".$var->{new_aa} . "\t" . $var->{degeneracy}
					}
				}
			}
		}

		if (!$mut) {
			my $feat = "";
			if ($var->{new_stop}||$var->{lost_stop}||$var->{lost_start}||$var->{new_start}||$var->{splicechange}||$var->{frameshift}) {
				if ($var->{new_stop}>0)		{ if ($feat ne "") { $feat.=","; } $feat .= "new_stop (".$var->{new_stop}.")"; }
				if ($var->{lost_stop}) 		{ if ($feat ne "") { $feat.=","; } $feat .= "lost_stop"; }
				if ($var->{new_start}) 		{ if ($feat ne "") { $feat.=","; } $feat .= "new_start"; }
				if ($var->{lost_start})		{ if ($feat ne "") { $feat.=","; } $feat .= "lost_start"; }
				if ($var->{splicechange})	{ if ($feat ne "") { $feat.=","; } $feat .= "splicechange"; }
				if ($var->{frameshift})	 	{ if ($feat ne "") { $feat.=","; } $feat .= "frameshift"; }
			}
	
			$output_string .= "\t$feat" if ($feat ne "");
		}

		$output_string .= "\n";


		my $FH = *SNPOUT;

	#	if (-z $snp_file && $mut) {
	#		print $FH "chr\tstart\tend\tsample\ttype\tBPchange/sequence/length\tgene\tgene_pos\tcds_pos\tAAchange\tdegeneracy\n";
	#	}
		print $FH $output_string;

		$oldpos = $var->{position};

	}

	# Print gene alignments in file GENEOUT;
	if ($cds_alignments) {
		foreach (sort keys %genes) {
			if (scalar(@{$genes{$_}->{coding_Vars}}) > 0) {
				print GENEOUT ">".$_."\n";
				print GENEOUT $genes{$_}->{align_ref}."\n";
				print GENEOUT $genes{$_}->{align_alt}."\n";
			}
		}
	}

	$next_chr = 0 if ($target_region);

} # for each chromosome

close SNPOUT if ($snp_file ne "");
close GENEOUT if ($cds_alignments);

#system("rm $var_file");


exit(0);


### Read command line parameters
sub GetCom {
  my @usage = ("\nUsage: $0

Mandatory: 
--snp      STRING      Input rbci file containing variants
--out	   STRING      Output file prefix (Output will be in <prefix>.rbca)

If chromosome names are not numerical, please provide a chromosome translation file:
--trans    STRING      Chromosome translation table (column1: numerical chrom.name,
		       column2: chrom.name as in all provided files)

Target region (optional, only used when all three options are specified):
--chrom    STRING      Chromosome of target region
--start    INT         Start of target region
--end      INT         End of target region

Functional SNP and indel annotation (only used if genome and gff are specified):
--genome   STRING      Reference sequence file (chromosome names have to be equal to SNP file)
--gff      STRING      Gene annotation in GFF format
--config   STRING      Config file for specific gff file format (default: O.sativa Genbank)

Output option:
--mut		       Output format is .mut (used for IGV)
--alleles  STRING      Output only 'major', 'minor' or 'all' alleles (defautl: all)
--filter   STRING[,..] List of filters in FILTER column of rbci file leading to exclusion
--cds_alignments       Outputs coding sequence alignments of reference and variant genes
		       (writes into <prefix>.genes)
\n");

#--all_isoforms         Report annotation for all listed isoforms (default: only isoform 1)

	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "snp=s", "chrom=s", "start=s", "end=s", "genome=s", "gff=s", "config=s", "iso", "out=s", "trans=s", "cds_alignments", "mut", "alleles=s", "filter=s");

	die("Please specify an rbci file and an output prefix\n") unless (defined($CMD{snp}) || defined($CMD{out}));

	$chromosome = 0;
	$start = 0;
	$end = 0;
	$gff = "";
	$refseq_file = "";
	$refseq_error_file = "";
	$snp_file = "";

	if (defined $CMD{chrom} || defined $CMD{start} || defined $CMD{end}) {
		die("Please specify chromosome of target region\n") unless defined($CMD{chrom});
		die("Please specify start of target region\n") unless defined($CMD{start});
		die("Please specify end of target region\n") unless defined($CMD{end});

		$target_region = 1;
	}

	if (defined $CMD{genome} || defined $CMD{gff}) {
		die("Please specify genome fasta file\n") unless defined($CMD{genome});
		die("Please specify gff file\n") unless defined($CMD{gff});
	}

	if(defined $CMD{chrom}) { $chromosome = $CMD{chrom}; }
	if(defined $CMD{start}) { $start = $CMD{start}; }
	if(defined $CMD{end}) { $end = $CMD{end}; }
	if(defined $CMD{gff}) { $gff = $CMD{gff}; }
	if(defined $CMD{genome}) { $refseq_file = $CMD{genome}; }
	if(defined $CMD{config}) { $config_file = $CMD{config}; }
	if(defined $CMD{iso}) { $all_isoforms = 1; }
	if(defined $CMD{snp}) { $snp_file = $CMD{snp}; }
	if(defined $CMD{trans}) { $trans_file = $CMD{trans}; }
	if(defined $CMD{out}) { $outfile_prefix = $CMD{out}; }
	if(defined $CMD{cds_alignments}) { $cds_alignments = 1; }
	if(defined $CMD{mut}) { $mut = 1; }

	if(defined $CMD{alleles}) {
		if ($CMD{alleles} ne "major" && $CMD{alleles} ne "minor" && $CMD{alleles} ne "all") {
			die "Values for option --alleles are only 'major','minor' or 'all'!\n";
		}
		else {
			$majmin_allele = $CMD{alleles};
		}
	}

	if(defined $CMD{filter}) {
		@filter = split/,/, $CMD{filter};
	}
}

