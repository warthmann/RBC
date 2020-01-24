#!/usr/bin/perl 

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;
#use Bio::Perl;
#use SNPlist;
use Varlist;

package GeneSNPlist;

use List::Util qw[min max];
sub new {
	my ($self, $name, $start, $end, $coding_features, $gene_features) = @_;
	$self = {
		coding_features=> {},
		gene_id        => $name,
		gene_features  => {},
		isoform        => 1,
		chromosome     => '',
		start          => $start,
		end            => $end,
		len            => ($end-$start+1),
		cds_length     => 0,
		protein_length => 0,
		orientation    => '',
		Varlist	       => [],
#		SNP_lists      => {},
#		Ins_lists      => {},
#		Del_lists      => {},
		coding_Vars    => [],
#		NS_changes     => {},
#		AA_changes     => {},
		CDSexon        => {},
		protein        => {},
		gene_seq_ref   => '',
		align_ref      => '',
		align_alt      => '',
		coding_seq_AllVar => '',
		debug	       => 0,
	};

	$self->{coding_features} = $coding_features;
	$self->{gene_features} = $gene_features;

	bless $self;
	return $self;
}

#############################################################################################
### Fill Reference coding sequence underlying gene
sub get_cds_seq {
#	my ($self, $chromosome, $gene_start, $gene_end, $orientation, $gene_id, $isoform, $gene_seq_ref, %cds_ann) = @_;
	my ($self, $chromosome, $orientation, $isoform, $gene_seq_ref, %cds_ann) = @_;
	$self->{chromosome}  = $chromosome;
#	$self->{start}       = $gene_start;
#	$self->{end}         = $gene_end;
#	$self->{len}	     = $gene_end-$gene_start+1;
	$self->{orientation} = $orientation;
#	$self->{gene_id}     = $gene_id;
	$self->{isoform}     = $isoform;
	$self->{gene_seq_ref}= $gene_seq_ref;

	# Get features of gene (range, CDS, RNA)
	foreach my $cds_start ( sort {$a <=> $b} keys %cds_ann ) {
		my $cds_seq = $cds_ann{$cds_start}[3];
		my $cds_end = $cds_ann{$cds_start}[1];
		
		$self->{CDSexon}{$cds_start} = $cds_end;
		$self->{cds_length} = $self->{cds_length} + $cds_end - $cds_start + 1;
		$self->{coding_seq_ref} .= $cds_seq;
	}

	# reverse complement if orientation = "-"
	#if($self->{orientation} eq "-") {
	#	$self->{coding_seq_ref} = reverse($self->{coding_seq_ref});
	#	$self->{coding_seq_ref} =~ tr/ACGTUacgtu/TGCAATGCAA/;
	#}

	$self->{coding_seq_SNP} = $self->{coding_seq_ref};

	# Protein length
	$self->{protein_length} = $self->{cds_length} / 3;

	return($self->{protein_length});
}


my %DEGENERACY = (
"AAA" => "002","AAC" => "002","AAG" => "002","AAT" => "002","ACA" => "004","ACC" => "004","ACG" => "004","ACT" => "004","AGA" => "202","AGC" => "002","AGG" => "202","AGT" => "002","ATA" => "002","ATC" => "002","ATG" => "000","ATT" => "002","CAA" => "002","CAC" => "002","CAG" => "002","CAT" => "002","CCA" => "004","CCC" => "004","CCG" => "004","CCT" => "004","CGA" => "204","CGC" => "004","CGG" => "204","CGT" => "004","CTA" => "204","CTC" => "004","CTG" => "204","CTT" => "004","GAA" => "002","GAC" => "002","GAG" => "002","GAT" => "002","GCA" => "004","GCC" => "004","GCG" => "004","GCT" => "004","GGA" => "004","GGC" => "004","GGG" => "004","GGT" => "004","GTA" => "004","GTC" => "004","GTG" => "004","GTT" => "004","TAA" => "022","TAC" => "002","TAG" => "022","TAT" => "002","TCA" => "004","TCC" => "004","TCG" => "004","TCT" => "004","TGA" => "020","TGC" => "002","TGG" => "000","TGT" => "002","TTA" => "202","TTC" => "002","TTG" => "202","TTT" => "002",
);


my %AA = (
"AAA" => "K","AAC" => "N","AAG" => "K","AAT" => "N","ACA" => "T","ACC" => "T","ACG" => "T","ACT" => "T","AGA" => "R","AGC" => "S","AGG" => "R","AGT" => "S","ATA" => "I","ATC" => "I","ATG" => "M","ATT" => "I","CAA" => "Q","CAC" => "H","CAG" => "Q","CAT" => "H","CCA" => "P","CCC" => "P","CCG" => "P","CCT" => "P","CGA" => "R","CGC" => "R","CGG" => "R","CGT" => "R","CTA" => "L","CTC" => "L","CTG" => "L","CTT" => "L","GAA" => "E","GAC" => "D","GAG" => "E","GAT" => "D","GCA" => "A","GCC" => "A","GCG" => "A","GCT" => "A","GGA" => "G","GGC" => "G","GGG" => "G","GGT" => "G","GTA" => "V","GTC" => "V","GTG" => "V","GTT" => "V","TAA" => "*","TAC" => "Y","TAG" => "*","TAT" => "Y","TCA" => "S","TCC" => "S","TCG" => "S","TCT" => "S","TGA" => "*","TGC" => "C","TGG" => "W","TGT" => "C","TTA" => "L","TTC" => "F","TTG" => "L","TTT" => "F",
"NNN" => "X", "ANN" => "X", "CNN" => "X", "GNN" => "X", "TNN" => "X", "NAN" => "X", "NCN" => "X", "NGN" => "X", "NTN" => "X", "NNA" => "X", "NNC" => "X", "NNG" => "X", "NNT" => "X", "AAN" => "X", "ACN" => "X", "AGN" => "X", "ATN" => "X", "CAN" => "X", "CCN" => "X", "CGN" => "X", "CTN" => "X", "GAN" => "X", "GCN" => "X", "GGN" => "X", "GTN" => "X", "TAN" => "X", "TCN" => "X", "TGN" => "X", "TTN" => "X", "ANA" => "X", "ANC" => "X", "ANG" => "X", "ANT" => "X", "CNA" => "X", "CNC" => "X", "CNG" => "X", "CNT" => "X", "GNA" => "X", "GNC" => "X", "GNG" => "X", "GNT" => "X", "TNA" => "X", "TNC" => "X", "TNG" => "X", "TNT" => "X", "NAA" => "X", "NAC" => "X", "NAG" => "X", "NAT" => "X", "NCA" => "X", "NCC" => "X", "NCG" => "X", "NCT" => "X", "NGA" => "X", "NGC" => "X", "NGG" => "X", "NGT" => "X", "NTA" => "X", "NTC" => "X", "NTG" => "X", "NTT" => "X",
);


#########################################################################################
### Calculates CDS and protein changes for all SNPs in "SNP_lists" and fills SNP object
sub get_protein_changes {
	my ( $self, $cds_name ) = @_;
#print $self->{gene_id}." ".$self->{ref_coding}."\n";
#exit if ($self->{ref_coding} eq "");

	my %last_feature=();
	my $last_feat_start = 0;
	my $varlist_idx = 0;
	my $cds_pos = 1;
	my $genomic_pos = 1;
	my $coding_indels = 0;
	my $last_var = 0;

	my $var = $self->{Varlist}[$varlist_idx];

print "@@@@@@@@@@@@@\nVAR: ".$var->{vartype}." @ ".$var->{position}."\n@@@@@@@@@@@@@\n" if ($self->{debug});
print "type: ".$var->{stype}."\n" if ($self->{debug});

	foreach my $feat_start (sort{$a<=>$b} keys %{$self->{gene_features}}) {
		my %feature = %{$self->{gene_features}{$feat_start}};
		$genomic_pos = $feat_start;
print "\n--------------\nFEATURE (".$self->{orientation}.") at pos ".$feat_start."-".$feature{"end"}."\n"."cds_pos $cds_pos genomic_pos $genomic_pos\n--------------\n\n" if ($self->{debug});

		# variation is before very first feature
		while ($last_feat_start == 0 && $var->{end} < $feat_start+($var->{vartype} eq "I")) {
			# next variation
			$varlist_idx++;
			if (!defined $self->{Varlist}[$varlist_idx]) {
				$last_var = 1;
				last;
			}
			$var = $self->{Varlist}[$varlist_idx];
print "@@@@@@@@@@@@@@\nVAR: ".$var->{vartype}." @ ".$var->{position}."\n@@@@@@@@@@@@@@\n" if ($self->{debug});
		}


		# variation in introns:
		while (!$last_var && exists $last_feature{"end"} && $var->{position} > $last_feature{"end"} && $var->{position} < $feat_start) {
print "> between gene features (intronic)\n" if ($self->{debug});
			#$var->{gene_id} = $self->{gene_id};

			# DELETION:
			if ($var->{vartype} eq "D") {
				if ($var->{global_start} != $var->{position} && $var->{stype} ne "") {
					# variation started earlier

					# if variation started in a coding feature (cds_pos_left!=0) and last feature was a coding feature -> splice site disrupted
					$var->{splicechange} = 1 if ($var->{cds_pos_left} != 0 && exists $self->{coding_features}{$last_feature{"type"}});
					$var->{stype} .= "-";
print "type: ".$var->{stype}."\n" if ($self->{debug});
				}
				else {
					# variation starts in intron
					#$var->{gene_pos} = ($self->{orientation} eq "+") ? $var->{position} - $self->{start} + 1 : $self->{end} - $var->{end} + 1;
					#if ($var->{gene_pos} == 0) {
					#	if ($var->{global_start} == $var->{position} && $self->{orientation} eq "+") {
					#		$var->{gene_pos} = $var->{position} - $self->{start} + 1;
					#	}
					#	elsif ($var->{end} <= $feature{"end"} && $self->{orientation} eq "-") {
					#		$var->{gene_pos} = $var->{end} - $self->{start} + 1;
					#	}
					#}
				}
					
				$var->{stype} .= "intronic/noncoding";
				$var->{stype} .= ",".(min($var->{end},$feat_start-1)-$var->{position}+1) if ($var->{end} >= $feat_start || $var->{global_start} != $var->{position});
				#$var->{stype} .= ($feat_start - $var->{position} + 1) if ($var->{end} >= $feat_start || $var->{global_start} != $var->{position});
print "type: ".$var->{stype}." EEEEEEEEND ".$var->{end}." FEATURE_START $feat_start\n" if ($self->{debug});

				# check if deletion spans breakpoint between intron and other features:
				if ($var->{end} >= $feat_start) {
					$var->{splicechange} = 1 if (%last_feature && exists $self->{coding_features}{$last_feature{"type"}} && exists $self->{coding_features}{$feature{"type"}});
					# update deletion position to the start of the following feature
					# (original position is in $var->{global_start}):
					$var->{position} = $feat_start;
				}
				else {
				# if not, deletion ends in intron
					# next variation:
					$varlist_idx++;
					if (!defined $self->{Varlist}[$varlist_idx]) {
						$last_var = 1;
						last;
					}
					$var = $self->{Varlist}[$varlist_idx];
print "@@@@@@@@@@@@@@\nVAR: ".$var->{vartype}." @ ".$var->{position}."\n@@@@@@@@@@@@@@\n" if ($self->{debug});
				}
print "var->position ".$var->{position}."\n" if ($self->{debug});
			}
			else {
			# SNP or INSERTION:
				$var->{stype} = "intronic/noncoding";
print "type: ".$var->{stype}."\n" if ($self->{debug});
				#$var->{gene_pos} = ($self->{orientation} eq "+") ? $var->{position} - $self->{start} + 1 : $self->{end}-$self->{start}+1 - $var->{position} + 1;

				# check if splice site is hit:
				if (
					# first 2 bases after coding feature or last 2 bases before coding feature:
					($var->{position} <= $last_feature{"end"}+1+($var->{vartype} eq "S") || $var->{position} >= $feat_start-2)
					&&
					# between coding features:
					(%last_feature && exists $self->{coding_features}{$last_feature{"type"}} && exists $self->{coding_features}{$feature{"type"}})
				) {
					$var->{splicechange} = 1;
					$var->{stype} = "splicesite";
print "type: ".$var->{stype}."\n" if ($self->{debug});

					# in case of SNP: parse splice motif
					if ($var->{vartype} eq "S") {
						if ($var->{position} <= $last_feature{"end"}+2) {
							$var->{codon} = substr( $self->{gene_seq_ref}, $last_feature{"end"}-$self->{start}+1, 2 );
							$var->{codon_pos} = $var->{position} - $last_feature{"end"};
						}
						else {
							$var->{codon} = substr( $self->{gene_seq_ref}, $feat_start-3-$self->{start}+1, 2 );
							$var->{codon_pos} = $feat_start - $var->{position} == 1 ? 2 : 1;
						}
						if ($self->{orientation} eq "-") {
							$var->{codon} = reverse($var->{codon});
							$var->{codon} =~ tr/ACGTUacgtu/TGCAATGCAA/;
							$var->{codon_pos} = $var->{codon_pos} == 1 ? 2 : 1;
						}
						substr($var->{codon}, $var->{codon_pos}%2, 1, lc(substr($var->{codon}, 0, 1)));
						$var->{codon_alt} = $var->{codon};
						substr($var->{codon_alt}, $var->{codon_pos}-1, 1, $var->{new_base});
					}
				}

				# next variation:
				$varlist_idx++;
				if (!defined $self->{Varlist}[$varlist_idx]) {
					$last_var = 1;
					last;
				}
				$var = $self->{Varlist}[$varlist_idx];
print "@@@@@@@@@@@@@@\nVAR: ".$var->{vartype}." @ ".$var->{position}."\n@@@@@@@@@@@@@@\n" if ($self->{debug});

			} # SNP / Insertion

		} # for ech variation in intron

		# variation in gene feature:
		while (!$last_var && $var->{position} >= $feat_start && $var->{position} <= $feature{"end"}) {
print "> in a gene feature\n" if ($self->{debug});

			$var->{gene_id} = $self->{gene_id};

			# in coding feature:
			if (exists $self->{coding_features}{$feature{"type"}}) {
print "> > coding\n" if ($self->{debug});

				push @{$self->{coding_Vars}}, $var;

				# SNP or INSERTION
				if ($var->{vartype} ne "D") {
					
					#my $cds_loc = ($self->{orientation} eq "+") ? $cds_pos + $var->{position} - $genomic_pos : $self->{cds_length}-($cds_pos + $var->{position}-$genomic_pos+1) + ($var->{vartype} eq "S");
					$var->{cds_pos} = $cds_pos + $var->{position} - $genomic_pos;
print "    cds_pos: ".$var->{cds_pos}."\n" if ($self->{debug});
#					$self->{coding_SNP}{$var->{$cds_pos}} = $var;

					my $seq = substr( $self->{coding_seq_ref}, $cds_pos-1, $var->{cds_pos}-$cds_pos+($var->{vartype} eq "I") );
					my $new_seq = $var->{vartype} eq "S" ? $var->{new_base} : $var->{seq};
					$self->{align_alt} .= $seq . $new_seq;
					$self->{coding_seq_AllVar} .= $seq . $new_seq;

					if ($var->{vartype} eq "I") {
						$self->{align_ref} .= $seq . ("-" x ( length($new_seq) ));
						$var->{alt_coding} = substr( $self->{coding_seq_ref}, 0, $var->{cds_pos} ) . $new_seq . substr( $self->{coding_seq_ref}, $var->{cds_pos}, $self->{cds_length}-$var->{cds_pos}+1 );
						#if ($self->orientation eq "-") {
						#	$var->{alt_coding} = reverse $var->{alt_coding};
						#	$var->{alt_coding} =~ tr/ACGTacgt/TGCAtgca/;
						#}
						#compare_proteins($self->orientation eq "-" ? $coding_seq_ref_reverse : $coding_ref_seq, $var->{alt_coding}, $var->{cds_pos});
						$self->compare_proteins( $var, $var->{cds_pos} );

						#$var->{cds_pos}++;

						#$coding_indels++;

						# will affect splice site if 1) varpos is last pos of feature, 2) coding feature is not last one
						$var->{splicechange} = 1 if ($var->{position} == $feature{"end"} && $var->{cds_pos} != $self->{cds_length});
						$var->{frameshift} = 1 if (length($new_seq) % 3 != 0);
						# check if insertion destroys start or stop codon
						if ($self->{orientation} eq "+") {
							$var->{lost_start} = 1 if ($var->{cds_pos} < 3);
							$var->{lost_stop}  = 1 if ($var->{cds_pos} != $self->{cds_length} && $var->{cds_pos} > $self->{cds_length}-3);
						}
						else {
							$var->{lost_stop}  = 1 if ($var->{cds_pos} < 3);
                                                        $var->{lost_start} = 1 if ($var->{cds_pos} != $self->{cds_length} && $var->{cds_pos} > $self->{cds_length}-3);
						}
					}
					else {
						$self->{align_ref} .= $seq . substr( $self->{coding_seq_ref}, $var->{cds_pos}-1, 1);
						#substr($self->{coding_seq_SNP}, $var->{cds_pos}-1, 1, $var->{new_base});

						#*****************************************************
						# Annotate Amino Acid:
						#*****************************************************
						$var->{codon} = substr( $self->{coding_seq_ref}, $var->{cds_pos} - ( ($var->{cds_pos}-1) % 3 ) - 1, 3 );
						$var->{codon_pos} = ($var->{cds_pos}-1) % 3 + 1;

						my $base = $new_seq;

#						if ($self->{orientation} eq "+") {
#							$var->{codon} = substr( $self->{coding_seq_ref}, $var->{cds_pos} - ( ($var->{cds_pos}-1) % 3 ) - 1, 3 );
#							$var->{codon_pos} = ($var->{cds_pos}-1) % 3 + 1;
#						}
#						else {
						if ($self->{orientation} eq "-") {
#							$var->{codon} = substr( $self->{coding_seq_ref}, ($self->{cds_length}-$var->{cds_pos}+1) - ( (($self->{cds_length}-$var->{cds_pos}+1)-1) % 3 ) - 1, 3 );
							$var->{codon} = reverse($var->{codon});
							$var->{codon} =~ tr/ACGTUacgtu/TGCAATGCAA/;
							$var->{codon_pos} = ((($self->{cds_length}-$var->{cds_pos}+1)-1) % 3 + 1);
							$base =~ tr/ACGTUacgtu/TGCAATGCAA/;
						}

						$var->{codon_alt} = $var->{codon};
						substr($var->{codon_alt}, $var->{codon_pos}-1, 1, $base);

						if (length($var->{codon}) < 3) {
							$var->{ref_aa} = "NNN";
						} else {
							$var->{ref_aa} = $AA{$var->{codon}};
						}
						$var->{new_aa} = ($var->{codon_alt} =~ m/[^ACGT]/) ? "X" : $AA{$var->{codon_alt}};

						$var->{degeneracy} = ($var->{codon_alt} =~ m/[^ACGT]/) ? "-" : substr($DEGENERACY{$var->{codon_alt}}, $var->{codon_pos}-1, 1);

						substr($var->{codon}, 0, 1, lc(substr($var->{codon}, 0, 1))) if ($var->{codon_pos} != 1);
						substr($var->{codon}, 1, 1, lc(substr($var->{codon}, 1, 1))) if ($var->{codon_pos} != 2);
						substr($var->{codon}, 2, 1, lc(substr($var->{codon}, 2, 1))) if ($var->{codon_pos} != 3);
						substr($var->{codon_alt}, 0, 1, lc(substr($var->{codon_alt}, 0, 1))) if ($var->{codon_pos} != 1);
						substr($var->{codon_alt}, 1, 1, lc(substr($var->{codon_alt}, 1, 1))) if ($var->{codon_pos} != 2);
						substr($var->{codon_alt}, 2, 1, lc(substr($var->{codon_alt}, 2, 1))) if ($var->{codon_pos} != 3);

						if ($var->{ref_aa} ne $var->{new_aa}) {
							if ($var->{new_aa} eq "*") { $var->{new_stop} = 1; }
							if ($var->{ref_aa} eq "*") { $var->{lost_stop} = 1; }
							if ($self->{orientation} eq "+" && $var->{cds_pos} <= 3 && $var->{new_aa} ne "M") { $var->{lost_start} = 1; }
							if ($self->{orientation} eq "-" && $var->{cds_pos} >= $self->{cds_length}-2 && $var->{new_aa} ne "M") { $var->{lost_start} = 1; }
						}
						#****************************************************
						#****************************************************
					}

					$cds_pos = $var->{cds_pos} + 1;
					$genomic_pos = $var->{position} + 1;
				}
				# DELETION
				else {

					if ($var->{cds_pos_left} == 0) {
						if ($var->{global_start} < $self->{start}) {
							# deletion started before this gene
							$var->{cds_pos_left} = 1;
						}
						else {
							# variation starts in this feature
							$var->{cds_pos_left} = $cds_pos + $var->{position} - $genomic_pos;

							# fill with sequence from last variation or feature start:
							if ($var->{position} - $genomic_pos > 0) {
								my $seq = substr( $self->{coding_seq_ref}, $cds_pos-1, $var->{position} - $genomic_pos );
								$self->{align_alt} .= $seq;
								$self->{coding_seq_AllVar} .= $seq;
							}
						}

						if ($self->{orientation} eq "+") {
							if ($var->{cds_pos_left} < 4) { $var->{lost_start} = 1; }
							if ($var->{cds_pos_left} >= $self->{cds_length}-3) { $var->{lost_stop} = 1; }
						}
						else {
							if ($var->{cds_pos_left} < 4) { $var->{lost_stop} = 1; }
							if ($var->{cds_pos_left} >= $self->{cds_length}-3) { $var->{lost_start} = 1; }
						}
					}

print "    cds_pos_left: ".$var->{cds_pos_left}."\n" if ($self->{debug});
					if ($var->{end} > $feature{"end"}) {
						# variation exceeds feature
						$var->{cds_pos_right} = $cds_pos + $feature{"end"}-$genomic_pos;
						#$cds_pos += $feature{"end"}-$genomic_pos+1;
						#$genomic_pos = $feature{"end"}+1;
					}
					else {
						# variation ends in feature
						$var->{cds_pos_right} = $cds_pos + $var->{end} - $genomic_pos;
						#$var->{alt_coding} = substr( $self->{coding_seq_ref}, $var->{cds_pos_right} );
						#$genomic_pos = $var->{"end"}+1;
						#$cds_pos = $var->{cds_pos_right}+1;
					}
print "    cds_pos_right: ".$var->{cds_pos_right}." ".($var->{cds_pos_right} - max($var->{cds_pos_left}, $cds_pos) + 1)."\n" if ($self->{debug});

					$self->{align_alt} .= "-" x ( $var->{cds_pos_right} - max($var->{cds_pos_left}, $cds_pos) + 1 );
					$self->{align_ref} .= substr( $self->{coding_seq_ref}, $cds_pos-1, $var->{cds_pos_right}-$cds_pos+1 );

					$genomic_pos = min($var->{"end"}, $feature{"end"}) + 1;
					$cds_pos = $var->{cds_pos_right} + 1;

					if ($self->{orientation} eq "+") {
						if ($var->{cds_pos_right} >= $self->{cds_length}-3) { $var->{lost_stop} = 1; }
						if ($var->{cds_pos_right} < 4) { $var->{lost_start} = 1; }
					}
					else {
						if ($var->{cds_pos_right} >= $self->{cds_length}-3) { $var->{lost_start} = 1; }
						if ($var->{cds_pos_right} < 4) { $var->{lost_stop} = 1; }
					}					
				} # deletion

			} # coding feature

			if ($var->{vartype} ne "D") {
				$var->{stype} = $feature{"type"};
				#$var->{gene_pos} = $var->{position}-$self->{start}+1;
			}
			else {
				#if ($var->{gene_pos} == 0) {
				#	if ($var->{global_start} < $self->{start} && $self->{orientation} eq "+") {
				#		$var->{gene_pos} = 1;
				#	}
				#	elsif ($var->{end} > $self->{end} && $self->{orientation} eq "-") {
				#		$var->{gene_pos} = 1;
				#	}
				#	elsif ($var->{global_start} == $var->{position} && $self->{orientation} eq "+") {
				#		$var->{gene_pos} = $var->{position} - $self->{start} + 1;
				#	}
				#	elsif ($var->{end} <= $feature{"end"} && $self->{orientation} eq "-") {
				#		$var->{gene_pos} = $var->{end} - $self->{start} + 1;
				#	}
				#}
				$var->{stype} .= "-" if ($var->{global_start} != $var->{position} && $var->{stype} ne "");
				$var->{stype} .= $feature{"type"};

				if ($var->{global_start} != $var->{position} || $var->{end} > $feature{"end"}) {
					$var->{stype} .= ",".(min($var->{end}, $feature{"end"}) - $var->{position} + 1);
					$var->{position} = $feature{"end"}+1 if ($var->{end} > $feature{"end"});
				}
			}
print "type: ".$var->{stype}." var->position ".$var->{position}." var->end ".$var->{end}."\n"if ($self->{debug});
#print "  gene_pos: ".$var->{gene_pos}."\n" if ($self->{debug});

			# next variation:
			if ($var->{vartype} ne "D" || $var->{end} <= $feature{"end"}) {
				$varlist_idx++;
				if (!defined $self->{Varlist}[$varlist_idx]) {
					$last_var = 1;
					last;
				}
				$var = $self->{Varlist}[$varlist_idx];
print "@@@@@@@@@@@@@@\nVAR: ".$var->{vartype}." @ ".$var->{position}."\n@@@@@@@@@@@@@@\n" if ($self->{debug});
			}

print "LAST_VAR? $last_var\n" if ($self->{debug});

		} # for each variation in feature

		# fill up coding sequences till the end of coding feature and update $cds_pos
		if (exists $self->{coding_features}{$feature{"type"}}) {
#print $self->{coding_seq_AllVar}."\t".$self->{coding_seq_ref}."\t".$cds_pos."\t".$feature{"end"}."\t".$genomic_pos."\n";
			if ($feature{"end"}-$genomic_pos+1 > 0) {
				my $seq = substr( $self->{coding_seq_ref}, $cds_pos-1, $feature{"end"}-$genomic_pos+1);
				$self->{coding_seq_AllVar} .= $seq;
				$self->{align_ref} .= $seq;
				$self->{align_alt} .= $seq;
			}
			$cds_pos += $feature{"end"}-$genomic_pos+1;
		}

		$last_feat_start = $feat_start;
		%last_feature = %feature;

	} # for each feature

	if (scalar(@{$self->{coding_Vars}}) == 0) {
print "no coding vars\n" if ($self->{debug});
		return(1);
	}

	# check deletions of this gene for premature stop and new stop codons
	foreach my $var (@{$self->{coding_Vars}}) {
print "cds del len: ".($var->{cds_pos_right}-$var->{cds_pos_left}+1)."\n" if ($self->{debug});
		if ($var->{vartype} eq "D" && ($var->{cds_pos_right}-$var->{cds_pos_left}+1) % 3 != 0) {
			$var->{frameshift} = 1;
			$var->{alt_coding} = substr( $self->{coding_seq_ref}, 0, $var->{cds_pos_left}-1 ) if ($var->{cds_pos_left} > 1);
			$var->{alt_coding} .= substr( $self->{coding_seq_ref}, $var->{cds_pos_right}) if ($var->{cds_pos_right} != $self->{cds_length});
			if ($var->{alt_coding} ne "") {
				$self->compare_proteins($var, ($self->{orientation} eq "+") ? $var->{cds_pos_left} : $var->{cds_pos_right});
			}
		}
	}




	if ($self->{orientation} eq "-") {
		$self->{coding_seq_ref} = reverse $self->{coding_seq_ref};
		$self->{coding_seq_ref} =~ tr/ACGTUacgtu/TGCAATGCAA/;

#                $self->{coding_seq_SNP} = reverse $self->{coding_seq_SNP};
#                $self->{coding_seq_SNP} =~ tr/ACGTUacgtu/TGCAATGCAA/;

                $self->{coding_seq_AllVar} = reverse $self->{coding_seq_AllVar};
                $self->{coding_seq_AllVar} =~ tr/ACGTUacgtu/TGCAATGCAA/;

                $self->{align_ref} = reverse $self->{align_ref};
                $self->{align_ref} =~ tr/ACGTUacgtu/TGCAATGCAA/;

                $self->{align_alt} = reverse $self->{align_alt};
                $self->{align_alt} =~ tr/ACGTUacgtu/TGCAATGCAA/;
	}

	# get protein sequences:
#	my $refseq = $self->{align_ref};
#	my $altseq = $self->{align_alt};
#	$refseq =~ tr/-//;
#	$altseq =~ tr/-//;
#	my $ref_aa; my $alt_aa;
#	for (my $p = 0; $p < max(length($refseq, $altseq)); $p += 3) {
#		if ($p < length($refseq)-2 && $ref_aa ne "*") {
#			$ref_aa = $AA{ substr($refseq, $p, 3) };
#			$self->{protein}{"ref"} .= $ref_aa;
#		}
#		if ($p < length($altseq)-2 && $alt_aa ne "*") {
#			$alt_aa = $AA{ substr($altseq, $p, 3) };
#			$self->{protein}{"alt"} .= $alt_aa;
#		}
#	}
	


	return(1);
}



sub compare_proteins
{
	my ( $self, $var, $cdspos ) = @_;
	
	my $altaa = "";

	if ($self->{orientation} eq "-") {
		$var->{alt_coding} = reverse $var->{alt_coding};
		$var->{alt_coding} =~ tr/ACGTUacgtu/TGCAATGCAA/;

		$cdspos = $self->{cds_length}-$cdspos+($var->{vartype} eq "D");
	}

	my $start = int(($cdspos-1) / 3) * 3;
print $var->{vartype}." ".$var->{global_start}." alt_coding: ".$var->{alt_coding}." len: ". length($var->{alt_coding})." (".$self->{cds_length}.") start: $start\n" if ($self->{debug});
	for (my $pos = $start; $pos < length($var->{alt_coding})-2; $pos+=3) {
		$altaa = $AA{ substr($var->{alt_coding}, $pos, 3) };
print "codon ".substr($var->{alt_coding}, $pos, 3)." ".$altaa." pos: $pos cds_pos: $cdspos\n" if ($self->{debug});
		#if ($pos < length($var->{alt_coding})-5 && $altaa eq "*") {
		if ($altaa eq "*") {
			$var->{new_stop} = ($pos-$start)/3+1;
			last;
		}
	}

	if ($altaa ne "" && $altaa ne "*") {
		$var->{lost_stop} = 1;
	}

	return(1);
}


#	$self->{protein}{'ref'} = Bio::Perl::translate_as_string($self->{coding_seq_ref});
##        $self->{protein}{'snp'} = Bio::Perl::translate_as_string($self->{coding_seq_SNP});
#        $self->{protein}{'var'} = Bio::Perl::translate_as_string($self->{coding_seq_AllVar});
#
#
#	for(my $j = 1; $j <= length($self->{protein}{'snp'}); $j++) {
#
#		my $ref_aa = substr($self->{protein}{'ref'}, $j - 1, 1);
#		my $new_aa = substr($self->{protein}{'snp'}, $j - 1, 1);
#
#		for (my $codon_pos = 0; $codon_pos < 3; ++$codon_pos) {
#
#			if( exists $self->{coding_SNP}{(3 * $j) - $codon_pos} ) {
#				my $genome_position = $self->{coding_SNP}{(3 * $j) - $codon_pos};
#				$self->{SNP_lists}{$genome_position}{codon_pos} = 3-$codon_pos;
#				$self->{SNP_lists}{$genome_position}{ref_aa} = $ref_aa;
#				$self->{SNP_lists}{$genome_position}{new_aa} = $new_aa;
#
#				$self->{SNP_lists}{$genome_position}{protein_ref} = $self->{protein}{'ref'};
#				$self->{SNP_lists}{$genome_position}{protein_alt} = $self->{protein}{'snp'};
#
#				$self->{SNP_lists}{$genome_position}{codon} = substr( $self->{eco_coding}, (3*$j)-3, 3 );
#				if (index($self->{SNP_lists}{$genome_position}{codon}, '-') < 0 && index($self->{SNP_lists}{$genome_position}{codon}, 'N') < 0) {
#	#print "genomepos ".$genome_position."\t".$self->{SNP_lists}{snps}{$genome_position}{codon}."\n";
#					if ($self->{orientation} eq "-") {
#						$self->{SNP_lists}{$genome_position}{codon} = reverse($self->{SNP_lists}{$genome_position}{codon});
#						$self->{SNP_lists}{$genome_position}{codon} =~ tr/acgtuACGTU/TGCATGCAA/;
#					}
#					$self->{SNP_lists}{$genome_position}{degeneracy} = substr($DEGENERACY{$self->{SNP_lists}{$genome_position}{codon}}, $self->{SNP_lists}{$genome_position}{codon_pos} - 1, 1);
#				}
#				else {
#					$self->{SNP_lists}{$genome_position}{degeneracy} = "-";
#				}
#
#				if( $ref_aa ne $new_aa ) {
#					$self->{AA_changes}{$j} = 1;
#					$self->{NS_changes}{(3 * $j) - $codon_pos} = $genome_position;
#					$self->{SNP_lists}{$genome_position}{ns_change} = 1;
#					if   ($new_aa eq '*') {$self->{SNP_lists}{$genome_position}{new_stop} = 1;}
#					elsif($ref_aa eq '*') {$self->{SNP_lists}{$genome_position}{lost_stop} = 1;}
#				}
#			}
#
#		}
#	}
#
#	return(1);
#}


# OLD:
				
#		if( exists $self->{coding_SNP}{(3 * $j) - 1} ) { 
#			my $genome_position = $self->{coding_SNP}{(3 * $j) - 1};
#			$self->{SNP_lists}{$genome_position}{codon_pos} = 2;
#			$self->{SNP_lists}{$genome_position}{ref_aa} = $ref_aa;
#			$self->{SNP_lists}{$genome_position}{new_aa} = $new_aa;
#
#			$self->{SNP_lists}{$genome_position}{protein_ref} = $self->{protein}{'ref'};
#			$self->{SNP_lists}{$genome_position}{protein_alt} = $self->{protein}{'alt'};
#
#			$self->{SNP_lists}{$genome_position}{codon} = substr( $self->{eco_coding}, (3*$j)-3, 3 );
#			if (index($self->{SNP_lists}{$genome_position}{codon},'-') < 0 && index($self->{SNP_lists}{$genome_position}{codon}, 'N') < 0) {
##print $self->{SNP_lists}{snps}{$genome_position}{codon}."\n";
#				if ($self->{orientation} eq "-") {
#					$self->{SNP_lists}{$genome_position}{codon} = reverse($self->{SNP_lists}{$genome_position}{codon});
#					$self->{SNP_lists}{$genome_position}{codon} =~ tr/acgtuACGTU/TGCATGCAA/;
#				}
#				$self->{SNP_lists}{$genome_position}{degeneracy} = substr($DEGENERACY{$self->{SNP_lists}{$genome_position}{codon}}, $self->{SNP_lists}{$genome_position}{codon_pos} - 1, 1);
#			}
#			else {
#				$self->{SNP_lists}{$genome_position}{degeneracy} = "-";
#			}
#
#
#			if( $ref_aa ne $new_aa ) {
#				$self->{AA_changes}{$j} = 1;
#				$self->{NS_changes}{(3 * $j) - 1} = $genome_position;
#				$self->{SNP_lists}{$genome_position}{ns_change} = 1;
#				if   ($new_aa eq '*') {$self->{SNP_lists}{$genome_position}{new_stop} = 1;}
#				elsif($ref_aa eq '*') {$self->{SNP_lists}{$genome_position}{lost_stop} = 1;}
#			}
#		}
#				
#		if( exists $self->{coding_SNP}{3 * $j} ) {
#			my $genome_position = $self->{coding_SNP}{3 * $j};
#			$self->{SNP_lists}{$genome_position}{codon_pos} = 3;
#			$self->{SNP_lists}{$genome_position}{ref_aa} = $ref_aa;
#			$self->{SNP_lists}{$genome_position}{new_aa} = $new_aa;
#
#			$self->{SNP_lists}{$genome_position}{protein_ref} = $self->{protein}{'ref'};
#			$self->{SNP_lists}{$genome_position}{protein_alt} = $self->{protein}{'alt'};
#
#			$self->{SNP_lists}{$genome_position}{codon} = substr( $self->{eco_coding}, (3*$j)-3, 3 );
#			if (index($self->{SNP_lists}{$genome_position}{codon}, '-') < 0 && index($self->{SNP_lists}{$genome_position}{codon}, 'N') < 0) {
##print $self->{SNP_lists}{snps}{$genome_position}{codon}."\n";
#				if ($self->{orientation} eq "-") {
#					$self->{SNP_lists}{$genome_position}{codon} = reverse($self->{SNP_lists}{$genome_position}{codon});
#					$self->{SNP_lists}{$genome_position}{codon} =~ tr/acgtuACGTU/TGCATGCAA/;
#				}
##if (!defined $DEGENERACY{$self->{SNP_lists}{snps}{$genome_position}{codon}} || !defined $self->{SNP_lists}{snps}{$genome_position}{codon}) {
##	print $self->{eco_coding}."\n".$self->{SNP_lists}{snps}{$genome_position}{codon}."\n".((3*$j)-3)."\n".$genome_position."\n";
##	exit;
##}
#				$self->{SNP_lists}{$genome_position}{degeneracy} = substr($DEGENERACY{$self->{SNP_lists}{$genome_position}{codon}}, $self->{SNP_lists}{$genome_position}{codon_pos} - 1, 1);
#			}
#			else {
#				$self->{SNP_lists}{$genome_position}{degeneracy} = "-";
#			}
#
#
#			if( $ref_aa ne $new_aa ) {
#				$self->{AA_changes}{$j} = 1;
#				$self->{NS_changes}{3 * $j} = $genome_position;
#				$self->{SNP_lists}{$genome_position}{ns_change} = 1;
#				if   ($new_aa eq '*') {$self->{SNP_lists}{$genome_position}{new_stop} = 1;}
#				elsif($ref_aa eq '*') {$self->{SNP_lists}{$genome_position}{lost_stop} = 1;}
#			}
#		}
##	}
#	
##	return(1);
##}


