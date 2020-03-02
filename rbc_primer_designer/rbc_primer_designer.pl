#!/usr/bin/perl

###############################################################################
#
# rbc_primer_designer -- Designs primers for given variants in rbci format
#
# Version: 0.1
#
# Author: 
#   Joffrey Fitz <joffrey.fitz@tuebingen.mpg.de>
#
# Copyright (c) 2012 Max Planck Institute for Developmental Biology, 
#   http://www.eb.tuebingen.mpg.de
#
# This program is free software; you can redistribute it and/or modify it                        
# under the terms of the GNU General Public License as published by the Free                     
# Software Foundation; either version 2, or (at your option) any later 
# version.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
# for more details.
#
###############################################################################

my $VERSION = "0.1";

#
# Global values
###############################################################################
my $VERBOSE    => 0;
use constant DEBUG      => 0;
use constant DEBUG_P3   => 0;
use constant DEBUG_MFEP => 0;


# Primer product sizes
my $PRODUCT_SIZE_RANGE_MIN = 280;
my $PRODUCT_SIZE_RANGE_MAX = 400; 
my $PRODUCT_OPT_SIZE = 300;

# Regions for primers
my $P_DIST_MAX     = 50;   # max distance from SNP
my $P_DIST_MIN     = 3;    # min distance from SNP


# Variants stuff
# skip secondary variants if indel length larger than this
my $MAX_SECONDARY_INDEL_LEN = 5; 

#
# 3rd party tools
###############################################################################

#
# Executables
#
my $fastacmd_bin = `which fastacmd`; chomp $fastacmd_bin;
my $p3_bin = `which primer3_core`; chomp $p3_bin;
my $mfeprimer_bin = `which MFEprimer.py`; chomp $mfeprimer_bin;
my $p3_conf = "";

my $mfeprimer_args_intern = "-T F -W 4 -a 2 --ppc_cutoff=0.255";
my $mfeprimer_args = "";

my $PYTHON2 = `which python2`; chomp $PYTHON2; #added by Norman Jan 2020
 
###############################################################################
# END OF CONFIGURATION
#


use strict;
use warnings;

use File::Temp;
use Getopt::Long qw(:config no_ignore_case);

# Accepted Primer3 tags are:
# PRIMER_LEFT_0
# PRIMER_RIGHT_0
# PRIMER_LEFT_0_SEQUENCE
# PRIMER_RIGHT_0_SEQUENCE
# PRIMER_PAIR_0_PRODUCT_SIZE
# PRIMER_LEFT_0_TM
# PRIMER_RIGHT_0_TM
# PRIMER_LEFT_0_GC_PERCENT	
# PRIMER_RIGHT_0_GC_PERCENT
# PRIMER_LEFT_0_SELF_ANY_TH
# PRIMER_RIGHT_0_SELF_ANY_TH
# PRIMER_LEFT_0_SELF_END_TH
# PRIMER_RIGHT_0_SELF_END_TH
# PRIMER_LEFT_0_HAIRPIN_TH
# PRIMER_RIGHT_0_HAIRPIN_TH
# PRIMER_LEFT_0_END_STABILITY
# PRIMER_RIGHT_0_END_STABILITY
# PRIMER_LEFT_0_TEMPLATE_MISPRIMING_TH
# PRIMER_RIGHT_0_TEMPLATE_MISPRIMING_TH
# PRIMER_PAIR_0_TEMPLATE_MISPRIMING_TH
# PRIMER_LEFT_EXPLAIN
# PRIMER_RIGHT_EXPLAIN
# PRIMER_LEFT_EXPLAIN
# PRIMER_PAIR_EXPLAIN
# PRIMER_ERROR

# Primer 3 result tags to be reported
my @p3_res_tags = qw/
PRIMER_LEFT_0
PRIMER_RIGHT_0
PRIMER_PAIR_0_PRODUCT_SIZE
PRIMER_LEFT_0_SEQUENCE
PRIMER_RIGHT_0_SEQUENCE
/;

my $prog_name = "rbc_primer_designer";


#
# Command line interface
#

my $run_mode = "design_primers";

my ($print_version, $print_help);
my ($loc_seq_db, $loc_rbci);
my ($chr, $start, $end);
my ($usr_filter, $restriction_seq) = ("PASS", "");
my $ann_major = 0; # use IUPAC base or major allel as annotation for primary SNP

my $result = GetOptions(
	"loc_seq_db|b=s"               => \$loc_seq_db,
	"loc_rbci|r=s"                 => \$loc_rbci,
	
	"chr|c=s"               	   => \$chr,
	"start|s=s"                    => \$start,
	"end|e=i"                      => \$end,
	"usr_filter|f=s"               => \$usr_filter,
	"restriction_seq|x=s"          => \$restriction_seq,
	"PRODUCT_SIZE_RANGE_MIN|Q=s"   => \$PRODUCT_SIZE_RANGE_MIN,
	"PRODUCT_SIZE_RANGE_MAX|P=s"   => \$PRODUCT_SIZE_RANGE_MAX,
	"PRODUCT_OPT_SIZE|O=s"         => \$PRODUCT_OPT_SIZE,
	"P_DIST_MIN|R=s"               => \$P_DIST_MIN,
	"P_DIST_MAX|S=s"               => \$P_DIST_MAX,
	"MAX_SECONDARY_INDEL_LEN|N=s"  => \$MAX_SECONDARY_INDEL_LEN,
	
	"p3_bin|p=s"                   => \$p3_bin,
	"p3_conf|d=s"                  => \$p3_conf,
	"mfeprimer_bin|m=s"            => \$mfeprimer_bin,
	"mfeprimer_args=s"             => \$mfeprimer_args,
	
	"major"                        => \$ann_major,
	
	"VERBOSE|v+"                   => \$VERBOSE,
	
	"version|V"                    => \$print_version,
	"help|h"                       => \$print_help,
);			

my $SNP_FLANKING = $PRODUCT_SIZE_RANGE_MAX;

if(defined $print_version) {
	print_version();
}

if(defined $print_help) {
	print_help();
}

# Check mandatory arguments
if($p3_conf eq "" || $loc_seq_db eq "") {
	print_help();
}

# Check executables
die "[ERROR] Executable for primer3 not found: '$p3_bin'\n" if (! -e $p3_bin);
die "[ERROR] Executable MFGprimer not found: '$mfeprimer_bin'\n" if (! -e $mfeprimer_bin);
die "[ERROR] Executable for fastacmd not found: '$fastacmd_bin'\n" if (! -e $fastacmd_bin);
die "[ERROR] '$p3_bin' not executable\n" if (! -x $p3_bin);
die "[ERROR] '$mfeprimer_bin' not executable\n" if (! -x $mfeprimer_bin);
die "[ERROR] '$fastacmd_bin' not executable\n" if (! -x $fastacmd_bin);

# Check data files
die "[ERROR] rbci input file does not exist: '$loc_rbci'\n" if(! -e $loc_rbci);
die "[ERROR] BLAST DB does not exist: '$loc_seq_db'\n" if(! -e $loc_seq_db);
die "[ERROR] Primer3 configuration does not exist: '$p3_conf'\n" if(! -e $p3_conf);

# Set mfeprimer_args
if($mfeprimer_args eq "") {
	$mfeprimer_args = $mfeprimer_args_intern . " -e " . 10*$PRODUCT_SIZE_RANGE_MAX;
}

if($VERBOSE) {
	print_params();
}

# init primary and secondary SNPs from rbci file
# use all PASS SNPs as primaries, all other as secondary
my ($prim, $prim_major_minors, $second, $second_types, $second_filters, $rbci) = init_rbci($loc_rbci, $chr);
my %prim_snplist = %$prim;
my %major_minors_list = %$prim_major_minors;
my %second_snplist = %$second;
my %rbci_vals = %$rbci;

#
# Run mode design_primers
#
# Input: rbci file
# Output: rbcp file
#
my $MAJOR_ANN_FLAG = "MAJOR";
if(!$ann_major) {
	$MAJOR_ANN_FLAG = "";
}
if($run_mode eq "design_primers") {
	
	my $snp_nr=1;
	foreach my $snp_pos (sort {$a<=>$b} keys %prim_snplist) {
		my $direction="fwd";
		
		my ($header, $ann_seq) = get_template_seq($snp_pos, $direction, "PRIMER_SEQ", "HEADER", "$MAJOR_ANN_FLAG");

		warn "[INFO] Trying $direction primers for SNP no $snp_nr, $chr:$snp_pos\n" if($VERBOSE);
	
		my $res = p3_res_2tab( run_primer3("$header", "$ann_seq"), "$header" );
		my @out = split("\t", $res );

		warn "[INFO]     Got primer3 results: @out\n" if($VERBOSE>1);
		
		my $found_primers = 0;
		if( not (($out[1] eq "N/A") or ($out[2] eq "N/A")) ) {
			# Run mfe_primer and check primer maps on genome
			warn "$out[1] $out[2]" if(DEBUG_MFEP);

			my $hits = run_mfeprimer($out[4], $out[5], $loc_seq_db);
			push @out, $hits;
			warn "[INFO]     MFP found $hits hits\n" if($VERBOSE);
			if($hits == 1) {
				warn "[INFO]     Found Primers ($direction, $hits off-target) for SNP no $snp_nr, $chr:$snp_pos.\n" if($VERBOSE);
				
				print rcbp_output($chr, $snp_pos, $direction, $out[4], $out[5],
					get_product_seq($ann_seq, $snp_pos, $out[1], $out[2], $direction, $restriction_seq)) ."\n";
						
				$found_primers = 1;
			}
		}
		if($found_primers != 1) { # no primers found
			my $direction="rev";
			my ($header, $ann_seq) = get_template_seq($snp_pos, $direction, "PRIMER_SEQ", "HEADER", "$MAJOR_ANN_FLAG");

			warn "[INFO] Trying $direction primers for SNP no $snp_nr, $chr:$snp_pos\n" if($VERBOSE);
			
			my $res = p3_res_2tab( run_primer3("$header", "$ann_seq"), "$header" );
			my @out = split("\t", $res );
		
			if( not (($out[1] eq "N/A") or ($out[2] eq "N/A")) ) {
				# Run mfe_primer and check primer maps on genome
				warn "$out[1] $out[2]" if(DEBUG_MFEP);

				my $hits = run_mfeprimer($out[4], $out[5], $loc_seq_db);
				push @out, $hits;
				warn "[INFO]     MFP found $hits hits\n" if($VERBOSE);
				if($hits == 1) {
					print rcbp_output($chr, $snp_pos, $direction, $out[4], $out[5],
						get_product_seq($ann_seq, $snp_pos, $out[1], $out[2], $direction, $restriction_seq)) ."\n";
							
					warn "[INFO] Found Primers ($direction, $hits off-target) for $chr:$snp_pos.\n" if($VERBOSE);
				} else {
					warn "[INFO] NO primers found for SNP no $snp_nr, $chr:$snp_pos.\n";
				}
			}
		}
		$snp_nr++;
	}

	exit;
}

exit;



#
# Starting subs
###############################################################################



#
# Translate given ACGT SNP (e.g. "AT") to IUPAC
#
# @return: $iupac_base, e.g. W
sub to_iupac {
	my ($snp) = shift;
	#IUPAC Codes
	#AT => W
	#CG => S
	#AC => M
	#GT => K
	#CT => Y
	#AG => R
	#CGT => B
	#AGT => D
	#ACT => H
	#ACG => V
	#ACGT => N

	$snp =~ s/AT|TA/W/;
	$snp =~ s/CG|GC/S/;
	$snp =~ s/AC|CA/M/;
	$snp =~ s/GT|TG/K/;
	$snp =~ s/CT|TC/Y/;
	$snp =~ s/AG|GA/R/;
	return $snp;
}



#
# Read SNP list in rbci format
# from disk and store in ds.
#
# @return: %variantlist, $pos => <ref>:<variant>, 
#	e.g. A:AACC (insertion)
#		 GTTA:G (deletion)
#		 A:G (SNP)
#
# and 
# %major_minor_list, $pos => <major>:<minor>, 
sub init_rbci {
	#rbci format: 
	#CHROM  POS     ID      REF     ALT     MAJOR   MINOR   QUAL    FILTER  TYPE    SEGR    HAPSC   FORMAT  rice_pools_res  rice_pools_sus
	my ($loc_snp_list, $chr) = @_;
	
	my %prim_snps;
	my %prim_major_minors;
	my %second_snps;
	my %second_types;
	my %second_filters;
	my %rbci;	# store QUAL, FILTER, TYPE, SEGR
	
	my $pass_cnt = 0;
	
	my @flags_ok = split ";", $usr_filter;
	warn "[INFO] Allowed filter flags: @flags_ok\n" if($VERBOSE);
	
	open(FH, "<$loc_snp_list") or die "Can not open SNP list $loc_snp_list: $!";
	while (<FH>) {
		next if(/^#/);
		chomp;

		my @l=split(/\t/, $_);
		next if($l[0] ne $chr);
		next if($l[1] <= $start);
		next if($l[1] >= $end);
				
		my ($key, $filter, $type) = @l[1,8,9];
		my @flags = split ";", $filter;
		my $filter_ok = 0;
				
		if(scalar(@flags_ok)>0) {
			my @union = (); my @isect = ();
			my %union = (); my %isect = ();
			
			foreach my $e (@flags, @flags_ok) { $union{$e}++ && $isect{$e}++ }
			@isect = keys %isect;
			
			warn "flags=(@flags), flags_ok=(@flags_ok), isect=(@isect)\n" if(DEBUG);
			
			if( (scalar(@flags)>scalar(@isect))  ) {
				warn "    ---> FAILED\n" if(DEBUG); 
			} else {
				$filter_ok = 1;
			}
		}
		if($filter_ok) { # potential primary SNP
			#1. get list of legal flags, e.g "HARD_TO_VALIDATE;MaxAC;NoSBtag"
			     			
			if($type eq "SNP") {
				$prim_snps{$key} = $l[3] . $l[4];
				$prim_major_minors{$key} = $l[5] . $l[6];
				$rbci{$key} = "$l[7]:$l[8]:$l[9]:$l[10]";
			}
			$pass_cnt++;
		} else { # secondary SNP
			#warn "[INFO] found type=$type, treating as secondary variant\n";
			if($type eq "indel") {
				#warn "[INFO] found indel as secondary variant\n";
				if(length($l[3]) >= $MAX_SECONDARY_INDEL_LEN 
					or length($l[4]) >= $MAX_SECONDARY_INDEL_LEN) {
						warn "[!!] Found secondary variant indel length >= $MAX_SECONDARY_INDEL_LEN at $chr:$key ($l[3]/$l[4]). Skipping.\n" if($VERBOSE>=2);
						next;
					} 
			}
			$second_snps{$key} = $l[3] . ":" . $l[4];
			$second_types{$key} = $type;
			$second_filters{$key} = $filter;
		}
	}
	
	if($VERBOSE) {
		warn "\n[INFO] Variant file loaded:\n";	
		warn "[INFO]     Passed variants: $pass_cnt\n";
		warn "[INFO]     Primary SNPs: Read ". scalar(keys(%prim_snps)) . " SNPs\n" ;
		warn "[INFO]     Primary SNPs: Read ". scalar(keys(%prim_major_minors)) . " major/minor\n";
		warn "[INFO]     Secondary variants: Read ". scalar(keys(%second_snps)) . " variants\n";
		warn "[INFO]     Secondary variants:\n";
		my %counts;
		my %unique;
		foreach my $key (sort keys %second_types) {
		    my $value = $second_types{$key}; 
		    if (not exists $counts{$value}) {
		        $unique{$key} = $value;
		    }
		    $counts{$value}++;
		}
		foreach my $key (keys %counts) {
			warn "[INFO]     $key => " . $counts{$key} . "\n";
		}
		warn "\n";
	}

	return (\%prim_snps, \%prim_major_minors, \%second_snps, \%second_types, \%second_filters, \%rbci);
}



#
# extract subsequences
#
sub get_seq {
	my ($loc_seq_db, $chr, $start, $end) = @_;
	#get_seq_cut($loc_seq_db, $chr, $start, $end);
	get_seq_fastacmd($loc_seq_db, $chr, $start, $end);
}



#
# Extract subsequences using fastacmd
#
sub get_seq_fastacmd {
	my ($loc_seq_db, $chr, $start, $end) = @_;
	
	warn "in get_seq_fastacmd()" if(DEBUG);
	my $cmd = "$fastacmd_bin -d $loc_seq_db -s $chr -L $start,$end";
	warn "cmd: $cmd" if(DEBUG);
	open P,"$cmd |" or die "error running command $cmd: $!";
	my $seq;
	while(<P>) {
		next if(/^>/);
		chomp;
		$seq .= $_;
	}
	close P;
	my $exit_value=$? >> 8;
	die "something went wrong: \n$exit_value\n$seq" if($exit_value != 0);
	return $seq;
}



#
# Annotate seq with primary and secondary SNPs 
#
sub get_template_seq {
	my ($snp_pos, $direction, $mapping_seq, $return_header, $major_allel, $secondary_no_n) = @_;
	$mapping_seq = 0 if (not defined $mapping_seq);
	$return_header = 0 if (not defined $return_header);
	$major_allel = 0 if (not defined $major_allel);
	$secondary_no_n = 0 if (not defined $secondary_no_n);
	
	warn "[INFO] Extracting sequence for SNP at $snp_pos, $snp_pos-$SNP_FLANKING ... $snp_pos+$SNP_FLANKING\n" if($VERBOSE>1);
	
	warn "checking primary snp at $snp_pos\n" if($VERBOSE>3);
	my $seq_start_pos = $snp_pos-$SNP_FLANKING;
	my $seq_end_pos = $snp_pos+$SNP_FLANKING;
	if($seq_start_pos <1) {
		warn "[!!] Skipping primary SNP $snp_pos,  \$seq_start_pos=$seq_start_pos <1";
		next;
	}

	my $snp_pos_rel = $SNP_FLANKING;
	warn "\n= found primary snp " . $prim_snplist{$snp_pos} . " at $snp_pos ($snp_pos_rel)\n" if(DEBUG);

	my $seq = &get_seq($loc_seq_db, $chr, $seq_start_pos, $seq_end_pos );
	my $ref_seq = $seq;
	my @seq_a = split('', $seq);
	warn "  - Got seq from $seq_start_pos .. $seq_end_pos\n" if(DEBUG);
	
	my $snp_str = "";
	if($major_allel eq "MAJOR") {
		$snp_str = $major_minors_list{$snp_pos};
	} else {	
		$snp_str = $prim_snplist{$snp_pos};
	}
	
	# We do not allow indels / tri-allelic SNPs
	if(length($snp_str)>2) {
		warn "[!!] Indel or tri-allelic SNP found at $chr,$snp_pos:'$snp_str'. Skipping.";
		return;
	}

	if($major_allel eq "MAJOR") {
		my ($major_allel, $minor_allel) = split '', $snp_str;
		warn "[INFO]     Annotating major allele $snp_pos_rel => $major_allel\n" if($VERBOSE>1);
		$seq_a[$snp_pos_rel] = $major_allel; 
	} else {
		my $iupac_snp = to_iupac( $snp_str );
		$seq_a[$snp_pos_rel] = $iupac_snp;
		warn "[INFO]     Annotating IUPAC allele $snp_pos_rel => $iupac_snp\n" if($VERBOSE>1);
	}
	
	if($mapping_seq ne "MAPPING_SEQ") {
		# Annotate secondary SNPs and indels
		foreach my $other_snp (sort {$a<=>$b} keys %second_snplist) {
			last if($other_snp > $snp_pos + $SNP_FLANKING);
			next if($other_snp < $snp_pos - $SNP_FLANKING);
			next if ($other_snp == $snp_pos);

			my $other_snp_rel = $other_snp-$seq_start_pos;
		
			my ($other_variant_ref, $other_variant_call) = split(":", $second_snplist{$other_snp});
			warn "[INFO]     Variant found at $other_snp: $other_variant_ref, $other_variant_call\n" if($VERBOSE>1);
	
			# Do we have an tri-allelic call (like A/T,C)?
			# Only care about first base, treat as SNP (e.g. A/T)
			if($other_variant_call =~ m/(.?),/) {
				warn "[!!]     Got an tri-allelic variant: $other_variant_call\n";
				$other_variant_call = $1;
				warn "[!!]         => Treating only $other_variant_call as SNP\n";
			}

			# Deletion found? Replace insertion seq with Ns
			if(length($other_variant_ref)>length($other_variant_call)) {
				my $ins_len = length($other_variant_ref);
				warn "[INFO]     Get_template_seq: Deletion found (len=$ins_len)\n" if($VERBOSE>1);
				for(my $i=$other_snp_rel; $i<$other_snp_rel+$ins_len; $i++) {
					$seq_a[$i] = "N";
					warn "$seq_a[$i]" if(DEBUG);
				}
			}
			if(length($other_variant_ref)<length($other_variant_call)) {
				warn "[INFO]     Get_template_seq: Insertion found\n\n\n" if($VERBOSE>1);
				# Annotate insertion with N. As insertions are annotated like T/TACC, N out the n+1st base
				$seq_a[$other_snp_rel+1] = "N";
			}
			else {
				# Annotate SNP one (1) N
				$seq_a[$other_snp_rel] = "N";	
			}
		}
	}

	my $ann_seq = join ("", @seq_a);
	
	my $header;
	if($direction eq "rev") {
		$ann_seq =~ tr/ACTGactgMRVHmrvhKYBDkybd/TGACtgacKYBDkybdMRVHmrvh/;
		$ann_seq = reverse($ann_seq);
		my $rev_snp_str = $snp_str;
		$rev_snp_str =~ tr/ACTGactgMRVHmrvhKYBDkybd/TGACtgacKYBDkybdMRVHmrvh/;
		my $rev_snp = $seq_a[$snp_pos_rel];
		$rev_snp =~ tr/ACTGactgMRVHmrvhKYBDkybd/TGACtgacKYBDkybdMRVHmrvh/;

		#$header = "SNP_${chr}:${snp_pos}_${direction}_rev | $rev_snp_str annotated as " . $rev_snp;
		$header = "${chr}_${snp_pos}_${direction}_${seq_start_pos}_${seq_end_pos}";
	} else {
		my $seq_len = length($ann_seq);
		#$header = "SNP_${chr}:${snp_pos}_${direction}_fwd ($seq_start_pos .. $seq_end_pos, len=$seq_len) | $snp_str annotated as " . $seq_a[$snp_pos_rel];
		$header = "${chr}_${snp_pos}_${direction}_${seq_start_pos}_${seq_end_pos}";
	}
	
	warn "[INFO]     ann_seq (${direction}): $ann_seq\n" if($VERBOSE>2);
	
	if($return_header eq "NO_HEADER") {
		return $ann_seq;
	} else {
		return ($header, $ann_seq);
	}
}



#
# Returns primer product, secondary SNPs annotated as Ns, with optional added restriction_seq
#
# params: $ann_seq, $snp_pos, $p3_start, $p3_end, $direction, $restriction_seq
# return: $rel_snp_pos, $product_seq, $N_cnt
#
sub get_product_seq {
	warn "[INFO] Getting product sequence\n" if($VERBOSE>1);
	
	my ($ann_seq, $snp_pos, $p3_start, $p3_end, $direction, $restriction_seq) = @_;
		
	#my $ann_seq = get_template_seq($snp_pos, $direction, "NO_HEADER", "MAJOR", "SECONDARY_NO_N");
	#my $ann_seq = get_template_seq($snp_pos, $direction, "MAPPING_SEQ", "NO_HEADER", "MAJOR");

	my ($primer_left_start, $primer_left_len) = split ",", $p3_start;
	my ($primer_right_start, $primer_right_len) = split ",", $p3_end;
	my $product_len = $primer_right_start-$primer_left_start+1;
	my $product_seq = substr $ann_seq, $primer_left_start, $product_len;
	
	if($direction eq "fwd") {
		$product_seq = $product_seq . $restriction_seq;
	} elsif($direction eq "rev") {
		$product_seq = $product_seq . $restriction_seq;
		$product_seq =~ tr/ACTGactgMRVHmrvhKYBDkybd/TGACtgacKYBDkybdMRVHmrvh/;
		$product_seq = reverse($product_seq);
		
	} else {
		die "Unknown direction '$direction'";
	}
	
	my $rel_snp_pos = "";
	if($direction eq "rev") {
		 $rel_snp_pos = $primer_right_start-$SNP_FLANKING+1+length($restriction_seq);
		warn "[INFO]     SNP relative position: $rel_snp_pos = $primer_right_start-$SNP_FLANKING+1+" . length($restriction_seq) . "\n" if($VERBOSE>1);
	} elsif($direction eq "fwd") {
		 $rel_snp_pos = $SNP_FLANKING-$primer_left_start+1;
		warn "[INFO]     SNP relative position: $rel_snp_pos = $SNP_FLANKING-$primer_left_start\n" if($VERBOSE>1);
	}
	
	my $N_cnt = grep(/N/, split("",$product_seq) );	 
	
	if($VERBOSE>1) {
		warn "[INFO]     Product seq length: " . length($product_seq) . "\n";
		warn "[INFO]     N count: $N_cnt\n";
	}
	warn "[INFO]     Product seq: $product_seq\n" if($VERBOSE>2);
	
	return($rel_snp_pos, $product_seq, $N_cnt);
}



#
# Run primer3
#
sub run_primer3 {
	my($seq_id, $seq_str) = @_;

	# my $PRODUCT_SIZE_RANGE_MIN = 100;
	# my $PRODUCT_SIZE_RANGE_MAX = 500; 
	# my $PRODUCT_OPT_SIZE = 300;
	my $start1 = $SNP_FLANKING - $P_DIST_MAX;
	my $len1   = $P_DIST_MAX-3;
	my $start2 = $start1 + $PRODUCT_SIZE_RANGE_MIN;
	my $len2   = $PRODUCT_SIZE_RANGE_MAX - $PRODUCT_SIZE_RANGE_MIN;
	
	my $out =<<OUT;
SEQUENCE_ID=$seq_id
SEQUENCE_TEMPLATE=$seq_str
SEQUENCE_TARGET=$SNP_FLANKING,1
SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=$start1,$len1,$start2,$len2
PRIMER_PRODUCT_SIZE_RANGE=$PRODUCT_SIZE_RANGE_MIN-$PRODUCT_SIZE_RANGE_MAX
PRIMER_PRODUCT_OPT_SIZE=$PRODUCT_OPT_SIZE
=
OUT

	my $cmd = qq(echo -n "$out" | $p3_bin -p3_settings_file=$p3_conf);
	warn "cmd: $cmd" if(DEBUG);
	
	open P,"$cmd |" or die "error running command $cmd: $!";

	my $output;
	while(<P>) {
		if(/PRIMER_ERROR=(.*)/) {
			#die "[!!] Primer3 Error: $1\nUsed input:\n$out\nDied";
		}
		$output .= $_;
	}
	my $exit_value=$? >> 8;
			
	die "something went wrong: \n$exit_value\n$output" if($exit_value != 0);
	
	warn "output:\n$output\n" if(DEBUG);
	
	if(DEBUG_P3) {
		$cmd .= " -format_output";
		warn "cmd: $cmd\n";
		open P,"$cmd |" or die "error running command $cmd: $!";
		while(<P>) {
			print STDERR $_;
		}
	}
	
	return $output;
}



#
# Parse and report primer3 results as tab delimited line.
#
# @return: string, tab delimited, new line terminated
#
sub p3_res_2tab {
	my ($p3_res, $header) = @_;
	
	my @line;
	push @line, $header;
	foreach my $tag (@p3_res_tags) {
		warn "checking tag '$tag'\n" if(DEBUG);
		if($p3_res =~/${tag}=(.*)/) {
			push @line, $1;
		}
		else {
			push @line, "N/A";
		}
	}
	my $res = join("\t", @line);
	return  $res;
}

sub p3_res_header_2tab {
	my @header_fileds=@p3_res_tags;
	unshift @header_fileds, "SNP";
	my $line = join("\t", @header_fileds);
	
	return "#" . $line; 
}



#
# Assemble rcbp output format
#
# 
sub rcbp_output {
	my ($chr, $snp_pos, $direction, $p1_seq, $p2_seq, $rel_snp_pos, $product_seq, $N_cnt) = @_;
	
	my $product_len = length($product_seq);
	my ($major_allel, $minor_allel) = split('',$major_minors_list{$snp_pos});
	my ($ref, $alt) = split('',$prim_snplist{$snp_pos});
	my ($qual, $type, $filter, $segr) = split(':', $rbci_vals{$snp_pos});
	
	return( join("\t", (
		$chr, $snp_pos, ".", $ref, $alt, $major_allel, $minor_allel, 
		$qual, $type, $filter, $segr,
		$direction, $p1_seq, $p2_seq, 
		$product_len, $N_cnt, $product_seq,
	)) );
}



#
# run_mfeprimer(1st_primer-seq, 2nd_primer_seq, loc_seq_db)
#
# Generates temporary multifasta sequence of given primers,
# invokes MFEPrimer and returns number of hits found as results
#
# @return: number_of_hits
#
sub run_mfeprimer {
	my ($seq_a, $seq_b, $loc_seq_db) = @_;
	my $fasta=">seq_a\n$seq_a\n>seq_b\n$seq_b\n";
	
	my $tmp = File::Temp->new( TEMPLATE => 'dcof-m_tempXXXXX',
		DIR => '/tmp',
		SUFFIX => '.fa',
		UNLINK => 1);
	my $fname = $tmp->filename;
	print $tmp $fasta;
	warn "Created tmp fasta file for MFEPrimer: $fname" if(DEBUG_MFEP);
	
	my $cmd = qq($PYTHON2 $mfeprimer_bin -i $fname -d $loc_seq_db $mfeprimer_args);
	warn "cmd=$cmd" if(DEBUG_MFEP);
	
	open P,"$cmd 2>/dev/null |" or die "error running command $cmd: $!";
	
	my $output;
	while(<P>) {
		if(/FATAL ERROR: (.*)/) {
			#die "[!!] Primer3 Error: $1\nUsed input:\n$out\nDied";
		}
		#chomp;
		
		# grep in MFEprimer output for
		# Distribution of 33 MFEprimer hits on the query primers
		
		chomp;
		if(m/Distribution of/) {
			$output = (split(" ", $_))[2];
		}
		#$output .= $_;
	}
	my $exit_value=$? >> 8;
	die "something went wrong: \n$exit_value\n$output" if($exit_value != 0);
	
	
	
	warn "output: '$output'\n" if($VERBOSE>2);
	
	return $output;
}



sub print_version {
	warn "$prog_name $VERSION\n\n";
	warn "Written 2012 by Joffrey Fitz <joffrey.fitz\@tuebingen.mpg.de>\n";
	warn "Copyright 2012 Max Planck Institute for Developmental Biology, Tübingen, Germany\n";
	warn "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n";
	warn "This is free software: you are free to change and redistribute it.\n";
	warn "There is NO WARRANTY, to the extent permitted by law.\n\n";
	exit 0;
}



sub print_help {
	print STDERR<<EOF;

$prog_name $VERSION

Usage: $0 [OPTIONS] -r rbcifile -b blastdb -c chr -s start -e end -p primer3_bin -d primer3_conf -m mfprimer_bin

Designs primers for SNPs defined in variants.rbci for specified region (chromosome, start, end). 
Only SNPs with user given flags (defaults to PASS) are considered.
Results are printed to STDOUT in rcbp format.

Mandatory input files:
  -i, --loc_rbci=FILE                (mandatory) use FILE as rbci variant input file
  -b, --loc_seq_db=FASTA             (mandatory) FASTA file within preprocessed BLAST DB.
  
SNP control:
  -c, --chr=ID                       (mandatory) use ID as chromosome identifier
  -s, --start=NUM                    (mandatory) look at region starting from NUM
  -e, --end=NUM                      (mandatory) in conjunction with -s; region ends at NUM      
  -f, --usr_filter=FLAGS             (=PASS)	 semicolon separated list of rbci variant filter flags
      --major                        (optional)  use major allele for SNP annotation instead of IUPAC base

Primer design control:
  -Q, --PRODUCT_SIZE_RANGE_MIN=NUM   (=$PRODUCT_SIZE_RANGE_MIN)	 set minimum primer product range to NUM
  -P, --PRODUCT_SIZE_RANGE_MAX=NUM   (=$PRODUCT_SIZE_RANGE_MAX)	 set maximum primer product range to NUM
  -O, --PRODUCT_OPT_SIZE=NUM         (=$PRODUCT_OPT_SIZE)	 set optimal primer product size to NUM
  -R, --P_DIST_MIN=NUM               (=$P_DIST_MIN)	 set minimal distance of primer to primery SNP to NUM
  -S, --P_DIST_MAX=NUM               (=$P_DIST_MAX)	 set maximal distance of primer to primery SNP to NUM
  -N, --MAX_SECONDARY_INDEL_LEN=NUM  (=$MAX_SECONDARY_INDEL_LEN)	 discard assay, if secondary indel exists longer than NUM

Helper tools:
  -d, --p3_conf=FILE                 (mandatory) use FILE as primer3 configuration
  -p, --p3_bin=BIN                   (optional)  set absolute path of primer3 binary to BIN.
                                                 Defaults to \$PATH/primer3_core
  -m, --mfeprimer_bin=BIN            (optional)  set absolute path of mfeprimer to BIN.
                                                 Defaults to \$PATH/MFEprimer.py
      --mfeprimer_args=ARGS          (optional)  use custom ARGS arguments for mfeprimer.
                                                 Defaults to "$mfeprimer_args_intern -e <10*PRODUCT_SIZE_RANGE_MAX>"

Miscellaneous:         
  -v, --verbose                      (optional)  print verbose status messages
  -h, --help                                     display this help and exit
  -V, --version                                  print version information and exit
  
Report bugs to: joffrey.fitz\@tuebingen.mpg.de
Contact for RBC questions: norman\@warthmann.com 
RBC home page: <http://rbc.weigelworld.org>

EOF

	exit 0;
}



sub print_params {
	warn "\n[INFO] Parameters:\n";
	warn "[INFO]     VERBOSE = $VERBOSE\n";
	warn "[INFO]     loc_seq_db = $loc_seq_db\n";
	warn "[INFO]     loc_rbci = $loc_rbci\n";
	warn "[INFO]     chr = $chr\n";
	warn "[INFO]     start = $start\n";
	warn "[INFO]     end = $end\n";
	warn "[INFO]     usr_filter = $usr_filter\n";
	warn "[INFO]     restriction_seq = $restriction_seq\n";

	warn "[INFO]     PRODUCT_SIZE_RANGE_MIN = $PRODUCT_SIZE_RANGE_MIN\n";
	warn "[INFO]     PRODUCT_SIZE_RANGE_MAX = $PRODUCT_SIZE_RANGE_MAX\n";
	warn "[INFO]     PRODUCT_OPT_SIZE = $PRODUCT_OPT_SIZE\n";

	# Regions for primers, this will result in
	warn "[INFO]     P_DIST = $P_DIST_MAX\n";
	warn "[INFO]     P_DIST_MIN = $P_DIST_MIN\n";

	# Flanking region around SNP
	warn "[INFO]     SNP_FLANKING = $SNP_FLANKING\n";

	# SNP annotation
	warn "[INFO]     Annotate primary SNP as major: $ann_major\n";
	# Variants stuff
	warn "[INFO]     MAX_SECONDARY_INDEL_LEN = $MAX_SECONDARY_INDEL_LEN\n";

	#
	# 3rd party tools
	####################################################################

	# Primer3 binary and config
	warn "[INFO]     p3_bin = $p3_bin\n";
	warn "[INFO]     p3_conf = $p3_conf\n";

	# MFEPrimer
	warn "[INFO]     mfeprimer_bin = $mfeprimer_bin\n";
	warn "[INFO]     mfeprimer_args = $mfeprimer_args\n";
	warn "\n";
}
