#!/usr/bin/perl -w

#
# VCF parser based on VCF.pm
#
# written 2012 by Joffrey Fitz <joffrey.fitz@tuebingen.mpg.de> and Joerg Hagmann
#
#
# Dependencies:
#	VCF.pm
#

use strict;
use warnings;

use Carp;
use Vcf;
use Getopt::Long;
use Data::Dumper;

# flags of the format string to process
my @flags_ok = qw(GT AD DP GQ GL);

# some pre-defined annotation marks (in case VCF format changes some day)
#my $filter_pass_annotation = "PASS";
my $filter_pass_annotation = ".";

# RBCI is tab delimited, columns 11+ contain the sample strings
my %rbci_h = (
	0 => "CHROM",
	1 => "POS",
	2 => "ID",
	3 => "REF",
	4 => "ALT",
	5 => "MAJOR",
	6 => "MINOR",
	7 => "QUAL",
	8 => "FILTER",
	9 => "TYPE",
	10 => "SEGR",
	11 => "HAPSC",
	12 => "FORMAT",
);


# initialize:
my $cmdline = $0." ".join(" ",@ARGV);
my %CMD;
my $vcf_file	= "";
my $ref_sample	= "";
my $pool1	= "";
my $pool2	= "";
my $parent1	= "";
my $parent2	= "";
my $pcr		= 0;
my %samples	= ();
# pcr vars:
my $sum_DPs=0;
my $nr_DPs=0;

# get user options:
GetCom();

my %pools = ("$pool1" => 1, "$pool2" => 1);
my %parents = ("$parent1" => 1, "$parent2" => 1);


# Parse VCF header
open F, "<$vcf_file";
my $vcf = Vcf->new(fh=>\*F);

# Print header
print "##fileformat=RCBIv0.1\n";
while (<F>) {
	print if (m/^##/);
}
print "##".$cmdline."\n";
seek F,0,0;

$vcf->parse_header();

# Get all samples available
my (@samples) = $vcf->get_samples();

# Parse through samples, check if the user-specified are present 
# and store them in a hash with (name,index) as (key,value)
foreach my $sample (@samples) {
	if (exists $samples{$sample}) {
		$samples{$sample} = $$vcf{has_column}{$sample} - 1;
	}
}

# Sanity check: number of samples must be 2 or 4!
if ( 	($parent1 ne "" && scalar(keys %samples) != 4) ||
	($parent1 eq "" && scalar(keys %samples) != 2) )
{
	die "[Error] Number of samples incorrect! Could not recognize all sample names. Sample names should be equal to the ones in the vcf file\n";
}


# Loop VCF lines
my $header_printed=0;
my $linenr=0;
while (my $x=$vcf->next_data_hash())
{
	if ( !$header_printed ) 
	{
		my @hf;

		# print rbci headers
		foreach my $h (sort{$a<=>$b} keys %rbci_h) {
			if ($pcr && $h>6) {
				push @hf, "FILTER";
				last;
			}
			push @hf, $rbci_h{$h};
		}

		# This prints header for sample columns
		foreach my $sample_id (sort {$samples{$a}<=>$samples{$b}} keys %samples)
		{
			if (!$pcr) {
				push @hf, $sample_id;
			}
			else {
				push @hf, "AD_MAJOR_$sample_id"; push @hf, "AD_MINOR_$sample_id"; push @hf, "DP_$sample_id";
				push @hf, "AD_MAJOR_$sample_id/DP_$sample_id"; push @hf, "AD_MINOR_$sample_id/DP_$sample_id";
			}
		}
		print "#" . join("\t", @hf) . "\n";

		$header_printed = 1;
	}

	++$linenr;
	
	my %rbci = (
		"CHROM"	=> $$x{CHROM},			# chromosome
		"POS"	=> $$x{POS},			# position
		"ID"	=> ".",				# sample name
		"REF"	=> $$x{REF},			# reference allele
		"ALT"	=> join(",", @{$$x{ALT}}),	# alternative allele
		"MAJOR"	=> $$x{REF},			# major allele
		"MINOR"	=> join(",", @{$$x{ALT}}),	# minor allele
		"QUAL"	=> $$x{QUAL},			# quality
		"FILTER"=> "", 				# filter
		"TYPE"	=> "", 				# type
		"SEGR"	=> "", 				# segr
		"HAPSC"	=> 0,				# haplotype score
		"FORMAT"=> join(":", @flags_ok),	# format string
		#sample_ID1, sample_ID2,... =>		# samples ...
	);

	my @segr=(0,0);
	my $filter = \@{$$x{FILTER}};
	my @parent_alleles;
	my $pool_allinfo = 1;



	# check presence of tags from @flags_ok
	foreach my $f (@flags_ok) {
		my $found = 0;
		foreach my $t (@{$$x{FORMAT}}) {
			$found = 1 if ($t eq $f);
		}
		if (!$found) {
			die "[Error] Required tag field $f does not exist in input entry number $linenr! Aborted";
		}
	}


	# Parse read counts of reference sample (to determine if alleles should be switched)
	my $switch = 0;
	my ($alleles_ref, $seps_ref, $is_phased_ref, $is_empty_ref) = $vcf->parse_haplotype($x,"$ref_sample");

	if (!$is_empty_ref && exists $$x{gtypes}{$ref_sample}{AD}) {
		my ($one, $two) = split ",", $$x{gtypes}{$ref_sample}{AD};
		if ($one < $two) {
			$switch = 1;
	                $rbci{"MAJOR"} = $rbci{"ALT"};
        	        $rbci{"MINOR"} = $rbci{"REF"};
		}
	}


	# compute type based on ref/alt observation:
	# indel, SNP, multi-allelic
	my $alt_str = join(",", @{$$x{ALT}});
	if( $alt_str =~ m/,/ || $$x{POS} =~m/,/)	{
		$rbci{"TYPE"} = "multi-allelic";
                $rbci{"MAJOR"} = "NA";
                $rbci{"MINOR"} = "NA";
	} elsif(length($alt_str)>=2 || length($$x{REF})>=2) {
		$rbci{"TYPE"} = "indel";
	} elsif($alt_str ne $$x{REF}) {
		$rbci{"TYPE"} = "SNP";
	} else {
		$rbci{"TYPE"} = ".";
	}
	

	# get haplotype score:
	$rbci{"HAPSC"} = (exists $$x{INFO}{HaplotypeScore}) ? $$x{INFO}{HaplotypeScore} : "NA";


	# filter absence of SB tag
	if (!exists $$x{INFO}{SB}) {
		if ($$filter[0] eq $filter_pass_annotation) {
			shift @$filter;
		}
		push @$filter, "NoSBtag";
	}


	# Loop through all available samples
	foreach my $sample (sort {$samples{$a}<=>$samples{$b}} keys %samples) {

		my ($alleles, $seps, $is_phased, $is_empty) = $vcf->parse_haplotype($x,"$sample");

		# if sample is "./."
		if ($is_empty) {
			$pool_allinfo = 0 if (exists $pools{$sample});
			$rbci{$sample} = "";
			foreach (keys %{$$x{gtypes}{$sample}}) {
				$rbci{$sample} .= $$x{gtypes}{$sample}{$_};
			}
			next;
		}


		# get phases (0/1, 0|1, 1|1, ...):

		# dirty workaround: if $$seps[0] equals '|' it has to be escaped
		my $splitchar = $$seps[0] eq "|" ? "\\\|" : $$seps[0];
		my ($phase1, $phase2) = split(/$splitchar/, $$x{gtypes}{$sample}{GT});

		# filter for SUS homozygotes
		if (exists $pools{$sample} && $sample ne $ref_sample && $phase1 eq $phase2) {
			if ($$filter[0] eq $filter_pass_annotation) {
				shift @$filter;
			}
			push @$filter, "NonPhasingPoolHomozygous";
		}

		# filter parental allele information
		if (exists $parents{$sample}) {
			# filter for parent heterozygotes
			if ($phase1 ne $phase2) {
				if ($$filter[0] eq $filter_pass_annotation) {
					shift @$filter;
				}
				push @$filter, "ParentHeterozygous";
			}
			else {
				# store parental alleles
				push @parent_alleles, $phase1;
			}
		}


		# get allele counts (AD):
		my ($allele_count1, $allele_count2) = split ",", $$x{gtypes}{$sample}{AD};
		my ($minor, $major) = $switch ? ($allele_count1, $allele_count2) : ($allele_count2, $allele_count1);

		# switch all necessary stuff:
		if ($switch && $rbci{"TYPE"} ne "multi-allelic") {

			# switch phases (GT) if applicable
			if (!($phase1 ne $phase2 && $$seps[0] eq "/")) {
				$$x{gtypes}{$sample}{GT} =
					($phase1 == 0 ? 1 : 0) . $$seps[0] .
					($phase2 == 0 ? 1 : 0);
			}

			# switch allele counts (AD) if applicable
			$$x{gtypes}{$sample}{AD} = $allele_count2.",".$allele_count1;

			# switch PL/GL flag
			my ($one, $two, $three) = split ",", $$x{gtypes}{$sample}{GL};
			$$x{gtypes}{$sample}{GL} = $three.",".$two.",".$one;

		}

		# store allele counts
		if (exists $pools{$sample}) {
			$segr[0] += $major;
			$segr[1] += $minor;
		}


		# assemble format string for each sample
		$rbci{$sample} = "";
		if (!$pcr) {
			foreach my $flag (@flags_ok) {
				$rbci{$sample} .= ":" if ($rbci{$sample} ne "");
				if (exists $$x{gtypes}{$sample}{$flag}) {
					$rbci{$sample} .= $$x{gtypes}{$sample}{$flag};
				}
				else {
					$rbci{$sample} .= "na";
				}
			}
		}
		else {
			$rbci{$sample} = $major . "\t" . $minor . "\t" . $$x{gtypes}{$sample}{DP};
			$rbci{$sample} .= "\t" . ($major/$$x{gtypes}{$sample}{DP}) . "\t" . ($minor/$$x{gtypes}{$sample}{DP}) if ($$x{gtypes}{$sample}{DP} > 0);
			$sum_DPs += $$x{gtypes}{$sample}{DP};
			$nr_DPs++ if ($$x{gtypes}{$sample}{DP} > 0);
		}

	} # for each sample


	# filter for segregating parents
	if (scalar(@parent_alleles) == 2 && $parent_alleles[0] eq $parent_alleles[1]) {
		if ($$filter[0] eq $filter_pass_annotation) {
			shift @$filter;
		}
		push @$filter, "ParentsNotSegr";
	}


	# Add segregation info, if any
	if($pool_allinfo) {
		$rbci{"SEGR"} = join ",", @segr;
	} else {
		$rbci{"SEGR"} = "NA";
	}
	
	# Add filter
	$rbci{"FILTER"} = join(";", @$filter);

	# Print output
	my @out;
	foreach my $c ( sort {$a<=>$b} keys %rbci_h) {
		if ($pcr && $c>6) {
			push @out, ".";
			last;
		}
		push @out, $rbci{$rbci_h{$c}};
	}
	foreach my $sample (sort {$samples{$a}<=>$samples{$b}} keys %samples) {
		push @out, $rbci{$sample};
	}
	print join("\t", @out);
	print "\n";
}

if ($pcr) {
	print STDERR "Average DP: ".($sum_DPs/$nr_DPs)."\n";
}



### Read command line parameters
sub GetCom {
	my @usage = ("\nUsage: $0 [options]

VCF to RBCI converter.

Mandatory arguments:
--vcf				STRING	Input VCF file
--phasing_sample		STRING	sample which will be used for phasing (sample name)
--phenotype_plus_pool		STRING	sample pool which shows the phenotype (sample name)
--phenotype_minus_pool		STRING	sample pool which does not show the phenotype (sample name) 

Optional arguments:
--phenotype_plus_parent		STRING	Parent which shows the phenotype (sample name)
--phenotype_minus_parent	STRING	Parent which does not show the phenotype (sample name)
--pcr					Analyze PCR products
\n");

	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "vcf=s", "phasing_sample=s", "phenotype_plus_pool=s", "phenotype_minus_pool=s", "phenotype_plus_parent=s", "phenotype_minus_parent=s", "pcr");

	die("Please specify a vcf file, a reference sample and both pool samples\n") unless (defined $CMD{vcf} && defined $CMD{phasing_sample} && defined $CMD{phenotype_plus_pool} && defined $CMD{phenotype_minus_pool});

	if (defined $CMD{phenotype_plus_parent} || defined $CMD{phenotype_minus_parent}) {
		die("Please specify both parents\n") unless (defined $CMD{phenotype_plus_parent} && defined $CMD{phenotype_minus_parent});
	}

	$vcf_file	= $CMD{vcf};
	$ref_sample 	= $CMD{phasing_sample};
	$pool1		= $CMD{phenotype_plus_pool};
	$pool2		= $CMD{phenotype_minus_pool};
	$samples{$pool1} = 0;
	$samples{$pool2} = 0;
	$pcr = 1 if (defined $CMD{pcr});

	if (defined $CMD{phenotype_plus_parent}) {
		$parent1 = $CMD{phenotype_plus_parent};
		$parent2 = $CMD{phenotype_minus_parent};
		$samples{$parent1} = 0;
		$samples{$parent2} = 0;
	}
}
