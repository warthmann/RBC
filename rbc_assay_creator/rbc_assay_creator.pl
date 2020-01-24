#!/usr/bin/perl

###############################################################################
#
# rbc_assay_creator.pl - Converts .rbcp to various formats
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

use strict;
use warnings;

my $verbose = 0;
my $i = "";
my $o = ".";
my $p = "out"; # better: infile without suffix
my @b = "";
my $force = 0;
my ($print_version, $print_help);
my $prog_name = "rbc_assay_creator";

use Getopt::Long qw(:config no_ignore_case);;

my $result = GetOptions(
	"file|i=s"              	=> \$i,
	"output_dir|o=s"			=> \$o,
	"output_file_prefix|p=s"	=> \$p,
	"barcodes|b=s"              => \@b,
	"force|f"					=> \$force,
	"verbose|v"					=> \$verbose,
	"version|V"                 => \$print_version,
	"help|h"                    => \$print_help,
);

if(defined $print_version) {
	print_version();
}

if(defined $print_help) {
	print_help();
}

# Check mandatory arguments
if($i eq "") {
	print_help();
}
	
my $rbcp_out = "$o/$p.assays.rbcp";
my $shore_out = "$o/$p.assays.shore";
my $vcf_out = "$o/$p.assays.vcf";

# check if output directory exists
if(not -e $o) {
	die "Output directory '$o' does not exist.\n";
}

# check the barcodes
my $b_str = join(",",@b);
$b_str =~ s/^^,//;
# Check if files exists which must not be overwritten
if(not $force) {
	my ($rbcp_out_e, $shore_out_e, $vcf_out_e) = ("","","");
	if (-e $rbcp_out) {
		$rbcp_out_e = $rbcp_out;
	}
	if ($b_str ne "") {
	 	if(-e $shore_out) {
			$shore_out_e = $shore_out;
		}
	}
	if (-e $vcf_out) {
		$vcf_out_e = $vcf_out;
	}
	
	if($rbcp_out_e ne "" || $shore_out_e ne "" || $vcf_out_e ne "") {
		my $existing_files_list = join(", ", ($rbcp_out_e, $vcf_out_e, $shore_out_e));
		$existing_files_list =~ s/, ,/,/g; $existing_files_list =~ s/^, //g;
		die "Output files already exist:\n$existing_files_list\nUse --force to overwrite\n";
	}
}

# Check user given barcodes
my %barcodes_h;

if($b_str ne "") {
	$b_str =~ s/\s//g;
	warn "[INFO] Checking barcodes '$b_str'\n";
	my @b_a = split(",", $b_str);
	die "[ERROR] Odd number of barcodes/sample_no\n" if(scalar(@b_a)%2);
	%barcodes_h = @b_a;	
	foreach my $b_key (keys %barcodes_h) {
		warn "[INFO]     $b_key => " . $barcodes_h{$b_key} . "\n";
	}
}

# Create handles
open(I, "$i" || die "Can't open infile '$i':$i");
open(RBCP, ">$rbcp_out") || die "Can't open $rbcp_out: $!\n";
	warn "[INFO] RBCP output file '$rbcp_out' created\n" if($verbose);
open(SHORE, ">$shore_out") || die "Can't open $rbcp_out: $!\n";
	warn "[INFO] shore output file '$shore_out' created\n" if($verbose);
open(VCF, ">$vcf_out") || die "Can't open $vcf_out: $!\n";
	warn "[INFO] vcf output file '$vcf_out' created\n" if($verbose);

# print headers for output files
print VCF output_vcf_header();

# Loop input, generate output
warn "[INFO] processing input file '$i'\n" if($verbose);
while(<I>) {
	chomp;
	next if(/^#/);
	my @c = split("\t");
	next if $c[2] eq ".";
	
	print SHORE output_shore(@c);
	print RBCP output_rbcp(@c) . "\n";
	print VCF output_vcf(@c) . "\n";
}

close RBCP;
close SHORE;
close VCF;

exit 0;


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

Usage: $0 [OPTIONS] --file|-i input.rbcp

Creates assay files for given rbcp input.

Mandatory input files:
  -i, --file=FILE           (mandatory) use FILE as rbcp input file
  -o, --output_dir=DIR      (=.) write output files to directory DIR
  -b, --barcodes=CSV        (optional) if given, shore barcode file is written.
                                       Multiple barcodes are supported.
                                       Example: -b "barcode,sample_no"
Miscellaneous:
  -f ,--force               (optional)  overwrite existing output files         
  -v, --verbose             (optional)  print verbose status messages
  -h, --help                            display this help and exit
  -V, --version                         print version information and exit
  
Report bugs to: joffrey.fitz\@tuebingen.mpg.de
Contact for RBC questions: norman\@warthmann.com 
RBC home page: <http://rbc.weigelworld.org>

EOF

	exit 0;
}


# SHORE
sub output_shore {
	my @c = @_;
	my $buf = "";

	my $ID = $c[2];
	my ($primer_left, $primer_right) = @c[12,13];
	
	foreach my $b_key (keys %barcodes_h) {
		# primer_left / primer_right (group 1)
		$buf .= "1	*	$primer_left	5prime	1	$barcodes_h{$b_key}_$ID\n";
		$buf .= "3	*	$primer_right	5prime	1	$barcodes_h{$b_key}_$ID\n";
		$buf .= "2	$b_key	*	5prime	1	$barcodes_h{$b_key}_$ID\n";
	
		# primer_right / primer_left (group 2)
		$buf .= "1	*	$primer_right	5prime	2	$barcodes_h{$b_key}_$ID\n";
		$buf .= "3	*	$primer_left	5prime	2	$barcodes_h{$b_key}_$ID\n";
		$buf .= "2	$b_key	*	5prime	2	$barcodes_h{$b_key}_$ID\n";
	}
	return $buf;
}



# RBCP
sub output_rbcp {
	my @c = @_;
	
	return join("\t", @c);
}


# VCF header
sub output_vcf_header {
	my ($y, $m, $d) = (localtime(time))[5,4,3];
	$y += 1900;
	
	my $header =<<EOH;
##fileformat=VCFv4.1
##fileDate=${y}${m}${d}
##source=$prog_name $VERSION
##source_file=$i
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
EOH

	return $header
}


# minimal VCF
sub output_vcf {
	my @c = @_;
	
	return "$c[0]\t$c[1]\t$c[2]\t$c[3]\t$c[4]\t.\t$c[8]\t" . ".";
}

