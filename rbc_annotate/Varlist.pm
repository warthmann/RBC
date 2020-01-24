#!/usr/bin/perl 

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;
use SNP;
use Indel;

package Varlist;

### Constructor
sub new {
	my $self = {
		chromosome   => "",
		region_start => 0,
		region_end   => 0,
		var_file     => "",
		vars	     => [],
		maxend	     => 0,
	};
	bless $self;
	return $self;
}
		

sub get_and_prioritize {

	### Init
	my ($self, $chromosome, $region_start, $region_end, $var_file, $peak_pos, $filepos, $vcf, $majmin_allele, $filter) = @_;
	$self->{chromosome}   = $chromosome;
	$self->{region_start} = $region_start;
	$self->{region_end}   = $region_end;
	$self->{var_file}     = $var_file;

	my @filter = @{$filter};

	open VARFILE, $var_file or die "Cannot open variation file: " . $var_file . "\n";
	seek VARFILE, $filepos, 0;
	while (<VARFILE>) {
		my @a = split/\t/;

		my ($type, $sample, $chr, $position, $end, $refbase, $snp, $length, $seq);

		chomp;
		my @fields = split "\t";

		if (substr($_,0,1) eq "S") {
			### SNP
			$type = $fields[0];
			$sample = $vcf ? $fields[3] : $fields[1];
			$chr = $fields[2-$vcf];
			$position = $fields[3-$vcf];
			$refbase = $fields[4];
			$snp = $fields[5];
			@fields = @fields[6..(scalar(@fields)-1)];

			#if ($snp =~ m/,/) {
			#	my @bases = split/,/,$snp;
			#	$snp = $bases[0];
			#}

			if ($snp eq "") { die "Unknown SNP base on chr $chr at pos $position!\n"; }
			chomp $snp;
		}
		else {
			### Indel
			if (!$vcf) {
				$type = $fields[0];
        	                $sample = $fields[1];
	                        $chr = $fields[2];
        	                $position = $fields[3];
                	        $end = $fields[4];
                        	$length = $fields[5];
				$seq = $fields[6];
        	                @fields = @fields[7..(scalar(@fields)-1)];
			}
			else {
				$type = $fields[0];
				$chr = $fields[1];
				$position = $fields[2];
				$sample = $fields[3];
				$refbase = $fields[4];
				$snp = $fields[5];
				@fields = @fields[6..(scalar(@fields)-1)];

				$length = length($refbase) > length($snp) ? length($refbase)-1 : length($snp)-1;
				$seq    = length($refbase) > length($snp) ? substr($refbase,1) : substr($snp,1);
				$end    = $position+$length-1 if ($type eq "D");
				$peak_pos = 1;
			}

			if ($type eq "I") { $end = $position+1; }
		}

		if ($sample eq ".") {
			if ($snp =~ m/,/) {
				$sample = "na";
			} else {
				$sample = "minor" if ($fields[1] =~ m/$snp/);
				$sample = "major" if ($fields[0] =~ m/$snp/);
			}
		}

		if ( $region_start == -1 || ($chr eq $chromosome && $position >= $region_start && $position <= $region_end) ) {
	
			last if ($chromosome ne "-1" && $chromosome ne $chr);

			# discard alleles which are not desired (major, minor or all)
			next if ($majmin_allele eq "major" && $snp eq $fields[1]);
			next if ($majmin_allele eq "minor" && $snp eq $fields[0]);

			# discard SNPs with given filters:
			my $filtered=0;
			foreach (@filter) {
				$filtered = 1 if ($fields[3] =~ m/$_/);
			}
			next if ($filtered);

			if ($type eq "S") {
				my @alt_bases = split/,/, $snp;
				foreach (@alt_bases) {
					$snp = $_;
					my $newsnp = new SNP();
					push @{$self->{vars}}, $newsnp;

					$newsnp->init( $sample, $chr, $position, $refbase, $snp, $type, $peak_pos, \@fields);
					$self->{maxend} = $self->{maxend} < $position ? $position : $self->{maxend};
				}
			}
			else {
				my $newindel = new Indel();
				push @{$self->{vars}}, $newindel;

				$newindel->init( $sample, $chr, $position, $end, $position, $end, $refbase, $snp, $seq, $type, $peak_pos, \@fields);
				$self->{maxend} = $self->{maxend} < $position ? $position : $self->{maxend}; 
			}

			$chromosome = $chr;

			$filepos = tell VARFILE;
		}
	}
	close VARFILE;

	return $filepos;
}

1;
