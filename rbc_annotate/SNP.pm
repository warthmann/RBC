#!/usr/bin/perl

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;

package SNP;

sub new {
	my $self = {
		ecotype       => '',
		chromosome    => 0,
		position      => 0,
		stype         => "",
		gene_id       => 'NA',
		gene_pos      => 0,
		cds_pos       => 0,
		codon_pos     => 0,
		new_stop      => 0,
		lost_start    => 0,
		lost_stop     => 0,
		splicechange  => 0,
		frameshift    => 0,
		ref_base      => '',
		new_base      => '',
		ref_aa        => '',
		new_aa        => '',
		degeneracy    => 0,
		codon         => "",
		codon_alt     => "",
		#support       => 0,
		#concordance   => 0,
		#quality       => 0,
		peak_distance => 999999999,
		protein_ref   => "",
		protein_alt   => "",
		vartype       => "",
		fields	      => [],
	};
	bless $self;
	return $self;
}

sub init
{
	#my ($self, $ecotype, $chromosome, $position, $ref_base, $new_base, $support, $concordance, $quality, $peak_pos, $vartype) = @_;
	my ($self, $ecotype, $chromosome, $position, $ref_base, $new_base, $vartype, $peak_pos, $fields) = @_;
	$self->{ecotype}     = $ecotype;
	$self->{chromosome}  = $chromosome;
	$self->{position}    = $position;
	$self->{end}	     = $position;
	$self->{ref_base}    = $ref_base;
	$self->{new_base}    = $new_base;
	#$self->{support}     = $support;
	#$self->{concordance} = $concordance;
	#$self->{quality}     = $quality;
	$self->{peak_distance} = abs($position - $peak_pos);
	$self->{vartype}     = $vartype;

	foreach (@{$fields}) {
		push @{$self->{fields}}, $_;
	}
}

1;
