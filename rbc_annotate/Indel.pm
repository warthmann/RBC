#!/usr/bin/perl

# --------------------------------------------------------------------
# ShoreMap extension to SHORE:
# Identification of causal mutations using IlluminaGA2 sequencing data
#
# Written by Stephan Ossowski and Korbinian Schneeberger
# --------------------------------------------------------------------

use strict;
use warnings;

package Indel;

sub new {
	my $self = {
		ecotype       => '',
		chromosome    => 0,
		position      => 0,
		end           => 0,
		global_start  => 0,
		global_end    => 0,
		vartype       => '',
		ref_base      => '',
		new_base      => '',
		seq           => '',
		stype         => "",
		gene_id       => 'NA',
		gene_pos      => 0,
		cds_pos       => 0,
		cds_pos_left  => 0,
		cds_pos_right => 0,
		codon_pos     => 0,
		new_stop      => 0,
		lost_start    => 0,
		lost_stop     => 0,
		splicechange  => 0,
		frameshift    => 0,
		alt_coding    => "",
		#support       => 0,
		#concordance   => 0,
		#quality       => 0,
		peak_distance => 999999999,
		fields	      => [],
	};
	bless $self;
	return $self;
}

sub init
{
	#my ($self, $ecotype, $chromosome, $global_begin, $global_end, $begin, $end, $seq, $support, $concordance, $quality, $peak_pos, $vartype) = @_;
	my ($self, $ecotype, $chromosome, $global_begin, $global_end, $begin, $end, $refbase, $newbase, $seq, $vartype, $peak_pos, $fields) = @_;
	$self->{ecotype}     = $ecotype;
	$self->{chromosome}  = $chromosome;
	$self->{global_start}= $global_begin;
	$self->{global_end}  = $global_end;
	$self->{position}    = $begin;
	$self->{end}         = $end;
	$self->{ref_base}    = $refbase;
	$self->{new_base}    = $newbase;
	$self->{seq}         = $seq;
	#$self->{support}     = $support;
	#$self->{concordance} = $concordance;
	#$self->{quality}     = $quality;
	$self->{peak_distance} = abs($peak_pos - $begin);
	$self->{vartype}     = $vartype;

	foreach (@{$fields}) {
		push @{$self->{fields}}, $_;
	}

}

1;
