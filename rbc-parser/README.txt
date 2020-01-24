VCF_WDIR=/ebio/abt6/joffrey/projects/rice/VCF

#
# Get VCFTools
#
wget http://sourceforge.net/projects/vcftools/files/vcftools_0.1.8a.tar.gz
tar xvzf vcftools_0.1.8a.tar.gz
cd vcftools_0.1.8a && make

# Setup P5
export PERL5LIB=/ebio/abt6/joffrey/projects/rice/VCF/vcftools_0.1.8a/perl

# Check 
$PERL5LIB/vcf-to-tab

#
# Paths
#
TDIR=/ebio/abt6/joffrey/projects/rice/
WDIR=/ebio/abt6_projects6/nobackup/joff/rice/



# Example VCF
VCFEX=/ebio/abt6_projects6/nobackup/Rice/asyl/bwa_mapping_results_010312/vcf_new/pools_on_9311.het01.SW.BOTH.phased.filterflagged.vcf

cat $VCFEX | ./vcf_parse_1st.pl