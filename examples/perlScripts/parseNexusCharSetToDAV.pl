#!/usr/bin/perl

# Reads a nexus file non-interleaved and with charset commands for a supermatrix.
# Writes a set of subsets of taxa according to the charset commands, and optionally some other stats

# USAGE: parseNexusCharSetToDAV.pl -f nexusfile [-v]

# The -v option will print out a bunch of other stats on the supermatrix

#############################

$name = '[\w\d\_\.]+|\'.*?\''; # notice the non-greedy match needed in the '...' format

use Getopt::Long;
my $nexFile;
my $verbose=0;
$result = GetOptions ("f=s"   => \$nexFile,      # string
                     "v"  => \$verbose);  # flag



open FH1, "<$nexFile";


while (<FH1>)
	{
	if (/ntax\s*=\s*(\d+)/i)
		{$nTax=$1}
	elsif (/nchar\s*=\s*(\d+)/i)
		{$nChar=$1}
	elsif (/matrix$/i) {$okToRead=1}
	elsif ($okToRead==1 && /;/) {$okToRead=0}
	elsif ($okToRead && /($name)\s+([A-Z\?\-]+)/)
		{
		push @taxa,$1;
		push @seqs,$2;
		#print "$1:$2\n";
		}
	elsif (/charset\s+(\w+)\s*=\s*(\d+)\s*\-\s*(\d+)\s*;/i)
		{
		push @locus,$1;
		push @start,$2;
		push @end, $3;
		if ($verbose) {print "charset $1:$2 - $3\n"}
		}

	}
if ($verbose) {print "nTax=$nTax nChar=$nChar\n\n"}
$nTaxArray=@taxa;
$nLoci = @locus;
die "Number of taxa in matrix ($nTaxArray) does not match number specified in header\n" if ($nTaxArray != $nTax);
for $tax (0..$nTax-1)
	{
	$seq = $seqs[$tax];
	for $loc (0..$nLoci-1)
		{
		$slength = $end[$loc]-$start[$loc]+1;
		$subseq = substr($seq,$start[$loc]-1,$slength);
		$DAV[$tax][$loc]=present($subseq);
		}
	}
	
if ($verbose)	
{
for $tax (0..$nTax-1)
	{
	printf ("%-20s\t",$taxa[$tax]); # left justifies taxon name in large enough field to display nicely
	$sumt=0;
	for $loc (0..$nLoci-1)
		{
		$dav=$DAV[$tax][$loc];
		print "$dav";
		$sumloci[$loc]+=$dav;
		$sumt+=$dav;
		}
	$sumtax[$tax]=$sumt;
	if ($sumtax[$tax]==$nLoci) 
		{ push @refTaxa,$taxa[$tax]  };
	print "\t\t$sumtax[$tax]\n";	
	}
print "Reference taxa in this matrix: @refTaxa\n";
print "\n\nNumber of taxa present per locus:\n";
$sumall=0;
for $loc (0..$nLoci-1)
	{
	print "$locus[$loc]\t$sumloci[$loc]\n";
	$sumall+=$sumloci[$loc];
	}
$coverage = $sumall/($nLoci*$nTax);
print "Coverage density is $coverage\n";

}

for $loc (0..$nLoci-1)
	{
	for $tax (0..$nTax-1)
		{
		$dav=$DAV[$tax][$loc];
		if ($dav ==1) {print "$taxa[$tax] "}
		}
	print "\n";
	}
#############################################################################################

sub present

# Determines if the sequence has ANY instances of the defined character class

{
my ($seq)=@_;
$DNAclass_noambig = "ACGT";
$PROTclass_noambig = "A-WYZ";

if ($seq =~ /[$DNAclass_noambig]/i)
#if ($seq =~ /[$PROTclass_noambig]/i)
	{return 1}
else
	{return 0}
}

