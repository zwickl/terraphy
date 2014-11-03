#!/usr/bin/perl

# Writes a list of all taxa in the input trees in the first row, 
# and a list of all triplets after that, one per row, where the 
# third field is the outgroup to the first two.


# USAGE: maketriplets.pl treefile

# Takes one file as a command line argument, which contains newick strings
# ASSUMES TREES ARE ALL BINARY! Dies otherwise

#######


open FH, "<$ARGV[0]";

while (<FH>)
	{
	chomp;
	$nwk =$_;
	$t = treeInit($nwk);
	push @trees,$t;	
	}
close FH;
%allTaxaH = allTaxa(@trees);
@allTaxaA = keys %allTaxaH;
print "@allTaxaA\n";
foreach $t (@trees)
	{
	makeTriplets($t->{ROOT});
	}

####################################################################
sub makeTriplets
	{
	my ($n) = @_;
	if (isLeaf($n)) {return};
	if (!isRoot($n)) 
		{
		my ($anc, $sib,$L0,$L1,$L2);
		($L1,$L2)=invLCA($n);
		$anc = $n->{ANC};
		if (  ${$anc->{DESC}}[0]== $n  ) 
			{    $sib = ${$anc->{DESC}}[1]    }
		else
			{    $sib = ${$anc->{DESC}}[0]    }  # Assumes a binary tree!
		$L0 = oneDesc($sib);
		print "$L1->{NAME}\t$L2->{NAME}\t$L0->{NAME}\n";
		}	
	for my $child (@{$n->{DESC}})
		{ makeTriplets($child); }
	}

####################################################################
sub treeInit
{
my ($newick)=@_;
my ($root,%treeH);
$name = '[\w\d\_\.]+|\'.*?\''; # notice the non-greedy match needed in the '...' format
@tokens = ($newick=~/($name|\,|\(|\)|\;)/g);
$tokens_ix=0;
parseCheck();
$tokens_ix=0;
#print "@tokens\n";

my $tok = next_tok();
if ($tok =~ /\(/)
	{
	$root = make_group();
	}
else
	{die "First token is incorrect: $tok\n";}

initIDTree(0,$root);
initNLEAVES($root);
my @tt = descLeafNodes($root);
my %leafH = makeLeafHash(@tt);
#foreach (keys %leafH) {print "$_\n"}
$treeH{ROOT}=$root;
$treeH{LEAVES}=[@tt];
$treeH{LEAFH}={%leafH};
$treeH{NLEAVES}= $root->{NLEAVES};
die ("Fatal error: tree is NOT binary\n") if (!isBinaryTree($root));
return \%treeH;
}

sub initIDTree
{
my ($startIX,$root)=@_;
$gIX=$startIX;
recurseIndexTree($root);
}
sub recurseIndexTree
{
my ($nodeRef)=@_;
$nodeRef->{ID}=$gIX++;
for $child (@{$nodeRef->{DESC}})
	{ recurseIndexTree($child); }
return ;
}
sub initNLEAVES
{
my ($nodeRef)=@_;
my $sum;
if (isLeaf($nodeRef)) {$sum=1} else {$sum=0}; 
for $child (@{$nodeRef->{DESC}})
	{ $sum += initNLEAVES($child); }
$nodeRef->{NLEAVES}=$sum;
return $sum;
}

sub oneDesc
# returns a single leaf descendant of a given internal node
{
my ($nd)=@_;
while (!isLeaf($nd))
        { $nd = ${$nd->{DESC}}[0]; }
return ($nd);

}


sub invLCA
# returns two leaf nodes of given internal node, such that the node is the lca of the two leaves
{
my ($nd)=@_;
if (isLeaf($nd)) {return (undef,undef)}
$n1 = ${$nd->{DESC}}[0];
$n2 = ${$nd->{DESC}}[1];

while (!isLeaf($n1))
	{ $n1 = ${$n1->{DESC}}[0]; }
while (!isLeaf($n2))
	{ $n2 = ${$n2->{DESC}}[0]; }
return ($n1,$n2);

}

####################################################################
sub makeLeafHash
{
my (%leafH)=();
foreach (@_)
	{$leafH{$_->{NAME}}=$_}
return %leafH;
}
####################################################################
sub descLeafNames
{
return map {$_->{NAME}} descLeafNodes(@_[0]);
}

sub descLeafNodes
{
my ($root)=@_;
undef @taxaT;
recursePush($root);
return @taxaT;
}

sub recursePush
{
my ($nodeRef)=@_;

if (isLeaf($nodeRef))
	 {push @taxaT, $nodeRef;}
for $child (@{$nodeRef->{DESC}})
	{ recursePush($child); }
return ;
}



####################################################################
sub make_group
{
my $rootRef=nodeNew("");
while (my $tok = next_tok())
	{
	if ($tok =~ /$name/)
		{
		my $nodeRef = nodeNew($tok);
		addChild($rootRef,$nodeRef);
		}
	elsif ($tok =~ /\(/)
		{
		$nodeRef = make_group();
		addChild($rootRef,$nodeRef);
		}
	elsif ($tok =~ /\)/)
		{
		return $rootRef;
		}
	elsif ($tok =~ /,/)
		{
		next;
		}
	}
}


# **********************************************************

sub next_tok
{
if ($tokens_ix >= $#tokens) {return 0}
return ($tokens[$tokens_ix++])
}

# **********************************************************

sub nodeNew
{
my ($name)=@_;
return {ID=>-1,NAME=>$name,DESC=>[],ANC=>-1,NLEAVES=>-1};
}
# **********************************************************
sub addChild 
{
my ($nodeRef,$childRef)=@_;
$childRef->{ANC}=$nodeRef;
push @{ ${$nodeRef}{DESC} },$childRef;
}


# **********************************************************
sub isBinaryTree
{
my ($nodeRef)=@_;
if (isLeaf($nodeRef)) {return 1};
if (scalar @{$nodeRef->{DESC}} != 2) {return 0};
for $child (@{$nodeRef->{DESC}})
	{ if (!isBinaryTree($child)) {return 0}; }
return 1 ;
}
# **********************************************************

sub recursePrint
{
my ($nodeRef)=@_;
print "$nodeRef->{ID}: num leaves:$nodeRef->{NLEAVES}:taxon name = $nodeRef->{NAME}\n";
if (!isLeaf($nodeRef) && exists $nodeRef->{LCA1})
	{
	print "\tTree 1 invLCA: $nodeRef->{LCA1}{ID}\n";
	print "\tTree 2 invLCA: $nodeRef->{LCA2}{ID}\n";
	}
for $child (@{$nodeRef->{DESC}})
	{ recursePrint($child); }
return ;
}

# **********************************************************

sub isLeaf
{
my ($nodeRef)=@_;
if (scalar @{$nodeRef->{DESC}} == 0)
	{return 1}
else
	{return 0}
}
sub isRoot
{
my ($nodeRef)=@_;
if ($nodeRef->{ANC} == -1)
	{return 1}
else
	{return 0}
}

#**********************************************************
sub parseCheck
{
for $tok (@tokens)
	{
	if ($tok =~ /\(/) {$left++};
	if ($tok =~ /\)/) {$right++};
	if ($tok =~ /$name/) {$nTax++};
	if ($tok =~ /\,/) {$commas++};
	}
#print "ntaxa=$nTax,left=$left, right=$right,commas=$commas\n";
die "Unmatched parens in newick string\n" if ($left != $right);
}
#**********************************************************
sub allTaxa
{
my %allTaxaH;
foreach (@_)
        {
        my %leafH = %{$_->{LEAFH}};
        foreach my $t (keys %leafH)
                {$allTaxaH{$t}=1}
        }
return %allTaxaH
}


