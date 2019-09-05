#!/usr/bin/perl -w
use warnings;
use strict;
use iGraph;
 
die "perl $0 <1:get mobile COGs 2:get distance matrix 3:get boostrap tree[newick format]> <output file name>" unless ($#ARGV == 1);

my $cog_dir = "./Adjustment/cog/";
my $gff_dir = "./Adjustment/gff/";
my $all = iGraph->new($cog_dir,$gff_dir,'gff');




if($ARGV[0] == 1)
{
	my %hash_info = &BuildCOGInfohash();
	open OUT,">$ARGV[1]";
	print OUT 'COG ID',"\tAnnotation","\tDifferent Degree\n";
	my $mobile = $all->getMobile(@{$all->getSpecies()});
	foreach my $key (sort{${$mobile}{$b} <=> ${$mobile}{$a}} keys %{$mobile}) 
	{
		print OUT $key,"\t",$hash_info{$key},"\t",${$mobile}{$key},"\n";
	}
	
}
if($ARGV[0] == 2)
{
	open OUT,">$ARGV[1]";
	print OUT $all->getDistanceMatrix("mega","SimiMatrixUn");
}

if($ARGV[0] == 3)
{
	$all->getBootstrapMatrix(1000,$ARGV[1]);
}

sub BuildCOGInfohash()
{
	open COG,"<whog";
	my %coginfo_hash =();
	while(my $line = <COG>)
	{
		if($line =~ /COG\d{4}/)
		{
			chomp $line;
			my @arr = split/\s+/,$line,3;
			$coginfo_hash{$arr[1]} = $arr[2];
		}
		else
		{
			next;
		}	
	}
	return %coginfo_hash;
}

