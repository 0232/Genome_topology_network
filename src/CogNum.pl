#!usr/bin/perl
use strict;

my @file=glob("Adjustment/topo/*.topo");
my @data;


foreach my $one (@file){
	$one=~/topo\/(.*)\.topo/;
	my $name=$1;
	my %hash=();
	
	open IN , "$one";
	while(<IN>){
		@data=split(/\t/,$_);
		$hash{$data[0]}++;
		}
	close IN;
	
	open OUT , ">Adjustment/topo/$name.cog.num";
	foreach my $two (keys %hash){
		print OUT "$two\t$hash{$two}\n";
		}
	close OUT;
	}