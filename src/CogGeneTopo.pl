#!usr/bin/perl
use strict;

my @file=glob("Adjustment/gff/*.gff");
my @data;
my $pro_id;
my $id;
my %cog;

mkdir "Adjustment/topo" unless -d "Adjustment/topo";
open COG , "./Adjustment/cog/cog.group.txt" or die;
while(<COG>){
	$_=~tr/\r\n//d;
	@data=split(/\s/,$_);
	my $temp=shift @data;
	$temp=~s/\://;
	
	foreach my $one (@data){
		$one=~/\|(.*)/;
		$cog{$1}=$temp;
		}
	}
close COG;


foreach my $one (@file){
	$one=~/gff\/(.*)\.gff/;
	my $spe=$1;
	open OUT , ">Adjustment/topo/$spe.topo";

	open GFF , "$one";
	while(<GFF>){
		@data=split(/\t/,$_);
		if(/protein_id\=(.*?)\;/){
			$_=~/protein_id\=(.*?)\;/;
			$pro_id=$1;
			$pro_id=~/__(.*)/;
			my $spe=$1;
			$id=$spe.'|'.$pro_id;
			if(!$cog{$pro_id}){
				next;
				}
			print OUT "$cog{$pro_id}\t$id\n";
			}
		}
	close GFF;
	close OUT;
	}