#!usr/bin/perl
use strict;

# ffn faa fna gff

if(-d "Adjustment"){
	system("rm -r Adjustment");
}
my %seq;
my $id;
my @title;

my @file=glob("./temp/filter/*.fna");
foreach my $one (@file){
	$one=~/filter\/(.*)\.fna/;
	my $id=$1;
	system ("perl src/Pre-Adjustment.pl $id ./temp/filter/$id.filter.ffn ./temp/filter/$id.filter.faa ./temp/filter/$id.fna ./temp/filter/$id.filter.gff");

	open IN , "./temp/filter/$id.filter.faa";
	while(<IN>){
		$_=~tr/\r\n//d;
		if(/>/){
			my @data=split(/\|/,$_);
			$id=$data[3];
			push @title , $id;
			next;
			}
		$seq{$id}.=$_;
		}
	close IN;
	}
my %cog;
my %list;

open COG , "myva" or die;
while(<COG>){
	$_=~tr/\r\n//d;
	if(/>(.*)/){
		$id=$1;
		next;
		}
	$cog{$id}.=$_;
	}
close COG;
open OUT , ">Adjustment/tmp/all.pep";

foreach(@title){
	print OUT ">$_\n$seq{$_}\n";
	}
	
foreach(keys %cog){
		print OUT ">$_\n$cog{$_}\n";
	}
close OUT;


