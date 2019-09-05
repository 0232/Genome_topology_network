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
	system ("perl src/Pre-Adjustment.pl $id ./temp/data/$id.ffn ./temp/data/$id.faa ./temp/filter/$id.fna ./temp/data/$id.gff");

	open IN , "./temp/data/$id.faa";
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


