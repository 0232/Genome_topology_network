#!usr/bin/perl
use strict;

my $reff;
my $qry;
my $path;
my @file=glob("./temp/filter/*.fna");
my $one;
my $two;
my $id;

foreach $one (@file){
	$one=~/(.\/temp\/filter\/)(.*)\.fna/;
	$path="./temp/data/";
	$qry=$2;
	
	foreach $two (@file){
		if($one eq $two){
			next;
			}
		$two=~/.\/temp\/filter\/(.*)\.fna/;
		$reff=$1;
		$id=$reff."#".$qry;
		system ("nucmer --prefix=$path$id $two $one");
		system ("show-coords -rcl $path$id.delta > $path$id.coord");
		}
	}
