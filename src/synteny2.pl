#!usr/bin/perl
use strict;

#查找每个菌株的共有的one hit坐标，过滤，输出在坐标上的基因的GFF "synteny_once_gene.txt"  只有gene无CDS ,输出唯一共线性坐标 “synteny_once_coordinate.txt”

my %synteny;
my %count;
my $id;
my %once;


my @fna=glob("./temp/filter/*.fna");
my $count=@fna-1;

open IN , "./temp/data/synteny_coordinate.txt" or die;
open OUT , ">./temp/data/synteny_once_gene.txt";
open COOR , ">./temp/data/synteny_once_coordinate.txt";
while(<IN>){
	$_=~tr/\r\n//d;
	if(/>(.*)/){
		$id=$1;
		$count{$id}++;
		next;
		}
	$synteny{$id}.="$_\t";
	}
close IN;

my @data;
my $one;
my $two;
my $three;
my $four;
my $xx;
my $i;
foreach $one (keys %synteny){
	@data=split(/\t/,$synteny{$one});
	my %hash=();
	my @coor=();
	
	foreach $two (@data){
		$two=~/(\d+)\.\.(\d+)/;

		foreach $three ($1..$2){
			$hash{$three}++;
			}
		}
	
	foreach $four (keys %hash){
		if($hash{$four}==$count){
      push @coor , $four;
			}
		}
	print COOR ">$one\n";
	@coor= sort {$a <=> $b} @coor;
	$xx=$coor[0];
	for($i=1;$i<=$#coor;$i++){
	 if($coor[$i]-$coor[$i-1]==1){
	  next;
	  }
	 if($coor[$i]-$coor[$i-1]>1){
	  if($xx and $coor[$i-1]){
	  	print COOR "$xx..$coor[$i-1]\t";
	  	$once{$one}.="$xx..$coor[$i-1]\t";
	  	}	  		
	  $xx=$coor[$i];
	  next;
	  }
	 }
	 if($xx and !($xx==$coor[-1])){
	   print COOR "$xx..$coor[-1]\n";
	   $once{$one}.="$xx..$coor[-1]\t";
	  }
	}
close COOR;


my %gff;
my %gene;
my @file=glob("./temp/data/*.gff");
my $name;

foreach $one (@file){
	$one=~/data\/(.*)\.gff/;
	$name=$1."#";
	open IN , "$one";
	while(<IN>){
		if(/#/){
			next;
			}
		@data=split(/\t/,$_);
		if($data[2]=~/gene/){
			$id=$data[2];
			my $sca=$name.$data[0];
			$gene{$sca}.="$_";
			}
		$gff{$id}.=$_;
		}
	close IN;
	}


my $xx;
my $yy;
my $aa;
my $bb;
my $gene_id;
foreach $one (keys %gene){
	my @gene=split(/\n/,$gene{$one});
	my @coor=split(/\t/,$once{$one});
	
	foreach $two (@gene){
		if($two=~/old_locus_tag/){
			$two=~/locus_tag\=(.*?)\;old/;
		  $gene_id=$1;
			}
		if(!($two=~/old_locus_tag/)){
			$two=~/locus_tag\=(.*?);{0,1}/;
		  $gene_id=$1;
			}
		@data=split(/\t/,$two);
		$xx=$data[3];
		$yy=$data[4];
		
		foreach $three (@coor){
			$three=~/(\d+)\.\.(\d+)/;
			$aa=$1;
			$bb=$2;
			if($xx>=$aa and $yy<=$bb){
				print OUT "$two\n";

				}
			}
		}
	}
close OUT;


