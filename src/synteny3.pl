#!usr/bin/perl
use strict;

##

my %filter;
open ONCE , "./temp/data/synteny_once_gene.txt" or die;
while(<ONCE>){
	$_=~tr/\r\n//d;
	$filter{$_}=1;
	}
close ONCE;




my @file=glob("./temp/data/*.gff");
my $one;
my $two;
my @data;
my @temp;
my $sca;
my $id;

foreach $one (@file){
	$one=~/temp\/data\/(.*)\.gff/;
	my $name=$1;
	my %gff=();
	my %ffn=();
	my %faa=();
	my @line=();
	my %check=();
	my $line="";
	my @title=();
	
	open GFF , ">./temp/filter/$name.filter.gff";
	open IN , "$one";
	while(<IN>){
		$_=~tr/\r\n//d;
		@data=split(/\t/,$_);
		if($data[2]=~/gene/){
			$id=$_;
			$gff{$id}.="$_\n";
			push @title,$id;
			}
		if($data[2]=~/CDS/){
			$gff{$id}.="$_\n";
#			push @title,$id;
			}
		}
	close IN;
	
	foreach $two (@title){
		if($filter{$two}){
			$line.=$gff{$two};
			@data=split(/\t/,$two);
			$id="$data[3]\t$data[4]";
			$check{$id}=1;
			print GFF "$gff{$two}";
#			push @line,$gff{$two};
			}
		}
	close GFF;
  

  open FAA , "./temp/data/$name.faa" or die;
  
  while(<FAA>){
  	$_=~tr/\r\n//d;
  	if(/>/){
  		@data=split(/\|/,$_);
  		$id=$data[3];
  		$faa{$id}="$_\n";
  		next;
  		}
  	$faa{$id}.=$_;
  	}
  close FAA ;
  
  
  open OUT , ">./temp/filter/$name.filter.faa";
  foreach $two (keys %faa){
  	if($line=~/$two/){
  		print OUT "$faa{$two}\n";
  		}
  	}
  close OUT;
  
  open FFN , "./temp/data/$name.ffn";
  while(<FFN>){
  	$_=~tr/\r\n//d;
  	if(/>/){
  		@temp=split(/\s/,$_);
  		@data=split(/\|/,$temp[0]);
  		if($data[-1]=~/c(\d+)\D+(\d+)/){
  			$id=$2."\t".$1;
  			}
  		if($data[-1]!~/c/){
  			$data[-1]=~/(\d+)\D+(\d+)/;
  			$id=$1."\t".$2;
  			}
  		$ffn{$id}="$_\n";
  		next;
  		}
  	$ffn{$id}.=$_;
  	}
  close FFN;
  
  open OUT , ">./temp/filter/$name.filter.ffn";
  foreach $two (keys %ffn){
  	@data=split(/\t/,$two);
  	$id="$data[0]\t$data[1]";
  	if($check{$id}){
  		print OUT "$ffn{$two}\n";
  		}
  	}
  close OUT;
	}
