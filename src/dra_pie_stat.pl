#!usr/bin/perl
use strict;


my @file=glob("Adjustment/pie/*.seq.filter.pie");
open OUT, ">fragment_connection_with_draft.info";
print OUT "branch\taverage_length \(KB\)\taverage_number\n";
my %tab;
open TAB, "Adjustment/pie/file2branch.tab" or die;
while(<TAB>){
	$_=~tr/\r\n//d;
	my @data=split(/\t/,$_);
	$tab{$data[0]}=$data[1];
	}
close TAB;

foreach my $c1 (@file) {
	my $geno_num;
	my $len;
	my $pie;
	$c1=~/.*\/(.*).seq.filter.pie/;
	my $name=$tab{$1};
	if($name!~/\+/){
		next;
	}
	open IN, "$c1";
	while(<IN>){
		$_=~tr/\r\n//d;
		$geno_num++;
		my @data=split(/\t/,$_);
		my $tmp=shift @data;
		$pie+=scalar @data;

		foreach my $c2 (@data){
#			$pie++;
			$c2=~/(\d+)\.\.(\d+)/;
			my $xx=$1;
			my $yy=$2;
			$len+=$yy-$xx+1;
			}
		}
	close IN;
	
	my $aver_len=$len/$geno_num;
	my $aver_pie=$pie/$geno_num;
	$aver_len=$aver_len/1000;
	$aver_len=~s/\.\d+//;
	print OUT "$name\t$aver_len\t$aver_pie\n";
	}
close OUT;
