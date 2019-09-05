#!usr/bin/perl
use strict;


my %para;
my %dd;
my @data;
my $id;
my $n;
my %str;
my %dd;
my %gene;
open PARA , "Adjustment/cog/cog.NoOut.txt" or die;

while(<PARA>){
	$_=~/(COG\d+)\:/;
	$id=$1;
	$para{$id}=($_=~s/\|/\|/g);
	@data=split(/\s/,$_);
	my %hash;
	my $num;
	my %tmp;
	my $gene;
	foreach my $one (@data){
		if($one=~/(.*)\|(.*)__/){
			$hash{$1}=1;
			$tmp{$2}=1;
			}
		}
	
	foreach my $one (keys %hash){
		$num++;
		}
	foreach my $one (keys %tmp){
		$gene++;
		}
	
	$str{$id}=$num;
	$gene{$id}=$gene;
	}
close PARA;

my $relative;

open DD , "$ARGV[0]" or die;
open OUT , ">relative_dd.txt";
print OUT "COG\tCOG_function\tDD\/str\tPara\/str\tRelative_dd\n";
while(<DD>){
	$_=~tr/\r\n//d;
	if(/COG ID	Annotation	Different Degree/){
		next;
		}
	@data=split(/\t/,$_);
	if($para{$data[0]}==0){
		next;
		}
	$relative=$data[2]/$para{$data[0]}*$str{$data[0]};
	my $para_str=$para{$data[0]}/$str{$data[0]};
	print OUT "$data[0]\t$data[1]\t$data[2]\t$para_str\t$relative\n";
	}
close DD;
close OUT;