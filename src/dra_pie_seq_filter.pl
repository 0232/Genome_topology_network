#!usr/bin/perl
use strict;

my %syn;
my $id;
open IN, "temp/data/synteny_once_coordinate.txt";
while(<IN>){
	$_=~tr/\r\n//d;
	if(/>.*\#(.*)/){
		$id=$1;
		next;
		}
	$syn{$id}.=$_;
	}
close IN;

my @file=glob("Adjustment/pie/*.seq.pie");
foreach my $c1 (@file){
	if($c1=~/seq.filter.pie/){
		next;
	}
	
	$c1=~/.*\/(.*).seq.pie/;
	my $name=$1;
	open OUT, ">Adjustment/pie/$name.seq.filter.pie";
	open IN, "$c1";
	
	while(<IN>){
		$_=~tr/\r\n//d;;
		my @data=split(/\t/,$_);
		my $genome=shift @data;
		print OUT "$genome";
		
		my %pie_coor;
		foreach my $c2 (@data){
			$c2=~/(.*)\|(\d+)\.\.(\d+)/;
			my $sca=$1;
			my $xx=$2;
			my $yy=$3;
			if(!$syn{$sca}){
				next;
				}
			$pie_coor{$sca}.="$xx\.\.$yy\t";
			}
		
		foreach my $c2 (keys %pie_coor){
			my %hash_syn;
			my @data_syn=split(/\t/,$syn{$c2});
			
			foreach my $c3 (@data_syn){
				$c3=~/(\d+)\.\.(\d+)/;
				foreach my $c4 ($1..$2){
					$hash_syn{$c4}=1;
					}
				}		
			my @data_pie=split(/\t/,$pie_coor{$c2});
			
			foreach my $c3 (@data_pie){
				$c3=~/(\d+)\.\.(\d+)/;
				my @res_coor;
				my $aa=$1;
				my $bb=$2;
				
				foreach my $c4 ($aa..$bb){
					if($hash_syn{$c4}){
						push @res_coor, $c4;
						}
					}
				my @res=&coor(@res_coor);
				
				foreach my $c4 (@res){
					print OUT "\t$c4";
					}
				}
			}
		print OUT "\n";
		}
	close OUT;
	}
	
sub coor{
	my @line;
	my @data=@_;
	my $start=$data[0];
	

	for(my $i=1;$i<$#data;$i++){
		if($data[$i]-$data[$i-1]>=2){
			my $coor="$start\.\.$data[$i-1]";
#			print "$coor\n";
			push @line, $coor;
			$start=$data[$i];
			next;
			}
		}
	push @line, "$start\.\.$data[-1]";
	return @line;
	}