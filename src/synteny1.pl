#!usr/bin/perl
use strict;



#查找每个coord的one hit 坐标

my @file=glob("./temp/data/*.coord");
my $one;
my $two;
my $three;
my $four;
my $five;
my $id;
my @data;
my @temp;
my $name;
my $xx;
my $yy;
my %count;
my %one_hit;

open OUT , ">./temp/data/synteny_coordinate.txt";
foreach $one (@file){
	$one=~/\#(.*)\.coord/;
	$id=$1;
	my %sca=();
	
	open IN , "$one";
	while(<IN>){
		$_=~tr/\r\n//d;
	  if(/\// or /NUCMER/ or !($_) or /\[E1\]/ or /\=/){
		  next;
		  }
		@data=split(/\s\|\s/,$_);
		@temp=split(/\|/,$data[6]);
		$data[1]=~/\s(\d+)\s+(\d+)\s/;
		if($1<$2){
		  $sca{$temp[7]}.="$1..$2\t";
		  }
		if($1>$2){
		  $sca{$temp[7]}.="$2..$1\t";
		  }
		}
	close IN;
	
	foreach $two (keys %sca){
		
		@data=split(/\t/,$sca{$two});
		my %hash=();
		my @coor=();
		
		foreach $three (@data){
			$three=~/(\d+)\.\.(\d+)/;
			foreach $four ($1..$2){
				$hash{$four}++;
				}
			}
		$name=$id."#".$two;
		$count{$name}++;
		print OUT ">$name\n";
	  
	  foreach $five (keys %hash) {
	  	if($hash{$five}==1){
	  		push @coor , $five;
	  		}
	  	}
	  
	  @coor= sort {$a <=> $b} @coor;
	  $xx=$coor[0];
	  for(my $i=1;$i<=$#coor;$i++){
	  	if($coor[$i]-$coor[$i-1]==1){
	  		next;
	  		}
	  	if($coor[$i]-$coor[$i-1]>1){
	  		if($xx and $coor[$i-1]){
	  			print OUT "$xx..$coor[$i-1]\t";
	  			}	  		
        $one_hit{$name}.="$xx..$coor[$i-1]\t";
	  		$xx=$coor[$i];
	  		next;
	  		}
	  	}
	  if($xx and $xx..$coor[-1]){
	    print OUT "$xx..$coor[-1]\n";
	   }
		}
	}
	
#foreach $one (keys %one_hit){
#	@data=split(/\t/,$one_hit{$one});
#	my %hash=();
#	my @coor=();
#	
#	foreach $two (@data) {
#		$two=~/(\d+)\.\.(\d+)/;
#		
#		foreach $three ($1..$2){
#			$hash{$three}++;
#			}
#		}
#	foreach $four (keys %hash){
#		if($hash{$four}==$count{$one}){
#			push @coor , $four;
#			}
#		}
#  
#  print OUT ">$one\n";
#  @coor= sort {$a <=> $b} @coor;
#  
#  for(my $i=1;$i<=$#coor;$i++){
#	  	if($coor[$i]-$coor[$i-1]==1){
#	  		next;
#	  		}
#	  	if($coor[$i]-$coor[$i-1]>1){
#	  		print OUT "$xx..$coor[$i-1]\t";
#	  		$xx=$coor[$i];
#	  		next;
#	  		}
#	  	}
#	 print OUT "$xx..$coor[-1]\n";
#	}
close OUT;
