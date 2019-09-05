#!usr/bin/perl
use strict;

chdir "Adjustment/tmp";
system("formatdb -i all.pep -p T");
system("blastall -p blastp -i all.pep -d all.pep -o all.pep.blast -m 9 -e 1e-5 -a $ARGV[0]");
system("mcxdeblast --m9 --line-mode=abc --score=r all.pep.blast \| mcl - --abc -o all.pep.clu");
my @clu;
my %tab;
my %pep_id;
my @data;
my $id;
open IN , "all.pep.clu" or die;               #mcl clu file
while(<IN>){
	$_=~tr/\r\n//d;
	push @clu , $_;
	}
close IN;

open PEP, "all.pep" or die;
while(<PEP>){
	$_=~tr/\r\n//d;
	if(/>(.*)/){
		$pep_id{$1}=1;
	}
}
close PEP;

open TAB, "../../src/cog.id.tab" or die;                      #cog family table
while(<TAB>){
	$_=~tr/\r\n//d;
	@data=split(/\t/,$_);
	$tab{$data[1]}=$data[0];
	}
close TAB;

open OUT , ">../cog/cog.group.txt";
my %result;
my $only_cog;
foreach my $one (@clu){
	my @temp=();
	my %hash=();
	my $line="";
	my $n=0;
	my $cog_count=0;
	my @gene=();
	my @all=();
	$id="";
	my $gene="";
	
	
	@data=split(/\t/,$one);
	
	foreach my $two (@data){
		if(!($tab{$two}=~/COG/) and !$pep_id{$two}){
			next;
			}
		if($tab{$two}){
			$hash{$tab{$two}}=1;
			next;
			}
		$hash{$two}=1;
		}
	
	foreach my $three (keys %hash){
		$line.="$three\t";
		}
	

	$cog_count=($line=~s/COG/COG/g);
	
	
	if($cog_count>=1){
		@temp=split(/\t/,$line);
		my $n=@temp;
		if($cog_count==$n){
			next;
			}
		if($cog_count==1){
			foreach my $five (@temp){
				if($five=~/COG/){
					$id=$five;
					next;
					}
				if($five=~/\.\d\__(.*)/){
					$five=~/\.\d\__(.*)/;
					$five=$1."|".$five;
					$gene.=" $five";				
				}			
			}
			if($gene){
				$result{$id}.="$gene";
			}		
		}
	}
}
foreach (keys %result){
	print OUT "$_:$result{$_}\n";
	}
close OUT;