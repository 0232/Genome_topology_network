#!usr/bin/perl
use strict;

my %ortho;
my %cog;


mkdir "Adjustment/pie" unless -d "Adjustment/pie";
mkdir "Adjustment/pie/topo" unless -d "Adjustment/pie/topo";

=pod
open ORTHO, "pie/ortho.group.txt";
while(<ORTHO>){
	$_=~tr/\r\n//d;
	my @data=split(/\s/,$_);
	my $id=shift @data;
	$id=~s/\://;

	foreach my $one (@data){
		$one=~s/(.*)\|//;
		$ortho{$one}=$id;
	}
}
close ORTHO;
=cut

open COG, "Adjustment/cog/cog.group.txt" or die;

while(<COG>){
	$_=~tr/\r\n//d;
	my @data=split(/\s/,$_);
	my $id=shift @data;
	$id=~s/\://;

	foreach my $one (@data){
		$one=~s/(.*)\|//;
		$cog{$one}=$id;
	}
}
close COG;
my @gff_file=glob("Adjustment/gff/*.gff");
my $name_order=1;

foreach my $one (@gff_file){
	$one=~/gff\/(.*)\.gff/;
	my $genome=$1;
	open OUT, ">Adjustment/pie/topo/$genome.pie.topo";
	open IN, "$one";
	my $ord=1;

	while(<IN>){
		my @data=split(/\t/,$_);

		if(/protein_id\=(.*?)\;/){
			my $gene=$1;
			if($cog{$gene}){
				if(!$ortho{$gene}){
					$ortho{$gene}='*';
				}
				print OUT "$cog{$gene}\t$data[0]\t$ord\t$gene\t$data[3]\t$data[4]\n";
				$ord++;
			}
		}
	}
	close IN;
	close OUT;
}