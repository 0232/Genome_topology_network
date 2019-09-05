#!usr/bin/perl
use strict;

chdir "Adjustment/cog";
system("mv cog.group.txt cog.group");
open IN, "cog.group" or die;
open OUT, ">cog.NoOut.txt";

if($ARGV[0]){
	while(<IN>){
		$_=~tr/\r\n//d;
		my @data=split(/\s/,$_);
		my $id=shift @data;
		print OUT "$id";

		foreach my $one (@data){
			$one=~/(.*)\|/;
			my $spe=$1;
			if(($spe eq $ARGV[0]) and $ARGV[0]){
				next;
			}
			print OUT " $one";
		}
		print OUT "\n";
	}
	close IN;
	close OUT;
}else{
	while(<IN>){
		print OUT $_;
	}
	close IN;
	close OUT;
}




