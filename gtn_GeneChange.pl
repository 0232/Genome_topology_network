#!usr/bin/perl
use strict;

# Usage: perl gtn_GeneChange.pl <nwk file> <1 (complete only)|2 (draft)> <outgroup>
# If there is no outgroup , then just: perl gtn_GeneChange.pl <nwk file> <1 complete only|2 draft> <outgroup>
# For exsample: perl gtn_GeneChange.pl test.nwk 1
# Xiao Deng



system("perl src/RelativeDdStat.pl Adjustment/tmp/mobile.list");

open IN, "$ARGV[0]" or die;    # bootstrap file
my $outgroup=$ARGV[2];    # outgroup
open OUT, ">cog.nwk";

while (<IN>){
	$_ =~ s/(\))(\d)/$1\:$2/g;
	print OUT $_;
}
close IN;

system("mv ./Adjustment/cog/cog.group ./Adjustment/cog/cog.group.txt");
system("mv ./Adjustment/cog/cog.NoOut.txt ./Adjustment/cog/cog.NoOut");
system("perl src/ChangeTb.pl ./Adjustment/cog ./Adjustment/gff gff > Adjustment/tmp/gene.change.txt");
system("perl src/CogGeneTopo.pl");
system("perl src/CogNum.pl");
system("perl src/CommonAncestor.pl Adjustment/tmp/gene.change.txt $outgroup");
system("perl src/CogUniqueNode.pl");
system("perl src/UniqueNodeCogGeneStat.pl");

if($ARGV[1]==1){
	system("perl src/pie_OrhtoCog.topo.pl");
	system("perl src/pie_Connection.pl");
	system("perl src/pie_stat.pl");
}

if($ARGV[1]==2){
	system("perl src/pie_OrhtoCog.topo.pl");
	system("perl src/pie_Connection.pl");
	system("perl src/dra_pie_seq_filter.pl");
	system("perl src/dra_pie_stat.pl");
}
