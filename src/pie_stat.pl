#!usr/bin/perl
use strict;


my %fra_tab;
open TAB, "Adjustment/pie/file2branch.tab" or die;
while(<TAB>){
	$_=~tr/\r\n//d;
	/(.*)\t(.*)/;
	$fra_tab{$1}=$2;
}
close TAB;

open INFO , ">fragment_connection.info";
print INFO "clade\taverage_number\taverage_length \(KB\)\n";
my @gene_pie_file=glob("Adjustment/pie/*.gene.pie");

foreach my $one (@gene_pie_file){
	$one=~/pie\/(.*)\.gene.pie/;
	my $file_name=$1;
	open OUT , ">Adjustment/pie/$fra_tab{$file_name}.gene.stat";
	my $gene_num=0;
	my $pie_num=0;
	my $genome_num=($fra_tab{$file_name}=~s/\+/\+/g)+1;

	open IN , "$one";
	while(<IN>){
		$_=~tr/\r\n//d;
		my @data=split(/\t/,$_);
		my $id=shift @data;
		print OUT "$id\t";
		$pie_num+=@data;

		foreach my $two (@data){
			$two=~/(\d+)\.\.(\d+)/;
			$gene_num+=$2-$1+1;
		}
		print OUT "$gene_num\n";
	}
	close IN;
	close OUT;
	my $aver_gene=$gene_num/$genome_num;
	$aver_gene=~s/\.\d+//;
	my $aver_pie=$pie_num/$genome_num;
	$aver_pie=~s/\.\d+//;
	if($fra_tab{$file_name}!~/\+/){
		next;
	}
	print INFO "$fra_tab{$file_name}\t$aver_pie";
	my $seq_long=0;

	open SEQ, "Adjustment/pie/$file_name.seq.pie" or die;
	while(<SEQ>){
		$_=~tr/\r\n//d;
		my @data=split(/\t/,$_);
		my $id=shift @data;

		foreach my $two (@data){
			$two=~/\|(\d+)\.\.(\d+)/;
			$seq_long+=$2-$1+1;
		}
	}
	close SEQ;
	my $aver_long=$seq_long/$genome_num;
	$aver_long=$aver_long/1000;
	$aver_long=~s/(\.\d\d)\d/$1/;
	print INFO "\t$aver_long\n";

}
close INFO;
