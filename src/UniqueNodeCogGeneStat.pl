#!usr/bin/perl
use strict;
use FileHandle;


my %unique_cog;
my %unique_gene;
my %miss_cog;
my %miss_gene;
my %delete_cog;
my %delete_gene;
my %insert_cog;
my %insert_gene;
my %clade;
my @title;
my %tmp;
open IN , "genes_in_unique_connection.txt" or die;

while(<IN>){
	$_=~tr/\r\n//d;

	if (/query_clade	reference_clade/) {
		next;
	}

	my @data=split(/\t/,$_);
	my $id="$data[0]\t$data[1]";

	if(!$tmp{$id}){
		push @title, $id;
		$tmp{$id}=1;
	}
	$clade{$id}=1;

	if($data[5]=~/(.*)\(delete/){
		$delete_cog{$id}.="$1\t";
		$delete_gene{$id}.="$1\|$data[7]\t";
	}elsif($data[8]=~/(.*)\(delete/){
		$delete_cog{$id}.="$1\t";
		$delete_gene{$id}.="$1\|$data[10]\t";
	}

	if($data[5]=~/(.*)\(insert/){
		$insert_cog{$id}.="$1\t";
		$insert_gene{$id}.="$1\|$data[7]\t";
	}elsif($data[8]=~/(.*)\(insert/){
		$insert_cog{$id}.="$1\t";
		$insert_gene{$id}.="$1\|$data[10]\t";
	}
	$data[5]=~s/\(.*\)//;
	$data[8]=~s/\(.*\)//;

	if($data[3] eq 'unique'){
		$unique_cog{$id}.="$data[5]\t$data[8]\t";
		$unique_gene{$id}.="$data[5]\|$data[7]\t$data[8]\|$data[10]\t";
	}elsif($data[3] eq 'miss'){
		$miss_cog{$id}.="$data[5]\t$data[8]\t";
		$miss_gene{$id}.="$data[5]\|$data[7]\t$data[8]\|$data[10]\t";
	}
}
close IN;
open UNICOG , ">unique_connection.cog.list";
print UNIGOG "query_clade\t reference_clade\tCOG\n";
open UNIGENE , ">unique_connection.gene.list";
print UNIGENE "query_clade\t reference_clade\tgene\n";
open MISCOG , ">miss_connection.cog.list";
print MISCOG "query_clade\t reference_clade\tCOG\n";
open MISGENE , ">miss_connection.gene.list";
print MISGENE "query_clade\t reference_clade\tgene\n";
open DELECOG , ">dele.cog.list";
print DELECOG "query_clade\t reference_clade\tCOG\n";
open DELEGENE , ">dele.gene.list";
print DELEGENE "query_clade\t reference_clade\tgene\n";
open INSECOG , ">insert.cog.list";
print INSECOG "query_clade\t reference_clade\tCOG\n";
open INSEGENE , ">inser.gene.list";
print INSEGENE "query_clade\t reference_clade\tgene\n";
open STAT , ">nodes_cog_gene.stat";
print STAT "clade\tclade\tunique_connection_number\tunique_gene_number\tmiss_cog_number\tmiss_gene_number\tdelete_cog_number\tdelete_gene_number\tinsert_cog_number\tinsert_gene_number\n";

foreach my $one (@title){
	my $id=$one;
	my $unique_cog_num=&hash($unique_cog{$id},$id,"UNICOG","COG");
	my $unique_gene_num=&hash($unique_gene{$id},$id,"UNIGENE","gene");
	my $miss_cog_num=&hash($miss_cog{$id},$id,"MISCOG","COG");
	my $miss_gene_num=&hash($miss_gene{$id},$id,"MISGENE","gene");
	my $dele_cog_num=&hash($delete_cog{$id},$id,"DELECOG","COG");
	my $dele_gene_num=&hash($delete_gene{$id},$id,"DELEGENE","gene");
	my $insert_cog_num=&hash($insert_cog{$id},$id,"INSECOG","COG");
	my $insert_gene_num=&hash($insert_gene{$id},$id,"INSEGENE","gene");
	print STAT "$one\t$unique_cog_num\t$unique_gene_num\t$miss_cog_num\t$miss_gene_num\t$dele_cog_num\t$dele_gene_num\t$insert_cog_num\t$insert_gene_num\n";
}
close UNICOG;
close UNIGENE;
close MISCOG;
close MISGENE;
close DELECOG;
close DELEGENE;
close INSECOG;
close INSEGENE;
close STAT;

sub hash{
	my %hash;
	my $strain=$_[0];
	chomp $strain;
	my $num=0;
	$_[2]->print("$_[1]\t");
	my @data=split(/\t/,$strain);
	foreach my $one (@data){
		$hash{$one}=1;
	}
	foreach my $one (keys %hash){
		$_[2]->print("$one ");
		$num++;
	}
	$_[2]->print("\n");
	return $num;
}