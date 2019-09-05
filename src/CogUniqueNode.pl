#!usr/bin/perl
use strict;

my %unique;
my %dup;
my %cds_annotation;
my %gene_start;
my %gene_end;
my %scaffold;

my @gff=glob("./Adjustment/gff/*.gff");

foreach my $one (@gff){
	open GFF , "$one";
	while(<GFF>){
		my @data=split(/\t/,$_);
		if($data[2] eq "CDS"){
			$_=~/product\=(.*?)\;/;
			my $anno=$1;
			$_=~/protein_id\=(.*?)\;/;
			my $cds=$1;
			$cds_annotation{$cds}=$anno;
			$gene_start{$cds}="$data[3]";
			$gene_end{$cds}="$data[4]";
			$scaffold{$cds}=$data[0];
			}
		}
	close GFF;
	}

my %cog_annotation;
open WHOG , "whog" or die;
while(<WHOG>){
	$_=~tr/\r\n//d;
	if(/\] (COG\d+) (.*)/){
		$cog_annotation{$1}=$2;
		}
	}
close WHOG;




my @title;
my %lie;
open TEST , ">genes_in_unique_connection.txt";
print TEST "query_clade\treference_clade\tconnection\tstatus\tgenome\tCOG1\tCOG1_function\tgene1\tCOG2\tCOG2_function\tgene2\n";
open VARIATION , "./Adjustment/tmp/ancenstor_unique_nodes.txt" or die;
while(<VARIATION>){
	$_=~tr/\r\n//d;
	my @data=split(/\t/,$_);
	my $id="$data[0]\t$data[1]";
	if(!$lie{$id}){
		push @title , $id;
		$lie{$id}=1;
		}
	if($data[3] eq "unique"){
		$data[2]=~/(.*)\-(.*)/;
		my $node;
		my $id1=$1;
		my $id2=$2;
		$id1=~/\D+(\d+)/;
		my $num1=$1;
		$id2=~/\D+(\d+)/;
		my $num2=$1;
		if($num1<=$num2){
			$node="$id1\-$id2";
		}else{
			$node="$id2\-$id1";
			}
		$unique{$id}.="$node\t";
		}
	}
close VARIATION;



my %unique_result;
my %miss_result;

foreach my $one (keys %unique){
	$one=~/(.*)\t(.*)/;
	my $other_id="$2\t$1";
	my @genome=split(/\+/,$1);
	my @other_genome=split(/\+/,$2);
	my @strain_node=split(/\t/,$unique{$one});
	my %cog_node=node(@genome);
	my %other_cog_num=cog_num(@other_genome);
  
#  foreach my $ii (keys %other_cog_num){
#  	print "$ii	$other_cog_num{$ii}\n";
#  	}
  
  foreach my $two (@strain_node){
  	if(!$cog_node{$two}){
  		next;
  		}

  	my @unique_cog=split(/\n/,$cog_node{$two});
  	
  	foreach my $four (@unique_cog){
  		my @temp=split(/\t/,$four);                  # $temp[0]\t$temp[1]\t$temp[2]\t$temp[3]
   		                                             #   cog        gene      cog      gene

  		$temp[1]=~/\|(.*)__(.*)/;
  		my $spe=$2;
  		my $gene1=$1;
  		$temp[3]=~/\|(.*)__(.*)/;
  		my $gene2=$1;
  		my $cog1=$temp[0];
  		my $cog2=$temp[2];
  		my $dele1=$cog1;
  		my $dele2=$cog2;
  		if(!$other_cog_num{$temp[0]}){
  			$cog1="$temp[0](insert)";
  			$dele1="$temp[0](delete)";
  			}
  		if(!$other_cog_num{$temp[2]}){
  			$cog2="$temp[2](insert)";
  			$dele2="$temp[2](delete)";
  			}
  		$unique_result{$one}.="$one\t$two\tunique\t$spe\t$cog1\t$cog_annotation{$temp[0]}\t$gene1\t$cog2\t$cog_annotation{$temp[2]}\t$gene2\n";
      $miss_result{$other_id}.="$other_id\t$two\tmiss\t$spe\t$dele1\t$cog_annotation{$temp[0]}\t$gene1\t$dele2\t$cog_annotation{$temp[2]}\t$gene2\n"
  		}
  	}
	}
foreach my $one (@title){
	print TEST "$unique_result{$one}$miss_result{$one}";
	}
close TEST;


sub node{
	my %hash;
	my $id;
	my $node;
	
	foreach my $one (@_){
  my @data=();
	open IN , "./Adjustment/topo/$one.topo" or die;
	while(<IN>){
		$_=~tr/\r\n//d;
		push @data , $_;
		}
	close IN;
	
	for(my $i=0;$i<$#data;$i++){
#		print TEST "$data[$i]\n";
		$data[$i]=~/(\D+(\d+))\t/;	
		my $id1=$1;
		my $num1=$2;
		$data[$i+1]=~/(\D+(\d+))\t/;
		my $id2=$1;
		my $num2=$2;

		if($num1<=$num2){
			$id="$id1\-$id2";
			$node="$data[$i]\t$data[$i+1]";
		}else{
			$id="$id2\-$id1";
			$node="$data[$i+1]\t$data[$i]";
			}
		$hash{$id}.="$node\n";
		}
	
		$data[-2]=~/(\D+(\d+))\t/;	
		my $id1=$1;
		my $num1=$2;
		$data[-1]=~/(\D+(\d+))\t/;
		my $id2=$1;
		my $num2=$2;
		if($num1<=$num2){
			$id="$id1\-$id2";
			$node="$data[-2]\t$data[-1]";
		}else{
			$id="$id2\-$id1";
			$node="$data[-1]\t$data[-2]";
			}

		$hash{$id}.="$node\n";
	}
		return %hash;
	}

sub cog_num{
	my %hash;
	
	foreach my $one (@_){
		open IN , "./Adjustment/topo/$one.cog.num" or die;
		while(<IN>){
			$_=~tr/\r\n//d;
			my @data=split(/\t/,$_);
			$hash{$data[0]}+=$data[1];
			}
		close IN;
		}
	return %hash;
	}