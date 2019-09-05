#!usr/bin/perl
use strict;


my %branch;
my @pair;
my $outgroup="$ARGV[1]";
my %sca;
my @gff_file=glob("./Adjustment/gff/*.gff");

foreach my $one (@gff_file){
	open IN , "$one";
	while(<IN>){
		my @data=split(/\t/,$_);
		if($data[2] eq 'CDS'){
			/protein_id\=(.*?)\;/;
			$sca{$1}=$data[0];
		}
	}
	close IN;
}

open CHA, "$ARGV[0]" or die;   # gene change table
while(<CHA>){
	$_=~tr/\r\n//d;
	my @data=split(/\t/,$_);
	push @pair, "$data[0]\t$data[1]";
	$branch{$data[0]}=1;
	$branch{$data[1]}=1;
}
close CHA;

open NODE, ">./Adjustment/tmp/ComAnc.nodes.txt";
my %branch_common_nodes;
foreach my $one (keys %branch){
	my %nodes_in_branch;
	my $node_strain;
	if($one eq $outgroup){
		next;
	}

	if($outgroup){
		if($one=~/$outgroup\+/){
			$one=~s/$outgroup\+//;
		}elsif($one=~/\+$outgroup/){
			$one=~s/\+$outgroup//;
		}
	}
	print NODE "$one\:";
	my @genome=split(/\+/,$one);
	foreach my $two (@genome){
		my %hash=&node($two,\%sca);
		foreach my $three (keys %hash){
			$nodes_in_branch{$three}++
		}
	}
	my $genome_num_branch=@genome;
	my $common_nodes_num;
	foreach my $two (keys %nodes_in_branch){
		if($nodes_in_branch{$two}==$genome_num_branch){
			$common_nodes_num++;
			$node_strain.="\t$two";
		}
	}
	print NODE "$node_strain\n";
	$branch_common_nodes{$one}=$node_strain;
}
close NODE;

open ANC, ">./Adjustment/tmp/ancenstor_unique_nodes.txt";
open ANCNUM , ">./Adjustment/tmp/ancenstor_nodes.stat";
print ANCNUM "query_clade\tsubject_clade\tcommon_nodes_number\tq_clade_unique_nodes\tq_clade_miss_node\n";
foreach my $one (@pair){
	$one=~/(.*)\t(.*)/;
	my $pair1=$1;
	my $pair2=$2;
	my %com_nodes;
	my %node1;
	my %node2;
	my @unique_node1;
	my $unique_node1_num;
	my @unique_node2;
	my $unique_node2_num;
	my $com_nodes_num;
	if(($pair1 eq $outgroup) or ($pair2 eq $outgroup)){
		next;
	}

	if($outgroup){
		if($pair1=~/$outgroup\+/){
			$pair1=~s/$outgroup\+//;
		}elsif($pair1=~/\+$outgroup/){
			$pair1=~s/\+$outgroup//;
		}
		if($pair2=~/$outgroup\+/){
			$pair2=~s/$outgroup\+//;
		}elsif($pair2=~/\+$outgroup/){
			$pair2=~s/\+$outgroup//;
		}
	}	
	my @node1=split(/\t/,$branch_common_nodes{$pair1});
	my @node2=split(/\t/,$branch_common_nodes{$pair2});
	foreach my $two (@node1){
		$node1{$two}=1;
	}
	foreach my $two (@node2){
		$node2{$two}=1;
	}
	foreach my $two (keys %node1){								# find common nodes
		if($node2{$two}){
			$com_nodes{$two}=1;
			$com_nodes_num++;
		}else{													# find unique nodes
			print ANC "$pair1\t$pair2\t$two\tunique\n";
			$unique_node1_num++;
		}
	}
	foreach my $two (keys %node2){
		if(!$node1{$two}){
			print ANC "$pair1\t$pair2\t$two\tmiss\n";						#find miss nodes
			$unique_node2_num++;
		}
	}
	print ANCNUM "$pair1\t$pair2\t$com_nodes_num\t$unique_node1_num\t$unique_node2_num\n";
}
close ANC;
close ANCNUM;


sub node{
	my %hash;
	my @cog;
	my $aa=$_[1];
	my %sca=%$aa;
	open IN , "./Adjustment/topo/$_[0].topo" or die;

	while(<IN>){
		$_=~tr/\r\n//d;
#		my @data=split(/\t/,$_);
		push @cog, $_;
	} 
	close IN;
	my $node;
	for(my $i=0;$i<$#cog; $i++){
		$cog[$i]=~/(.*\D(\d+))\t.*\|(.*)/;
		my $cog1=$1;
		my $xx=$2;
		my $gene1=$3;
		$cog[$i+1]=~/(.*\D(\d+))\t.*\|(.*)/;
		my $cog2=$1;
		my $yy=$2;
		my $gene2=$3;
		if(!($sca{$gene1} eq $sca{$gene2})){
			next;
		}
		if($xx<=$yy){
			$node="$cog1\-$cog2";
		}elsif($xx>$yy){
			$node="$cog2\-$cog1";
		}
		$hash{$node}=1;
	}
	$cog[-2]=~/(.*\D(\d+))\t.*\|(.*)/;
	my $cog1=$1;
	my $xx=$2;
	my $gene1=$3;
	$cog[-1]=~/(.*\D(\d+))\t.*\|(.*)/;
	my $cog2=$1;
	my $yy=$2;
	my $gene2=$3;
	if($xx<=$yy){
		$node="$cog1\-$cog2";
	}elsif($xx>$yy){
		$node="$cog2\-$cog1";
	}
	if($sca{$gene1} eq $sca{$gene2}){
		$hash{$node}=1;
	}
	return %hash;
}
