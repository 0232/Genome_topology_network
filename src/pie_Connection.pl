#!usr/bin/perl
use strict;


my %ance;

open COM, "Adjustment/tmp/ComAnc.nodes.txt" or die;
while(<COM>){
	$_=~tr/\r\n//d;
	my @data=split(/\t/,$_);
	my $id=shift @data;
	$id=~s/\://;
	my $line=join("\t",@data);
	$ance{$id}=$line;
}
close COM;


my $file_order=1;
open TAB, ">Adjustment/pie/file2branch.tab";

foreach my $one (keys %ance){
	my %node;
	my @node=split(/\t/,$ance{$one});
	
	foreach my $two (@node){
		$node{$two}=1;
		}
	my @genome=split(/\+/,$one);
	open GENE_PIE, ">Adjustment/pie/$file_order.gene.pie";
	open SEQ_PIE, ">Adjustment/pie/$file_order.seq.pie";
	print TAB "$file_order\t$one\n";
	$file_order++;
	foreach my $two (@genome){
		my @topo;
		my %gene_coor;
		
		open TOPO, "Adjustment/pie/topo/$two.pie.topo" or die;
		print GENE_PIE "$two";
		print SEQ_PIE "$two";
		
		while(<TOPO>){
			$_=~tr/\r\n//d;
			push @topo, $_;
			my @data_tmp=split(/\t/,$_);
			my $id="$two\t$data_tmp[2]";
			$gene_coor{$id}="$data_tmp[-2]\.\.$data_tmp[-1]";
			}
		close TOPO;
		my %sca;
		
		for(my $i=1;$i<=$#topo;$i++){
			my @data1=split(/\t/,$topo[$i-1]);
			my @data2=split(/\t/,$topo[$i]);
			my $node1=$data1[0].'-'.$data2[0];
			my $node2=$data2[0].'-'.$data1[0];
			
			if($node{$node1} and !($node1 eq $node2) and $data1[1] eq $data2[1]){
				$sca{$data1[1]}.="$data1[2]\t$data2[2]\t";
			}elsif($node{$node2} and !($node1 eq $node2) and $data1[1] eq $data2[1]){
				$sca{$data1[1]}.="$data1[2]\t$data2[2]\t";
			}elsif($node{$node2} and $node1 eq $node2 and $data1[1] eq $data2[1]){
				$sca{$data1[1]}.="$data1[2]\t$data2[2]\t";
				}
			}
		
		foreach my $c3 (keys %sca){
				my @order=split(/\t/,$sca{$c3});
				my %hash;
				
				foreach my $c4 (@order){
					$hash{$c4}=1;
					}
				my @sort_ord;
				
				foreach my $c4 (keys %hash){
					push @sort_ord, $c4;
					}
				@sort_ord=sort {$a<=>$b} @sort_ord;
				my @line=&coor(@sort_ord);
				
				foreach my $c4 (@line){
					print GENE_PIE "\t$c4";
					$c4=~/(\d+)\.\.(\d+)/;
					my $xx=$1;
					my $yy=$2;
					my $id1="$two\t$xx";
					my $id2="$two\t$yy";
					$gene_coor{$id1}=~/(\d+)\.\./;
#					print "$gene_coor{$id1}\n";
					my $coor1=$1;
					$gene_coor{$id2}=~/\.\.(\d+)/;
					my $coor2=$1;
					print SEQ_PIE "\t$c3\|$coor1\.\.$coor2";
					}
				}	
		print GENE_PIE "\n";
		print SEQ_PIE "\n";
		}
		close GENE_PIE;
		close SEQ_PIE;
	}
close TAB;
sub coor{
	my @line;
	my @data=@_;
	my $start=$data[0];
	

	for(my $i=1;$i<$#data;$i++){
		if($data[$i]-$data[$i-1]>=2){
			my $coor="$start\.\.$data[$i-1]";
			push @line, $coor;
			$start=$data[$i];
			next;
			}
		}
	push @line, "$start\.\.$data[-1]";
	return @line;
	}