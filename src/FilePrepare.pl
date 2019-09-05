#!usr/bin/perl
use strict;




#copy complete genome data to temp/data
# creat ffn file

my $one;
my $two;
my $name;
my @com_file=glob("./data_complete/*");
my $xiabiao="__";
foreach $one (@com_file){
	$one=~/complete\/(.*)/;
	my $id=$1;

  my @fna=glob("$one/*.fna");
  my @gff=glob("$one/*.gff");
  my @faa=glob("$one/*.faa");
  
	open GFF , ">./temp/data/$id.gff";
	open FNA , ">./temp/filter/$id.fna";
	open FAA , ">./temp/data/$id.faa";
	
	
  foreach $two (@gff) {
		  open IN , "$two";
		  while(<IN>){
		  	if(/protein_id\=/){
		  		$_=~s/protein_id\=(.*?)\;/protein_id\=$1$xiabiao$id\;/;
		  		}
			  print GFF "$_";
		  	}
		  close IN;
		  }
	
	foreach $two (@fna){
		open IN , "$two";
		while(<IN>){
                        $_=~tr/\r\n//d;
			if(/>(.*?)\s(.*)/){
				print FNA ">gi\|unkown\|ref\|$1\| $2\n";
				next;
				}
			print FNA "$_\n";
			}
		close IN;
		}
	
	foreach $two (@faa){
		open IN , "$two";
		while(<IN>){
                         $_=~tr/\r\n//d;
			if(/>(.*?)\s(.*)/){
				print FAA ">gi\|unkown\|ref\|$1$xiabiao$id\| $2\n";
				next;
				}
			print FAA "$_\n";
			}
		close IN;
		}
	close GFF;
	close FNA;
	close FAA;

	}



if($ARGV[0]==2){
	#move draft genome data to temp/data 
my @dra_file=glob("./data_draft/*");

my $name;

foreach $one (@dra_file) {
	$one=~/data_draft\/(.*)/;
	$name=$1;
#	my @tgz=glob("$one/*.gz");
#	
#	foreach $two (@tgz) {
#		system ("gunzip -f $two $one");         
#		}
	
	my @gff=glob("$one/*.gff");
	my @fna=glob("$one/*.fna");
	my @faa=glob("$one/*.faa");
	
	open GFF , ">temp/data/$name.gff";
	open FNA , ">temp/filter/$name.fna";
	open FAA , ">temp/data/$name.faa";
	foreach $two (@gff) {
		open IN , "$two";
		while(<IN>){
			if(/protein_id\=/){
		  		$_=~s/protein_id\=(.*?)\;/protein_id\=$1$xiabiao$name\;/;
		  		}
			print GFF "$_";
			}
		close IN;
		}
	
	foreach $two (@fna){
		open IN , "$two";
		while(<IN>){
      $_=~tr/\r\n//d;
			if(/>(.*?)\s(.*)/){
				print FNA ">gi\|unkown\|ref\|$1\| $2\n";
				next;
				}
			print FNA "$_\n";
			}
		close IN;
		}

	foreach $two (@faa){
		open IN , "$two";
		while(<IN>){
      $_=~tr/\r\n//d;
			if(/>(.*?)\s(.*)/){
				print FAA ">gi\|unkown\|ref\|$1$xiabiao$name\| $2\n";
				next;
				}
			print FAA "$_\n";
			}
		close IN;
		}
	close GFF;
	close FNA;
	close FAA;
#        system ("rm $one/*.gff");
#        system ("rm $one/*.fna");
##        system ("rm $one/*.ffn");
#        system ("rm $one/*.faa");
	}
}




#creat ffn


my $id;
my %seq;
my @data;
my @seq=glob("./temp/filter/*.fna");

foreach $one (@seq){
	open IN , "$one";
	while(<IN>){
		$_=~tr/\r\n//d;
		if(/>/){
                        my @data=split(/\|/,$_);
			$id=$data[3];
			next;
			}
		$seq{$id}.=$_;
		}
	close IN;
	}

my @gff=glob("./temp/data/*.gff");
foreach $one (@gff){
	$one=~/temp\/data\/(.*)\.gff/;
	$name=$1;
	open OUT , ">./temp/data/$name.ffn";
	open IN , "$one";
	
	while(<IN>){
		if($_=~/pseudo\=true/){
			next;
			}
		@data=split(/\t/,$_);
		if($data[2] eq "CDS"){
			if($data[6] eq "+"){
			  print OUT ">gi\|unkown\|ref\|$data[0]\|\:$data[3]\-$data[4] $name\n";
			  my $offset=$data[3]-1;
			  my $len=$data[4]-$data[3]+1;
			  my $seq=substr($seq{$data[0]},$offset,$len);
			  print OUT "$seq\n";
		    }
		   if($data[6] eq "-"){
		   	print OUT ">gi\|unkown\|ref\|$data[0]\|\:c$data[4]\-$data[3] $name\n";
		   	my $offset=$data[3]-1;
			  my $len=$data[4]-$data[3]+1;
			  my $seq=substr($seq{$data[0]},$offset,$len);
			  $seq=~tr/ATCG/TAGC/;
			  $seq=reverse $seq;
			  print OUT "$seq\n";
		   	}
			}
		}
	close IN;
	close OUT;
	}
