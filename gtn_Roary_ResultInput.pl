#!usr/bin/perl
use strict;


# This script takes cluster result from Roary tool as gene family assignment. This method can markedly reduce running time.
# Input data: file 'clustered_proteins' and correct gff files which are used as input data in Roary are required. Put all these files into one folder. If your roary result contains gff files in "fixed_input_files", please use them rather than the original gff files.
# Usage: perl gtn_Roary_ResultInput.pl <roary_result_folder_path> <1:complete data | 2: draft data> <outgroup>
# For example: perl GTN_Roary_ResultInput.pl /home/dx/test/roary_result 1 S.out
# An error message might appear : "Use of uninitialized value in print at src/GetMobileGOGs.pl line 23, <COG> line 101336.", ignore it.
# Xiao Deng



my $roary_result_path=$ARGV[0];

if ( ! -d "roary_temp" ) {
	mkdir "roary_temp";
} else {
	system ("rm -r roary_temp/*");
}

my %species;
my @gff_file = glob ("$roary_result_path/*.gff");


foreach my $c1 (@gff_file){

	my $genome_name = $c1;
	$genome_name =~ s/.*\///;
	$genome_name =~ s/\.gff//;

	if ( ! -d "roary_temp/$genome_name" ) {
		mkdir "roary_temp/$genome_name";
	} else {
		system ("rm roary_temp/$genome_name/*");
	}
	
	open GFF_OUT, ">roary_temp/$genome_name/$genome_name.gff";
	open FNA_OUT, ">roary_temp/$genome_name/$genome_name.fna";
	open FAA_OUT, ">roary_temp/$genome_name/$genome_name.faa";
	
	open GFF, "$c1";
	my $swi=1;
	my $seq_id;
	my %seq_hash;

	while(<GFF>){

		$_=~tr/\r\n//d;

		if ( /\#\#FASTA/ ){
			$swi = 2;
			next;
		}

		if ( /\#/ ) {
			next;
		}

		my @data=split(/\t/,$_);

		if ( $data[2]=~/CDS/ and $swi==1 ){
			/\tID\=(.*?)\;/;
			my $pro_id = $1;
			$species{$pro_id} = $genome_name;
			$_ =~ s/protein_id\=(.*?)\;/protein_id\=${pro_id}\;/;
			print GFF_OUT $_."\n";
			next;
		}

		if ( $swi == 1 ) {
			print GFF_OUT $_."\n";
		}
		
		if ( $swi==2  and  />(.*)/ ) {				# genome sequence id
			$seq_id = "$1 \[roary result\]";
			next;
		}

		if ( $swi==2 and $_ !~/>/  and  $_ ) {
			$seq_hash{$seq_id} .= $_;
		}

	}

	close GFF;
	close GFF_OUT;

	foreach my $c2 (keys %seq_hash) {
		my $seq_line = $seq_hash{$c2};
		$seq_line =~ tr/ATCGN//d;
		my $seq_length = length $seq_line;

		if ($seq_length > 10) {
			print FAA_OUT ">$c2\n$seq_hash{$c2}\n";
		} else {
			print FNA_OUT ">$c2\n$seq_hash{$c2}\n";
		}

	}

	close FNA_OUT;
	close FAA_OUT;
}


open CLUSTER, "$roary_result_path/clustered_proteins" or die "No clustered_proteins file !\n";
open CLU_OUT, ">roary_temp/roary.group.txt";
my $clu_num=1;

while(<CLUSTER>){

	$_ =~ tr/\r\n//d;
	my @data = split (/\s/, $_);
	my $tmp = shift @data;

	print CLU_OUT "cluster${clu_num}\:";

	foreach my $c1 (@data) {
		print CLU_OUT " $species{$c1}\|${c1}__$species{$c1}";
	}

	print CLU_OUT "\n";
	$clu_num++;

}

close CLUSTER;
close CLU_OUT;




# make distance file

system ("mkdir data_complete");
mkdir "temp" unless -d "temp";

if ( -d "temp/data" ){
	system ("rm -rf temp/data/*");
} else {
	mkdir "temp/data";
}

if ( -d "temp/filter" ){
	system ("rm -rf temp/filter/*");
} else {
	mkdir "temp/filter";
}

if ( $ARGV[1] == 1 ) {

	my @roary_file = glob ("roary_temp/*");

	foreach my $c1 (@roary_file) {

		if ($c1 =~ /roary.group.txt/) {
			next;
		}

		system ("cp -r $c1 data_complete");
	}

	system ("perl src/FilePrepare.pl 1");
	system ("perl src/ComToPre.pl");
	system ("cp roary_temp/roary.group.txt Adjustment/cog/cog.group.txt");
	system ("perl src/GetMobileGOGs.pl 3 bootstrap");
	system ("perl src/OutFilter.pl $ARGV[2]");
	system ("perl src/GetMobileGOGs.pl 1 Adjustment/tmp/mobile.list");

} elsif ( $ARGV[1] == 2 ) {

	system ("mkdir data_draft");
	my @roary_file = glob ("roary_temp/*");
	my $tmp = & file_cp (@roary_file);

	system ("perl src/FilePrepare.pl 2");
	system ("perl src/SyntenyMummer.pl");
	system ("perl src/synteny1.pl");
	system ("perl src/synteny2.pl");
	system ("perl src/synteny3.pl");
	system ("perl src/DraToPre.pl");
	system ("cp roary_temp/roary.group.txt Adjustment/cog/cog.group.txt");
	system ("perl src/GetMobileGOGs.pl 3 bootstrap");
	system ("perl src/OutFilter.pl $ARGV[2]");
	system ("perl src/GetMobileGOGs.pl 1 Adjustment/tmp/mobile.list");
}




sub file_cp {

	foreach my $c1 (@_) {

		if ($c1 =~ /roary.group.txt/) {
			next;
		}

		my @fna_file = glob ("$c1/*.fna");

		if ( ! @fna_file) {
			die "Error: No fna file in $c1 !\n";
		}

		my $contig_num;

		open FNA, "$fna_file[0]";
		while (<FNA>) {
			if (/>/) {
				$contig_num++;
			}
		}
		close FNA;

		if ( $contig_num == 1 ) {
			system ("cp $c1 data_complete");
		} else {
			system ("cp $c1 data_draft");
		}
	}

	return 1;
}