#!usr/bin/perl
use strict;


my $threads = $ARGV[0];

&get_all_pep('./myva');
system ("cd-hit -i Adjustment/tmp/all.pep -c 0.9 -g 1 -d 60 -T $threads -M 0 -o Adjustment/tmp/cdhit.out");


my %clu = &get_clu("Adjustment/tmp/cdhit.out.clstr");
my %seq = &get_seq("Adjustment/tmp/all.pep");

open CLU, ">Adjustment/tmp/all_pep.cdhit.clu";
open PEP, ">Adjustment/tmp/cluster_represent.pep";

foreach my $c1 (keys %clu) {
	if ($clu{$c1} =~/\s/) {
		my @data = split(/\s/, $clu{$c1});
		$data[0] =~s/.*\|//;
		print PEP ">$c1\n$seq{$data[0]}\n";
		print CLU "$c1: $clu{$c1}\n";
	}
}
close CLU;
close PEP;


system ("diamond makedb  --in Adjustment/tmp/COG.pep  -d Adjustment/tmp/COG.pep");
system ("diamond blastp -p $threads --sensitive -d Adjustment/tmp/COG.pep -q Adjustment/tmp/cluster_represent.pep -o Adjustment/tmp/all_pep.m8");

&get_cog_table(%clu);
system("cp Adjustment/tmp/diamond.cog.txt Adjustment/cog/cog.group.txt");





sub get_all_pep {
	my @pep_file = glob("Adjustment/pep/*.pep");
	open OUT, ">Adjustment/tmp/all.pep";
	open COG, ">Adjustment/tmp/COG.pep";

	foreach my $c1 (@pep_file) {
		open IN, "$c1";
		while (<IN>) {
			if (/>/) {
				my @data = split(/\|/, $_);
				print OUT ">$data[3]\n";
				next;
			}
			print OUT $_;
		}
		close IN;
	}

	my %cog_tab;
	open COG_TAB, "src/cog.id.tab" or die "Error: file cog.id.tab is not found!\n";
	while (<COG_TAB>) {
		$_ =~tr/\r\n//d;
		my @data = split(/\t/, $_);
		$cog_tab{$data[1]} = $data[0];
	}
	close COG_TAB;

	my $num = 1;
	my $id;
	my %cog;
	open MYVA, $_[0] or die "Error: file myva is not found!\n";
	while (<MYVA>) {
		$_ =~tr/\r\n//d;
		if (/>(.*)/) {
			$id = ">$cog_tab{$1}_$num";
			$num++;
			next;
		}
		$cog{$id} .= $_;
	}
	close MYVA;

	foreach my $c1 (keys %cog) {
		if ($c1 =~/COG/) {
			print COG "$c1\n$cog{$c1}\n";
		}
	}
	close COG;
	close OUT;
}



sub get_seq {
	my %hash;
	my $id;

	open SEQ, $_[0] or die "Error: all.pep file is not found!\n";
	while (<SEQ>) {
		$_ =~tr/\r\n//d;
		if (/>(.*)/) {
			$id = $1;
			next;
		}
		$hash{$id} .= $_;
	}
	close SEQ;

	return %hash;
}


sub get_clu {
	my %clu;
	my $id;

	open IN, $_[0] or die "Error: cdhit result file is not found!\n";
	while (<IN>) {
		$_ =~tr/\r\n//d;
		if ($_ =~/>.*\s(\d+)/ and $_ !~/\,/) {
			$id = "cluster".$1;
			next;
		}
		$_ =~/>(.*__(.*))\.\.\./;
		my $gene = $2.'|'.$1;
		if (!$clu{$id}) {
			$clu{$id} = $gene;
		}else {
			$clu{$id} .= " ".$gene;
		}
	}
	close IN;

	return %clu;
}


sub get_cog_table {
	my %clu = @_;
	my %anno;

	open IN, "Adjustment/tmp/all_pep.m8" or die;
	while (<IN>) {
		$_ =~tr/\r\n//d;
		my @data = split(/\t/, $_);

		if (!$anno{$data[0]}) {
			$anno{$data[0]} = $data[1];
		}else {
			$anno{$data[0]} .= "\t".$data[1];
		}
	}
	close IN;

	open OUT, ">Adjustment/tmp/diamond.cog.txt";
	foreach my $c1 (keys %clu) {
		if ($anno{$c1}) {
			if ($anno{$c1} =~/\t/) {
				my @data_gene = split(/\t/, $anno{$c1});
				$data_gene[0] =~s/_\d+//;
				print OUT "$data_gene[0]: $clu{$c1}\n";
			}else {
				$anno{$c1} =~s/_\d+//;
				print OUT "$anno{$c1}: $clu{$c1}\n";
			}
		}
	}
	close OUT;
}
