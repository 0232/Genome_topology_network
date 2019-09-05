#!/usr/bin/perl -w
use warnings;
use strict;
use File::Copy;

 
if($#ARGV != 4)
{
	print "Usage:perl $0 <species identity> <mRNA file> <protein file> <genome file> <GFF file>\n";
	exit;
}

for my $file(@ARGV[1..3])
{
	if(!&checkPara($file))
	{
		print $file,"not exists![Error!]\n";
		exit;
	}
}

my %genome_id = &checkFiles($ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4]);

if(not -d 'Adjustment')
{
	mkdir 'Adjustment';
}
chdir 'Adjustment';
if(not -d 'gff')
{
	mkdir 'gff';
}
if(not -d 'mRNA')
{
	mkdir 'mRNA';
}
if(not -d 'pep')
{
	mkdir 'pep';
}
if(not -d 'genome')
{
	mkdir 'genome';
}
if(not -d 'tmp')
{
	mkdir 'tmp';
}
if(not -d 'cog')
{
	mkdir 'cog';
}
chdir "..";
copy("$ARGV[1]","./Adjustment/mRNA/$ARGV[0].mRNA");
copy("$ARGV[2]","./Adjustment/pep/$ARGV[0].pep");
copy("$ARGV[3]","./Adjustment/genome/$ARGV[0].fa");
copy("$ARGV[4]","./Adjustment/gff/$ARGV[0].gff");

open REL,">>./Adjustment/species_id.tab";
foreach (keys %genome_id)
{
print REL "$_\t$ARGV[0]\n";
}
print "Adjustment done!\n";



sub checkPara()
{
	my $file = shift;
	if(-e $file)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub checkFiles()
{
	my $mRNA = shift;
	my $pep = shift;
	my $genome = shift;
	my $gff = shift;
	open GENOME,"<$genome";
	my $genome_id = '';
	my %genome_id=();
	
	while(my $line = <GENOME>)
	{
		if($line =~ />/)
		{
			my @eles = split/\|/,$line;
			$genome_id{$eles[3]} =1 ;
		}
	}
	close GENOME;
	open GFF,"<$gff";
	my %hash_gff = ();
	my %hash_pep = ();
	while(my $line = <GFF>)
	{
		if($line=~/pseudo\=true/){
			next;
			}
		if($line =~ /^#/)
		{
			next;
		}
		my @eles = split/\t/,$line;
		if( ! $genome_id{$eles[0]})
		{
#			die "Genome sequence ID is not identical to GFF file![Error!]\n";
			print "Genome sequence ID is not identical to GFF file![Error!]\n";
		}
		if($eles[2] eq 'CDS' and $line=~/pseudo\=true/){                      # some CDS is pseudo , skip these
			next;
			}
		if($eles[2] eq 'CDS')
		{
			$eles[8] =~ /protein_id=(.+?)\;/;
			$hash_pep{$1} = 1;
			if($eles[6] eq '+')
			{
				$hash_gff{$eles[3]}{$eles[4]} = 1;
			}
			else
			{
				$hash_gff{$eles[4]}{$eles[3]} = 1;
			}			
		}
	}
	close GFF;
	open PEP,"<$pep";
	while(my $line = <PEP>)
	{
		if($line !~ /^>/)
		{
			next;
		}
		else
		{
			my @eles = split/\|/,$line;
			my $protein_id = $eles[3];
			if(!exists $hash_pep{$protein_id})
			{
				print "Proetin[$protein_id] is not contained in GFF![Warning!]\n";
			}
		}
	}
	close PEP;
	open MRNA,"<$mRNA";
	while(my $line = <MRNA>)
	{
		if($line !~ /^>/)
		{
			next;
		}
		my @eles = split/\|/,$line;
		$eles[-1] =~ /(\d+-\d+)/;
		my @position = split/\-/,$1;
		if(!exists $hash_gff{$position[0]}{$position[1]})
		{
			print "$line is not contained in GFF![Warning!]\n";
		}
	}
	close MRNA;
	return %genome_id;
}