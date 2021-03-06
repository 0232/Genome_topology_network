#!usr/bin/perl
use strict;

# Usage: perl gtn_CompleteOnly.pl <thread number> <cluster method> <outgroup>
# If there is no outgroup, then just: perl gtn_CompleteOnly.pl <thread number> <cluster method>
# <cluster method> should be either 'blast' or 'cdhit'
# Xiao Deng

my $threads = $ARGV[0];
my $meth = $ARGV[1];
my $outgroup = $ARGV[2];

if ( !($meth eq 'blast') and !($meth eq 'cdhit')) {
	print "$meth\n";
	print "Error! <cluster method> must be either blast or cdhit ! \n";
	die;
}

if(not -d "data_complete"){
	print "Error! Complete genome data file \"data_complete\"does not exist!";
	die;
}
my @data_check=glob("data_complete/*");

if(!@data_check){
	print "Error! No data in \"data_complete\" !";
	die;
}
mkdir "temp" unless -d "temp";

if(-d "temp/data"){
	system("rm -rf temp/data/*");
}else{
	mkdir "temp/data";
}

if(-d "temp/filter"){
	system("rm -rf temp/filter/*");
}else{
	mkdir "temp/filter";
}


system("perl src/FilePrepare.pl 1");
system("perl src/ComToPre.pl");
if ($meth eq 'blast') {
	system("perl src/FamCluster_blast.pl $threads");
}elsif ($meth eq 'cdhit') {
	system("perl src/FamCluster_cdhit.pl $threads");
}
system("perl src/GetMobileGOGs.pl 3 bootstrap");
system("perl src/GetMobileGOGs.pl 2 distance.meg");
system("perl src/OutFilter.pl $outgroup");
system("perl src/GetMobileGOGs.pl 1 Adjustment/tmp/mobile.list");


