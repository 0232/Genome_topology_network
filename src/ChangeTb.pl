
use iGraph;

die "perl $0 <Ortholog annotation file><GFF files Dir><gff>" unless ($#ARGV == 2);

# useage :perl ChangeTb.pl ./Adjustment/cog ./Adjustment/gff gff

my $a = iGraph->new($ARGV[0],$ARGV[1],$ARGV[2]); 
my @out = @{$a->getVaration('cog.nwk')};
open VAR,">./Adjustment/tmp/variation.out";
for my $k(@out)
{
	print ${$k}[0],"\t",${$k}[1],"\t",${$k}[4],"|",${$k}[5],"|",${$k}[2],"|",${$k}[3],"\n";
	my %tmp = ();
	my %var = %{${$k}[-1]};
	for my $k1(keys %var)
	{
		for my $k2(keys %{$var{$k1}})
		{
			if(exists $tmp{$k1}{$k2})
			{
				next;
			}
			else
			{
				for my $k3(keys %{$var{$k1}{$k2}})
				{
					print VAR ${$k}[0],"\t",${$k}[1],"\t",$k1,'-',$k2,"\t",$k3,"\t",$var{$k1}{$k2}{$k3},"\n";
				}
				$tmp{$k1}{$k2} = 1;
				$tmp{$k2}{$k1} = 1;
			}
		}
	}
	
}
exit;
#my @arr = qw/COG5651,COG2801,COG2963,COG1309,COG1960,COG0318,COG0596,COG1028,COG1024,COG0500,COG3321,COG0477,COG2141,COG0277,COG2409,COG2124,COG2114,COG3315,COG3547,COG3328,COG0515,COG1680,COG0657,COG1848/;
#my @files = glob "./gff/*";
#for my $f(@files)
#{
#	if($f =~ /add/)
#	{
#		next;
#	}
#	$f =~ s/\.\/gff\///;
#	$a->getFile("$f");##获取物种的GTN网络图,生成dot文件
#}
#exit;

#$a->getFile("m.H37Rv.gff");
#exit;
#print join("\n",@{$a->getSpecies()}),"\n"; ##获取物种列表，返回数组指针
#exit;
#print $a->getDistanceMatrix("mega","SimiMatrixUn");##输出物种间的距离矩阵
#exit;
#my @arr = ('m.H37Ra.gff','m.H37Rv.gff');
#my $test = $a->getMobile(@arr); ##得到mobile元件Hash的地址
#my $test = $a->getCore(@{$a->getSpecies()}); ##得到共有元件Array的地址
#my $test = $a->getCore("m.CDC1551.gff","m.F11.gff");
#my $test = $a->getStable("m.CDC1551.gff","m.F11.gff");##得到保守（不动）元件Array的地址
#my $test = $a->getStable(@{$a->getSpecies()});
#print join("\n",@{$test}),"\n";
#my $test = $a->getVaration(@{$a->getSpecies()});
my @sp = reverse @{$a->getSpecies()};
#for my $k(@{$a->getSpecies()})
#{
	#print $k,"\t";
	#for my $k1(keys %{${$test}{$k}})
#	{
		#print $k1,"\t";
		#my $a = abs(${${$test}{$k}{$k1}}{'A'} - ${${$test}{$k}{$k1}}{'B'})/${${$test}{$k}{$k1}}{'Both'};
		#my $b = abs(${${$test}{$k}{$k1}}{'Diffa'} -${${$test}{$k}{$k1}}{'Diffb'})/${${$test}{$k}{$k1}}{'Both'};
		#print $k,"\t",$k1,"\t",${${$test}{$k}{$k1}}{'A'},"\t",${${$test}{$k}{$k1}}{'B'},"\t",${${$test}{$k}{$k1}}{'Both'},"\t",${${$test}{$k}{$k1}}{'Diffa'},"\t",${${$test}{$k}{$k1}}{'Diffb'},"\n";
		#print sqrt($a**2 + $b**2);
		#print "\t";
##	}
	#print "\n";
#}
for(my $i=0;$i<=$#sp;$i++)
{
	print $sp[$i],"\t";
	for(my $j=0;$j<$i;$j++)
	{
		#print $sp[j];
		my %tmp = %{$a->getVaration_Bi_Un($a->getGraph($sp[$i]),$a->getGraph($sp[$j]))};
		print $tmp{'A'}+$tmp{'B'}+$tmp{'C1'}+$tmp{'C2'},"\t";	
		
	}
	print "\n";
}
exit;

#foreach my $key (keys %{$test}) 
#{
#	print $key,"\t",${$test}{$key},"\n";
#}

#exit;

##my %tmp = %{$a->getNeighbor("a.RSK2980.gff","COG0582")};###输出指定物种，指定orthology的邻接orthology

#foreach my $key (keys %tmp) {
#	print $key,"\n";
#}

#print "\n";

#exit;

#%tmp = %{$a->dumpNeighbor("m.F11.gff")};###输出指定物种，所有orthology的邻接orthology
#foreach my $key (keys %tmp) {
#	foreach my $key_1 (keys %{$tmp{$key}}) {
#		print $key,"\t",$key_1,"\t",$tmp{$key}{$key_1},"\n";
#	}
#}

