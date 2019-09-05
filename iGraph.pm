#!c:/perl64/bin/perl.exe
#! /usr/bin/perl -w
package iGraph;
=head1 NAME

    iGraph

=head1 DESCRRIPTION

=head1 AUTHOR

       Jian-lei Gu <jianleigu@gmail.com>

=head1 COPYRIGHT

=head1 LICENSE

=cut
use strict;
use Exporter;
use vars qw{@ISA};
use Carp;

sub new{
	my $object = shift;
	my $class = ref($object) || $object;
	my $self = {};
	bless $self, $class;
	$self->_initialize(@_);
	return $self;
}

sub _initialize{
	my $self = shift;
	$self -> copy(load($self,@_));
}

sub copy () {
	my $self = shift;
	my $add = shift;
	%{$self->{'GTN'}} = ();
	foreach my $ispecie (keys %{$add}) {
		foreach my $forward (keys %{${$add}{$ispecie}}) {
			foreach my $backward (keys %{${$add}{$ispecie}{$forward}}) {
				$self->{'GTN'}{$ispecie}{$forward}{$backward} = ${$add}{$ispecie}{$forward}{$backward};
				$self->{'GTN'}{$ispecie}{$backward}{$forward} = ${$add}{$ispecie}{$backward}{$forward};
			}
		}
	}
	return 1;
}

sub load () {
	my $self = shift;
	my $ortholog_file = shift;
	my $gff_dir = shift;
	my $file_type = shift;
	my %gff = ();
	my $ortholog_add = &ortholog_load($ortholog_file);
	$self->{'Orth'} = $ortholog_add;
	my @gff_files = @{&Readingdir($gff_dir,$file_type)};
	if ($gff_dir !~/\/$|\/$/) {
		if ($gff_dir=~/\\/) {
			$gff_dir = $gff_dir."\\";
		}else{
			$gff_dir = $gff_dir."/";
		}
	}
	map {
		&gff_load($gff_dir.$_,\%gff,$ortholog_add)
		} @gff_files;
	my $network = &gff_rank(\%gff);
	my $graph = &Node_map($network,$ortholog_add);
	return $graph;
}

sub refreshSimiMatrix { #################
	my $self = shift;
	my $comm = shift || 0;
	my @species = keys %{$self->{'GTN'}};
	my $falg = 0;
	for (my $i = 0; $i < $#species; $i++) {
		for (my $j = $i +1; $j <= $#species; $j++) {
			if (not (exists $self->{'SimiMatrixUn'}{$species[$i]}{$species[$j]})) {
				$self->{'SimiMatrixUn'}{$species[$i]}{$species[$j]} = &unWeight($self,$species[$i],$species[$j]);
				$self->{'SimiMatrixUn'}{$species[$j]}{$species[$i]} = $self->{'SimiMatrixUn'}{$species[$i]}{$species[$j]};
				$self->{'SimiMatrixEn'}{$species[$i]}{$species[$j]} = &Weight($self,$species[$i],$species[$j]);
				$self->{'SimiMatrixEn'}{$species[$j]}{$species[$i]} = $self->{'SimiMatrixEn'}{$species[$i]}{$species[$j]};
				$falg= 1;
			}else{
				next unless ($comm == 1);
				$self->{'SimiMatrixUn'}{$species[$i]}{$species[$j]} = &unWeight($self,$species[$i],$species[$j]);
				$self->{'SimiMatrixUn'}{$species[$j]}{$species[$i]} = $self->{'SimiMatrixUn'}{$species[$i]}{$species[$j]};
				$self->{'SimiMatrixEn'}{$species[$i]}{$species[$j]} = &Weight($self,$species[$i],$species[$j]);
				$self->{'SimiMatrixEn'}{$species[$j]}{$species[$i]} = $self->{'SimiMatrixEn'}{$species[$i]}{$species[$j]};
				$falg= 1;
			}
		}
	}
	return $falg;
}


sub refreshBootstrapMatrix { #################
	my $self = shift;
	my $comm = shift || 0;
	my @species = keys %{$self->{'GTN'}};
	my $falg = 0;
	my $tmp = &getOrtho($self);
	my @ortho = @{$tmp};
	my @rand_ortho = ();
	for(my $i=0;$i<=$#ortho;$i++)
	{
		my $j = int(rand($#ortho+1));
		push @rand_ortho,$ortho[$j];
	}
	for (my $i = 0; $i < $#species; $i++) {
		for (my $j = $i +1; $j <= $#species; $j++) {
			if (not (exists $self->{'BootstrapMatrix'}{$species[$i]}{$species[$j]})) {
				$self->{'BootstrapMatrix'}{$species[$i]}{$species[$j]} = &BootStrapUnWeight($self,$species[$i],$species[$j],@rand_ortho);
				$self->{'BootstrapMatrix'}{$species[$j]}{$species[$i]} = $self->{'BootstrapMatrix'}{$species[$i]}{$species[$j]};
				$falg= 1;
			}else{
				next unless ($comm == 1);
				$self->{'BootstrapMatrix'}{$species[$i]}{$species[$j]} = &BootStrapUnWeight($self,$species[$i],$species[$j],@rand_ortho);
				$self->{'BootstrapMatrix'}{$species[$j]}{$species[$i]} = $self->{'BootstrapMatrix'}{$species[$i]}{$species[$j]};
				$falg= 1;
			}
		}
	}
	return $falg;
}


sub unWeight { #################
	my $self = shift;
	my $genome1 = shift;
	my $genome2 = shift;
	my %graph = %{&icombine($self,$genome1,$genome2)};
	my @tmp = ();
	foreach my $forward (keys %graph) {
		my $Both = 0;my $A1 = 0;my $A0 = 0;
		map {						
					if ($graph{$forward}{$_} eq "A") {
						$A1 += 1;
					}elsif ($graph{$forward}{$_} eq "B") {
						$A0 += 1;
					}elsif ($graph{$forward}{$_} eq "C") {
						$Both += 1;
					}else {
						warn "UnKnown Species $graph{$forward}{$_}.\n";
					}			
			} keys %{$graph{$forward}};
				my $out = (2*$Both)/((2*$Both) + $A0 +$A1);
				push @tmp, $out;

			
	}
	return (1 - &mean(@tmp));
	#return &sum(@tmp);
}


sub BootStrapUnWeight()
{
	my $self = shift;
	my $genome1 = shift;
	my $genome2 = shift;
	my @ortho = @_;
	my %graph = %{&icombine($self,$genome1,$genome2)};
	my @ortho_refine = ();
	my %hash_refine = ();
	my @tmp = ();
	for my $k(@ortho)
	{
		if(exists $graph{$k})
		{
			push @ortho_refine,$k;
		}
		else
		{
			next;
		}
	}
	
	foreach my $forward (@ortho_refine) {
		my $Both = 0;my $A1 = 0;my $A0 = 0;
		map {					
					if ($graph{$forward}{$_} eq "A") {
						$A1 += 1;
					}elsif ($graph{$forward}{$_} eq "B") {
						$A0 += 1;
					}elsif ($graph{$forward}{$_} eq "C") {
						$Both += 1;
					}else {
						warn "UnKnown Species $graph{$forward}{$_}.\n";
					}			
			} keys %{$graph{$forward}};
				my $out = (2*$Both)/((2*$Both) + $A0 +$A1);
				push @tmp, $out;			
	}
	return (1 - &mean(@tmp));	
}




sub RandomDistance { #################
	my $self = shift;
	my $genome1 = shift;
	my $genome2 = shift;
	my %graph = %{&icombine($self,$genome1,$genome2)};
	my @tmp = ();
	foreach my $forward (keys %graph) {
		my $Both = 0;my $A1 = 0;my $A0 = 0;
		map {						
					if ($graph{$forward}{$_} eq "A") {
						$A1 += 1;
					}elsif ($graph{$forward}{$_} eq "B") {
						$A0 += 1;
					}elsif ($graph{$forward}{$_} eq "C") {
						$Both += 1;
					}else {
						warn "UnKnown Species $graph{$forward}{$_}.\n";
					}			
			} keys %{$graph{$forward}};
				my $out = (2*$Both)/((2*$Both) + $A0 +$A1);
				push @tmp, $out;

			
	}
	my @num = ();
	for(my $i=0;$i<1000;$i++)
	{
		my $j = rand($#tmp+1);
		push @num,$tmp[$j];
	}
	return (1 - &mean(@num));
	#return &sum(@tmp);
}

sub PartUn()
{
	my $self = shift;
	my $genome1 = shift;
	my $genome2 = shift;
	my $ele = shift;
	my @arr = @{$ele};
	my %graph = %{&icombine($self,$genome1,$genome2)};
	my @tmp = ();
	foreach my $forward (@arr) {
		my $Both = 0;my $A1 = 0;my $A0 = 0;
		map {
				if ($graph{$forward}{$_} eq "A") {
					$A1 += 1;
				}elsif ($graph{$forward}{$_} eq "B") {
					$A0 += 1;
				}elsif ($graph{$forward}{$_} eq "C") {
					$Both += 1;
				}else {
					warn "UnKnown Species $graph{$forward}{$_}.\n";
				}
			} keys %{$graph{$forward}};
		push @tmp, (2*$Both)/((2*$Both) + $A0 +$A1);
	}
	return (1 - &mean(@tmp));
}



sub Weight () {
	my $self = shift;
	my $genome1 = shift;
	my $genome2 = shift;
	my %graph = %{&icombine($self,$genome1,$genome2)};
	my %graph1= %{&getGraph($self,$genome1)};
	my %graph2= %{&getGraph($self,$genome2)};
	my @tmp = ();
	foreach my $forward (keys %graph) {
		my $Both = 0;my $A1 = 0;my $A0 = 0;my $Div = 0;
		map {
				if ($graph{$forward}{$_} eq "A") {
					$A1 += $graph1{$forward}{$_};
				}elsif ($graph{$forward}{$_} eq "B") {
					$A0 += $graph2{$forward}{$_};
				}elsif ($graph{$forward}{$_} eq "C") {
					$Both += $graph1{$forward}{$_};
					$Both += $graph2{$forward}{$_};
					$Div += &min($graph1{$forward}{$_},$graph2{$forward}{$_});				

				}else {
					warn "UnKnown Species $graph{$forward}{$_}.\n";
				}
			} keys %{$graph{$forward}};
		push @tmp, $Div*2/($Both + $A0 +$A1);
	}
	return (1 - &mean(@tmp));
}

sub getMobile () {
	my $self = shift;
	my %out = ();
	if (scalar(@_) >= 2) {
		my $xy = 0;
		for (my $i = 0; $i < $#_; $i++) {
			for (my $j = $i+1; $j <= $#_; $j++) {
				$xy++;
				my %tmp = %{$self->getMobile_Bi($_[$i],$_[$j])};
				foreach my $key (keys %tmp) {
					if (not (exists $out{$key})) {
						$out{$key} = $tmp{$key};
					}else{
						$out{$key} += $tmp{$key};
					}
				}
			}
		}
		map {
			$out{$_} = int ($out{$_}/$xy);
		}keys %out;
		for my $k(keys %out)
		{
			if($out{$k} == 0)
			{
				$out{$k} = 1;
			}
		}
		return \%out;
	}else{
		return 0;
	}
}

sub getStable () {
	my $self = shift;
	my %out = ();
	my @output;
	if (scalar(@_) >= 2) {
		my $xy = 0;
		for (my $i = 0; $i < $#_; $i++) {
			for (my $j = $i+1; $j <= $#_; $j++) {
				$xy++;
				my %tmp = %{$self->getStable_Bi($_[$i],$_[$j])};
				foreach my $key (keys %tmp) {
					if (not (exists $out{$key})) {
						$out{$key} = 1;
					}else{
						$out{$key} ++;
					}
				}
			}
		}
		map {
			if ($out{$_} == $xy) {
				push @output, $_;
			}
		}keys %out;
	}else{
		push @output, 0;
	}
	return \@output;
}

sub getCore () {
	my $self = shift;
	my %graph = %{&Xcombine($self,@_)};
	my %tmp = ();
	map {$tmp{$_}=1} @_;
	my $num = scalar(keys %tmp);
	%tmp = ();
	foreach my $forward (keys %graph) {
		my $temp = scalar(keys %{$graph{$forward}});
		if ($temp ==  $num) { #######################÷ÿ∏¥≥È—˘Bug
			$tmp{$forward} = 1;
		}
	}
	my @out = keys %tmp;
	return \@out;
}


sub getVaration()
{
	my $self = shift;
	my $nwk = shift;
	open NWK,"<cog.nwk" or die $!;
	my @out = ();
	while(my $line = <NWK>)
	{
		chomp $line;
		$line =~ s/\:\d+\.\d+//g;
		my $ele = $line;
		$ele =~ s/\;|\r//g;
		my @eles = split/\,/,$ele;
		my @stack = ();
		for my $e(@eles)
		{
			my %tmp = ();
			my @l = $e =~ /\(/g;
			my @r = $e =~ /\)/g;
			$e =~ s/\(|\)//g;
			$tmp{'l'} = $#l + 1;
			$tmp{'r'} = $#r + 1;
			$tmp{'id'} = $e;
			push @stack,\%tmp;
		}
		my $s = shift @stack;
		my @queue = ();
		push @queue,$s;
		my %comput = ();
		my $j = 1;
		while($#stack != -1 && $#queue != -1)
		{
			my $pre = pop @queue;
			my $node = shift @stack;
			if(${$node}{'l'} == 0 && ${$node}{'r'} == 0)
			{
				$comput{$j}{'A'} = ${$pre}{'id'};
				$comput{$j}{'B'} = ${$node}{'id'};
				$j++;
				${$pre}{'id'} .= '+'.${$node}{'id'};			
				push @queue,$pre;
			}
			elsif(${$pre}{'l'} > 0 && ${$node}{'r'} >0)
			{
				${$pre}{'l'} -= 1;
				${$node}{'r'} -= 1;
				if(${$pre}{'l'} == 0 && ${$node}{'r'} != 0)
				{
					$comput{$j}{'A'} = ${$pre}{'id'};
					$comput{$j}{'B'} = ${$node}{'id'};
					$j++;		
					${$node}{'id'} .= '+'.${$pre}{'id'};
					unshift @stack,$node;
				}
				elsif(${$node}{'r'} == 0 && ${$pre}{'l'} != 0)
				{
				
					$comput{$j}{'A'} = ${$pre}{'id'};
					$comput{$j}{'B'} = ${$node}{'id'};
					$j++;	
					${$pre}{'id'} .= '+'.${$node}{'id'};
					push @queue,$pre;
				}	
				elsif(${$pre}{'l'} ==0 && ${$node}{'r'} == 0)
				{
					if($#queue == -1 && $#stack == -1)
					{
						$comput{$j}{'A'} = ${$pre}{'id'};
						$comput{$j}{'B'} = ${$node}{'id'};
						$j++;	
						${$pre}{'id'} .= '+'.${$node}{'id'};
						push @queue,$pre;
					}
					else
					{
						$comput{$j}{'A'} = ${$pre}{'id'};
						$comput{$j}{'B'} = ${$node}{'id'};
						$j++;
						my $t = pop @queue;
						$comput{$j}{'A'} = ${$pre}{'id'}.'+'.${$node}{'id'};
						$comput{$j}{'B'} = ${$t}{'id'};
						$j++;
						${$t}{'id'} .= '+'.${$node}{'id'}.'+'.${$pre}{'id'};
						push @queue,$t;
					}
				}
				
			}
			elsif(${$pre}{'l'} >0 && ${$node}{'l'} >0)
			{
				push @queue,$pre;
				push @queue,$node;
			}
		}
		
		for my $k(sort{$a<=>$b} keys %comput)
		{
			my @left = split/\+/,$comput{$k}{'A'};
			my @right = split/\+/,$comput{$k}{'B'};
			my %g1 = %{&GetCoreGraph($self,\@left)};
			my %g2 = %{&GetCoreGraph($self,\@right)};
			my %tmp = %{&getVaration_Bi_Un($self,\%g1,\%g2)};
			my @out1 = ($comput{$k}{'A'},$comput{$k}{'B'},$tmp{'A'},$tmp{'B'},$tmp{'C1'},$tmp{'C2'},$tmp{'G1'});
			my @out2 = ($comput{$k}{'B'},$comput{$k}{'A'},$tmp{'B'},$tmp{'A'},$tmp{'C2'},$tmp{'C1'},$tmp{'G2'});
			push @out,\@out1;
			push @out,\@out2;
			
		}
		
	}
	return \@out;
}
sub getStable_Bi () {
	my $self = shift;
	my $genome1 = shift;
	my $genome2 = shift;
	my %graph = %{&icombine($self,$genome1,$genome2)};
	my %graph1= %{&getGraph($self,$genome1)};
	my %graph2= %{&getGraph($self,$genome2)};
	my %tmp = ();
	foreach my $forward (keys %graph) {
		my $Both = 0;my $A1 = 0;my $A0 = 0;
		map {
				if ($graph{$forward}{$_} eq "A") {
					$A1 += $graph1{$forward}{$_};
				}elsif ($graph{$forward}{$_} eq "B") {
					$A0 += $graph2{$forward}{$_};
				}elsif ($graph{$forward}{$_} eq "C") {
					$Both += abs( $graph1{$forward}{$_} - $graph2{$forward}{$_} );
				}else {
					warn "UnKnown Species $graph{$forward}{$_}.\n";
				}
			} keys %{$graph{$forward}};
		my $temp = $A0+$A1+$Both;
		if ($temp == 0) {
			$tmp{$forward} = $temp;
		}
	}
	return \%tmp;
}

sub getMobile_Bi () {
	my $self = shift;
	my $genome1 = shift;
	my $genome2 = shift;
	my %graph = %{&icombine($self,$genome1,$genome2)};
	my %graph1= %{&getGraph($self,$genome1)};
	my %graph2= %{&getGraph($self,$genome2)};
	my %tmp = ();
	foreach my $forward (keys %graph) {
		my $Both = 0;my $A1 = 0;my $A0 = 0;
		map {
				if ($graph{$forward}{$_} eq "A") {
					$A1 += $graph1{$forward}{$_};
				}elsif ($graph{$forward}{$_} eq "B") {
					$A0 += $graph2{$forward}{$_};
				}elsif ($graph{$forward}{$_} eq "C") {
					$Both += abs( $graph1{$forward}{$_} - $graph2{$forward}{$_} );#
				}else {
					warn "UnKnown Species $graph{$forward}{$_}.\n";
				}
			} keys %{$graph{$forward}};
		my $temp = $A0+$A1+$Both;
		if ($temp != 0) {
			$tmp{$forward} = $temp;
		}
	}
	return \%tmp;
}


sub getVaration_Bi_Un()
{
	my $self = shift;
	my $genome1 = shift;
	my $genome2 = shift;
	my %graph1= %{$genome1};
	my %graph2= %{$genome2};
	my $a = 0;
	my $b = 0;
	my $c1 = 0;
	my $c2 = 0;
	my %g1 = ();;
	my %g2 = ();
	for my $forward(keys %graph1)
	{
		for my $backward(keys %{$graph1{$forward}})
		{
			#$t1 += $graph1{$forward}{$backward};
			if(exists $graph2{$forward}{$backward})
			{
				#$a += abs($graph2{$forward}{$backward} - $graph1{$forward}{$backward});
				if($graph1{$forward}{$backward} > $graph2{$forward}{$backward})
				{
					$a += $graph1{$forward}{$backward} - $graph2{$forward}{$backward};
					$g1{$forward}{$backward}{'duplicate'} = $graph1{$forward}{$backward} - $graph2{$forward}{$backward};
					$g2{$forward}{$backward}{'delete'} = $graph1{$forward}{$backward} - $graph2{$forward}{$backward};
				}
				elsif($graph1{$forward}{$backward} == $graph2{$forward}{$backward})
				{
					next;
				}
				else
				{
					$b += $graph2{$forward}{$backward} - $graph1{$forward}{$backward};
					$g1{$forward}{$backward}{'delete'} = $graph2{$forward}{$backward} - $graph1{$forward}{$backward};
					$g2{$forward}{$backward}{'duplicate'} = $graph2{$forward}{$backward} - $graph1{$forward}{$backward};
				}
			}
			else
			{
				$c1 += $graph1{$forward}{$backward};
				$g1{$forward}{$backward}{'unique'} = $graph1{$forward}{$backward};
				$g2{$forward}{$backward}{'miss'} = $graph1{$forward}{$backward};
			}
		}
	}
	for my $forward(keys %graph2)
	{
		for my $backward(keys %{$graph2{$forward}})
		{
			#$t2 +=  $graph2{$forward}{$backward};
			if(exists $graph1{$forward}{$backward})
			{
				next;
			}
			else
			{
				$c2 += $graph2{$forward}{$backward};
				$g1{$forward}{$backward}{'miss'} = $graph2{$forward}{$backward};
				$g2{$forward}{$backward}{'unique'} = $graph2{$forward}{$backward};
			}
		}
	}
	my %tmp = ();
	$tmp{'A'} = $a/2; $tmp{'B'} = $b/2; $tmp{'C1'} = $c1/2;$tmp{'C2'} = $c2/2;$tmp{'G1'} = \%g1;$tmp{'G2'} = \%g2;
	return \%tmp;	
	
}

sub GetCoreGraph()
{
	my $self = shift;
	my $sp = shift;
	my $first = shift @{$sp};
	my %core = %{&getGraph($self,$first)};
	if($#{$sp} != -1)
	{
		for my $s(@{$sp})
		{
			my %graph = %{&getGraph($self,$s)};
			my %tmp = ();
			for my $forward(keys %graph)
			{
				for my $backward(keys %{$graph{$forward}})
				{
					if(exists $core{$forward}{$backward})
					{
						$tmp{$forward}{$backward} = &min($core{$forward}{$backward},$graph{$forward}{$backward});
					}
				}
			}
			
			%core = %tmp;
		}
	}
	return \%core;
}

sub Xcombine () {
	my $self = shift;
	my %graph;
	foreach my $ispecie (@_) {
		foreach my $forward (keys %{$self->{'GTN'}{$ispecie}}) {
			$graph{$forward}{$ispecie} = 1;
			foreach my $backward (keys %{$self->{'GTN'}{$ispecie}{$forward}}) {
				$graph{$backward}{$ispecie} = 1;
			}
		}
	}
	return \%graph;
}

sub getNeighbor () {
	my $self = shift;
	my $ispecie = shift;
	my $iNode = shift;
	if (exists $self->{'GTN'}{$ispecie}{$iNode}) {
		return \%{$self->{'GTN'}{$ispecie}{$iNode}};
	}else{
		return 0;
	}
}

sub dumpNeighbor () {
	my $self = shift;
	my $ispecie = shift;
	if (exists $self->{'GTN'}{$ispecie}) {
		return \%{$self->{'GTN'}{$ispecie}};
	}else{
		return 0;
	}
}

sub getSpecies () {
	my $self = shift;
	my @ele = keys %{$self->{'GTN'}};
	return \@ele;
}

sub getOrtho()
{
	my $self = shift;
	my %hash_ortho = ();
	my %hash_tmp = %{$self->{'Orth'}};
	for my $k(keys %hash_tmp)
	{
		for my $k1(keys %{$hash_tmp{$k}})
		{
			if(exists $hash_ortho{$hash_tmp{$k}{$k1}})
			{
				next;
			}
			else
			{
				$hash_ortho{$hash_tmp{$k}{$k1}} = 1;
			}
		}
	}
	my @ortho = keys %hash_ortho;
	return \@ortho;
}

sub dumpDisDegree () {
	my $self = shift;
	my $ispecie = shift;
	my %tmp = ();
	if (exists $self->{'GTN'}{$ispecie}) {
		foreach my $iNode (keys %{$self->{'GTN'}{$ispecie}}) {
			$tmp{$iNode} = scalar(keys %{$self->{'GTN'}{$ispecie}{$iNode}});
		}
		return \%tmp;
	}else{
		return 0;
	}
}

sub icombine () {
	my $self = shift;
	my $genome1 = shift;
	my $genome2 = shift;
	my %graph1= %{&getGraph($self,$genome1)};
	my %graph2= %{&getGraph($self,$genome2)};
	my %graph = ();
	foreach my $forward (keys %graph1) {
		foreach my $backward (keys %{$graph1{$forward}}) {
			$graph{$forward}{$backward} = "A";
			#$graph{$backward}{$forward} = "A";
		}
	}
	foreach my $forward (keys %graph2) {
		foreach my $backward (keys %{$graph2{$forward}}) {
			if (not (exists $graph{$forward}{$backward})) {
				$graph{$forward}{$backward} = "B";
				#$graph{$backward}{$forward} = "B";
			}elsif ($graph{$forward}{$backward} eq "A") {
				$graph{$forward}{$backward} = "C";
				#$graph{$backward}{$forward} = "C";
			}
		}
	}
	return \%graph;
}

sub iMerge () {
	my $self = shift;
	my @ispecie = @_;
	my %graph = ();
	if (scalar(@ispecie) >= 2) {
		%graph= %{&getGraph($self,$ispecie[0])};
		for (my $i = 1; $i <= $#ispecie; $i ++) {
			my $tmp = &getGraph($self,$ispecie[$i]);
			foreach my $forward (keys %{$tmp}) {
				foreach my $backward (keys %{${$tmp}{$forward}}) {
					$graph{$forward}{$backward} = 1;
					$graph{$backward}{$forward} = 1;
				}
			}
		}
		return \%graph;
	}else{
		return &getGraph($self,$ispecie[0]);
	}
}

sub getGraph () {
	my $self = shift;
	my $Genome = shift;
	die "do not exists $Genome!\n" unless (exists $self->{'GTN'}{$Genome});
	return \%{$self->{'GTN'}{$Genome}};
}

sub getFile () {
	my $self = shift;
	my $ispecie = shift;
	my $ifile = shift || $ispecie.".dot";
	my %tmp = ();
	return 0 unless (exists $self->{'GTN'}{$ispecie});
	open (OUT, ">".$ifile ) or die "$!";
	print $ispecie,scalar keys %{$self->{'GTN'}{$ispecie}},"\n";
	foreach my $forward (keys %{$self->{'GTN'}{$ispecie}}) {
		foreach my $backward (keys %{$self->{'GTN'}{$ispecie}{$forward}}) {
			next unless (not (exists $tmp{$forward}{$backward}));
			print OUT ($forward," =",$self->{'GTN'}{$ispecie}{$forward}{$backward},"= ",$backward,"\n");
			$tmp{$forward}{$backward} = 1;
			$tmp{$backward}{$forward} = 1;
		}
	}
	close OUT or die "$!";
	return 1;
}

sub Node_map () {
	my $network = shift;##gff: gene position.
	my $ortholog = shift;
	my %graph = ();
	foreach my $ispecie (keys %{$network}) {
		#my %tmp = ();
		foreach my $forward (keys %{${$network}{$ispecie}}) {
			foreach my $backward (keys %{${$network}{$ispecie}{$forward}}) {
				if ( (exists ${$ortholog}{$ispecie}{$forward}) && (exists ${$ortholog}{$ispecie}{$backward}) ) {## forward,backward both have been annotated by othology
					if (not (exists $graph{$ispecie}{${$ortholog}{$ispecie}{$forward}}{${$ortholog}{$ispecie}{$backward}})) {
						$graph{$ispecie}{${$ortholog}{$ispecie}{$forward}}{${$ortholog}{$ispecie}{$backward}} = ${$network}{$ispecie}{$forward}{$backward};
						#$graph{$ispecie}{${$ortholog}{$ispecie}{$backward}}{${$ortholog}{$ispecie}{$forward}} = ${$network}{$ispecie}{$backward}{$forward};
					}else{

							$graph{$ispecie}{${$ortholog}{$ispecie}{$forward}}{${$ortholog}{$ispecie}{$backward}} += ${$network}{$ispecie}{$forward}{$backward};
							#$graph{$ispecie}{${$ortholog}{$ispecie}{$backward}}{${$ortholog}{$ispecie}{$forward}} += ${$network}{$ispecie}{$backward}{$forward};
						
					}
				}else{
					die "Graph::Node_map:$ispecie $forward $backward seems lost in ortholog group.\n";
				}
			}
		}
	}
	return \%graph;
}

sub ortholog_load () {##orthology{species}{protein id}
	my $ifile = shift;
	my @files = glob "$ifile/*";
	my %ortholog = ();
	for my $file(@files)
	{
		my $line;
		open (IN,$file) or die "$!";
		while ($line = <IN>) {
			chomp($line);
#			my @ele = split (" ",$line);
			my @ele = split (" ",$line);
			foreach my $igene (@ele[1..$#ele]) {
				do {
					warn $_." UnIdentified.\n"
#					}unless ($igene =~/^(.+?)\|(.+?)$/);
					}unless ($igene =~/(.*?)\|(.*)/);
				#my $species = $1.".gff";
				my $species = $1;
				my $protein = $2;
				$ele[0] =~s/\://g;
				$ortholog{$species}{$protein} = $ele[0];
			}
		}
		close IN or die "$!";
	}
	return \%ortholog;
}

sub gff_load () {## gff{species}{genome id}{protein id} = (start,end,direction)
	my $ifile = shift;
	my $gff = shift;
	my $ortholog = shift;
	my $ispecies;
	if($ifile =~ /add/)
	{
		open (IN,$ifile) or die "$!";
		while (my $line = <IN>)
		{
			chomp($line);
			my @ele = split ("\t",$line);
			#$ispecies = $ele[1].'.gff';
			$ispecies = $ele[1];
			next unless (exists ${$ortholog}{$ispecies}{$ele[5]});
			@{${$gff}{$ispecies}{$ele[0]}{$ele[5]}} = ($ele[2],$ele[3],$ele[4]);
		}
		close IN or die "$!";
		return;
	}
	my $line;
	if ($ifile =~/\\/) {
		my @tmp = split/\\/,$ifile;		
		$ispecies = $tmp[-1];
	}elsif ($ifile =~/\//) {
		my @tmp = split/\//,$ifile;
		$ispecies = $tmp[-1];;
	}
	$ispecies =~ s/\.gff//;
	open (IN,$ifile) or die "$!";
	while ($line = <IN>) {
		chomp($line);
		next unless ($line =~/\tCDS\t/);
		next unless ($line =~/protein_id=(.+?);/);############
		my $id = $1;
		my @ele = split ("\t",$line);
		next unless ($line !~ /^\#/);
		next unless (exists ${$ortholog}{$ispecies}{$id});

		die "Unrecognized GTF or GTF format. start > end.\n" unless ($ele[3] < $ele[4]);
		@{${$gff}{$ispecies}{$ele[0]}{$id}} = ($ele[3],$ele[4],$ele[6]);
		#print $ispecies,"\n";sleep(1);
	}
	close IN or die "$!";
}

sub getDistanceMatrix () {
	my $self = shift;
	my $iformat = shift;
	my $model = shift || "SimiMatrixUn";
	if ($iformat eq "mega") {
		return $self->getDistanceMatrix_mega($model);
	}elsif ($iformat eq "phylip") {
		return $self->getDistanceMatrix_phy($model);
	}else{
		return $self->getDistanceMatrix_mega($model);
	}
}

sub getBootstrapMatrix()
{
	my $self = shift;
	my $value = shift;
	my $outputname = shift;
	for(my $i=0;$i<$value;$i++)
	{
		$self->refreshBootstrapMatrix(1);
		my @species = keys %{$self->{'BootstrapMatrix'}};
		open MAT,">>./Adjustment/tmp/tmp.distance";
		print MAT "sp\t",join("\t",@species),"\n";
		my $out = "";
		for my $sp(@species)
		{
			
			$out .= "$sp\t";
			for my $s(@species)
			{
				if($sp eq $s)
				{
					$out .= "0\t";
				}
				else
				{
					$out .= "$self->{BootstrapMatrix}{$sp}{$s}\t";
				}
			}
			$out =~ s/\t$//;
			print MAT $out;
			$out = '';
			print MAT "\n";
		}
		if(not -e "./Adjustment/tmp/distance.R")
		{
			open Rscript,">./Adjustment/tmp/distance.R";
			print Rscript "library(ape)\n";
			print Rscript 'mt<-read.table("./Adjustment/tmp/tmp.distance",head=T,sep="\t")',"\n";
			print Rscript 'mt<-as.matrix(mt)',"\n";
			print Rscript 'rownames(mt)<-mt[,1]',"\n";
			print Rscript 'mt<-mt[,-1]',"\n";
			print Rscript 'tr<-bionj(mt)',"\n";
			print Rscript 'write.tree(tr,file = "./Adjustment/tmp/tr.nwk")',"\n";
			close Rscript;
		}
		
		system("R --vanilla --slave <./Adjustment/tmp/distance.R");
		system("cat ./Adjustment/tmp/tr.nwk >> $outputname.nwk");
		system("rm ./Adjustment/tmp/tmp.distance ./Adjustment/tmp/tr.nwk");
	}
}

sub getDistanceMatrix_phy () {
	my $self = shift;
	my $model = shift || "SimiMatrixUn";
	my $ifile = "DistanceMatrix.Dis";
	die "Phylip format output was not available.\n";
}

sub getDistanceMatrix_mega () {
	my $self = shift;
	my $model = shift || "SimiMatrixUn";
	my $ifile = "DistanceMatrix.meg";
	$self->refreshSimiMatrix();
#	open (OUT,">".$ifile) or die "$!";
	my @species = keys %{$self->{$model}};
	my $out = "#mega
!TITLE  Generated by GTN Module;
!Format DataType=distance;
!Description
     GTN $model Method;
";
	$out .= "\n\#".join("\n\#",@species)."\n\n\n\n";
	for (my $i = 0; $i <= $#species; $i ++) {
		for (my $j = 0; $j <= $#species; $j++) {
			next unless ($j < $i);
			$out.=" ".$self->{$model}{$species[$i]}{$species[$j]};
		}
	$out .="\n";
	}
#	close OUT or die "$!";
	return $out;
}

sub gff_rank () {##gff sort by gene position network{species}{pre}{next}
	my $gff = shift;
	my %network = ();
	foreach my $ispecie (keys %{$gff}) {
		my %replicates = ();
		foreach my $icontig (keys %{${$gff}{$ispecie}}) {
			my %tmp = ();
			foreach my $iprotein (keys %{${$gff}{$ispecie}{$icontig}}) {
				if (not (exists $tmp{${${$gff}{$ispecie}{$icontig}{$iprotein}}[0]})) {
					$tmp{${${$gff}{$ispecie}{$icontig}{$iprotein}}[0]} = $iprotein;
				}else{
					warn $ispecie." ".$iprotein." ".$tmp{${${$gff}{$ispecie}{$icontig}{$iprotein}}[0]}." position redun.\n";
					$replicates{$iprotein}{$tmp{${${$gff}{$ispecie}{$icontig}{$iprotein}}[0]}} = 1;
					$replicates{$tmp{${${$gff}{$ispecie}{$icontig}{$iprotein}}[0]}}{$iprotein} = 1;
				}
			}
			my $temporal = "0";
			foreach my $irank (sort{$a <=> $b} keys %tmp) {
				if ($temporal eq "0") {
					$temporal = $tmp{$irank};
				}else{
					if ( (not (exists $replicates{$temporal})) && (not (exists $replicates{$tmp{$irank}})) ) {
						if (not (exists $network{$ispecie}{$temporal}{$tmp{$irank}})) {
							$network{$ispecie}{$temporal}{$tmp{$irank}} = 1;
							$network{$ispecie}{$tmp{$irank}}{$temporal} = 1;
						}else{
							$network{$ispecie}{$temporal}{$tmp{$irank}} ++;
							$network{$ispecie}{$tmp{$irank}}{$temporal} ++;
						}
					}else{
						if (exists $replicates{$temporal}) { ################
							foreach my $i (keys %{$replicates{$temporal}}) {
								if (not (exists $network{$ispecie}{$i}{$tmp{$irank}})) {
									$network{$ispecie}{$i}{$tmp{$irank}} = 1;
									$network{$ispecie}{$tmp{$irank}}{$i} = 1;
								}else{
									$network{$ispecie}{$i}{$tmp{$irank}} ++;
									$network{$ispecie}{$tmp{$irank}}{$i} ++;
								}
							}
						}elsif (exists $replicates{$tmp{$irank}}) {
							foreach my $i (keys %{$replicates{$tmp{$irank}}}) {
								if (not (exists $network{$ispecie}{$i}{$temporal})) {
									$network{$ispecie}{$i}{$temporal} = 1;
									$network{$ispecie}{$temporal}{$i} = 1;
								}else{
									$network{$ispecie}{$i}{$temporal} ++;
									$network{$ispecie}{$temporal}{$i} ++;
								}
							}
						}else{
							die "Graph::gff_rank erro.\n";
						}
					}
					$temporal = $tmp{$irank};
					#print $ispecies,"\n";
				}
			}
		}
	}
	return \%network;
}

sub Readingdir() {
	my $dir = shift;
	my $hit = shift;
	opendir(DIR, $dir) or die "$!";
		my @list = readdir(DIR);
		@list = grep {$_=~/($hit)/} @list;
	closedir DIR or die "$!";
	return \@list;
}

sub max {
	if ($_[0] > $_[1]) {
		return $_[0];
	}else{
		return $_[1];
	}
}

sub min {
	if ($_[0] > $_[1]) {
		return $_[1];
	}else{
		return $_[0];
	}
}

sub sum {
	my $isum = 0;
	map {$isum += $_} @_;
	return $isum;
}

sub mean {
	my $num = scalar(@_);
	return &sum(@_)/$num;
}

sub sd {
	my $mean = &mean(@_);
	my $sd = 0;
	for my $key (@_) {
		$sd += ($key-$mean)^2;
	}
	$sd = $sd/2;
	return sqrt($sd);
}


1;