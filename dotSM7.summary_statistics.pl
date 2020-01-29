#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use File::Basename;
use Data::Dumper;
use FileHandle;
use Sort::Naturally;
use Statistics::Descriptive;
use Ebinf::Utils qw(Open_FH Dependency_check dim_1_hash);
$| = 1; 
sub usage {
	my $message = shift;
	$message = '' if (not defined $message);
	my $usage = qq/
usage: $0 [OPTIONS]

-g: <cluster_group_file>
-s: <cluster_singletons_file>
-m: <minimum_cluster_size>
-o: <output_directory>

This script parses main outputs of dotSM analysis to produce files for
data visualization in R\n\n/;
	die($usage, $message);
}

main: {
	my %opts;
	getopt('g:s:m:o', \%opts);
	Opts_check(\%opts);
	my ($groupsfile, $singletonsfile, $MINSIZE, $outdir) = ($opts{'g'}, $opts{'s'}, $opts{'m'}, $opts{'o'});
	my $header = "#output by dotSM5.extract_summary_statistics.pl using $groupsfile and $singletonsfile\n";
	my $countbysizefile = "$outdir/dotSM5.groupcount_by_groupsize.mincl$MINSIZE.txt";
	my ($graph1out) = Open_FH($countbysizefile, $header);
	my $memberbygroupfile = "$outdir/dotSM5.membercount_by_groupsize.mincl$MINSIZE.txt";
	my ($graph2out) = Open_FH($memberbygroupfile, $header);
	my $countbySMclassfile = "$outdir/dotSM5.clustercount_by_SMclass.mincl$MINSIZE.txt";
	my ($graph3out) = Open_FH($countbySMclassfile, $header);
	my $clustercountfile = "$outdir/dotSM5.SMclusters_per_ome.matrix";
	my ($graph4out) = Open_FH($clustercountfile, $header);
	my $clustersizefile = "$outdir/dotSM5.SMcluster_size_per_ome.table";
	my ($graph5out) = Open_FH($clustersizefile, $header);
	
	#this hardcoded file contains all mcl pairs with significant cooccurrences
	my $mclfile = "/fs/project/PAS1046/projects/dotSM/Results/nonsigSM_cluster_cooc_500001/dotSM3.all_cooccurrence_pvalues.txt";
	my $ome2orderfile = "/fs/project/PAS1046/projects/dotSM/Metadata/ome2order";
	Dependency_check($mclfile, $ome2orderfile);
	my $ome2order = dim_1_hash($ome2orderfile, "\t", "0:1");

	
	#1. Graph1: Number of homologous cluster groups vs. group size, grouped by clustergroup and singletons
	my ($groupstats, $singlestats) = Parse_cluster_file_graph1($groupsfile, $singletonsfile, $MINSIZE);
	Print_graph1($groupstats, $singlestats, $graph1out);
	
	#2. Graph2: Number of members per clustergroup, by groupsize
	($groupstats) = Parse_cluster_file_graph2($groupsfile, $MINSIZE);
	Print_graph2($groupstats, $graph2out);

	#3. Graph3: Number of clustergroups containing a PKS etc.., grouped by clustergroup and singletons
	($groupstats, $singlestats) = Parse_cluster_file_graph3($groupsfile, $singletonsfile, $MINSIZE);
	Print_graph3($groupstats, $singlestats, $graph3out); 
	
	#4. Graph4: total number of clusters per genome: SMURF and newSMURF
	#retrieve total number of predicted SMURF clusters per genome
	my $clusterdir = "/fs/project/PAS1046/projects/dotSM/SecondaryMetabolism";
	my @files = glob("$clusterdir/*.smclusters");
	my ($ome2smurfcount, $ome2SMURFsizes);
	foreach my $clusterfile (@files) {
		($ome2smurfcount, $ome2SMURFsizes) = SMURF_parser($clusterfile, $ome2smurfcount, $ome2SMURFsizes); #update ome2smurf with each iteration
	}
	#retrieve total number of newSMURF clusters per genome, only counting clusters with signature genes, and separately counting clusters without
	my ($ome2newcount, $ome2ndcount);
	($ome2newcount, $ome2ndcount) = Count_new_clusters($groupsfile, $ome2newcount, $ome2ndcount);
	($ome2newcount, $ome2ndcount) = Count_new_clusters($singletonsfile, $ome2newcount, $ome2ndcount);
	#print out counts 
	Print_cluster_counts($ome2smurfcount, $ome2newcount, $ome2ndcount, $graph4out);
	
	#5. Graph5: distribution of SMURF and newSMURF cluster sizes (in unique MCL groups) per genome, with order information
	#structure @{$ome} = [sizes]
	my ($ome2newSMURFsizes) = Parse_cluster_file_graph5($groupsfile, $singletonsfile);
	Print_graph5($ome2SMURFsizes, $ome2newSMURFsizes, $ome2order, $graph5out);
	
}

sub Print_graph5 {
	my ($ome2SMURFsizes, $ome2newSMURFsizes, $ome2order, $out) = @_;
	print $out "ome\torder\tcluster_type\tsize\n";
	foreach my $ome (nsort keys %{$ome2SMURFsizes}) {
		foreach my $size (@{$ome2SMURFsizes->{$ome}}) {
			print $out "$ome\t$ome2order->{$ome}\tSMURF\t$size\n";
		}
		foreach my $size (@{$ome2newSMURFsizes->{$ome}}) {
			print $out "$ome\t$ome2order->{$ome}\tnewSMURF\t$size\n";
		}
	}
}



sub Parse_cluster_file_graph5 {
	my ($groupsfile, $singletonsfile) = @_;
	my (%newSMURFcounts);
	open(my $groupin, '<', $groupsfile) or usage("Error: could not open $groupsfile for reading\n");
	while (my $line = <$groupin>) {
		next if ($line =~ m/^#/);
		chomp $line;
		my ($groupname, $clustername, $SMstring, $MCLstring) = split/\t/, $line;
		$clustername =~ m/^(.+?)_[^_]+$/; #ome is defined by up to last _
		my $ome = $1;
		my @MCLs = split/,/, $MCLstring;
		my %MCLs = map {$_ => 1} @MCLs;
		my $mclcount = scalar keys %MCLs;		
		push @{$newSMURFcounts{$ome}}, $mclcount;
	}
	open(my $singlein, '<', $singletonsfile) or usage("Error: could not open $singletonsfile for reading\n");
	while (my $line = <$singlein>) {
		next if ($line =~ m/^#/);
		chomp $line;
		my ($groupname, $clustername, $SMstring, $MCLstring) = split/\t/, $line;
		$clustername =~ m/^(.+?)_[^_]+$/; #ome is defined by up to last _
		my $ome = $1;
		my @MCLs = split/,/, $MCLstring;
		my %MCLs = map {$_ => 1} @MCLs;
		my $mclcount = scalar keys %MCLs;		
		push @{$newSMURFcounts{$ome}}, $mclcount;
	}
	return(\%newSMURFcounts);
}

sub Print_cluster_counts {
	my ($ome2smurfcount, $ome2newcount, $ome2ndcount, $out) = @_;
	print $out "ome\tSMURF\tnewSMURF\tn.d.\n";
	foreach my $ome (nsort keys %{$ome2smurfcount}) {
		my $ndcount = 0;
		$ndcount = $ome2ndcount->{$ome} if (exists $ome2ndcount->{$ome});
		print $out "$ome\t$ome2smurfcount->{$ome}\t$ome2newcount->{$ome}\t$ndcount\n";
	}
}


sub Count_new_clusters {
	my ($groupsfile, $ome2newcount, $ome2ndcount) = @_;
	open(my $groupin, '<', $groupsfile) or usage("Error: could not open $groupsfile for reading\n");
	while (my $line = <$groupin>) {
		next if ($line =~ m/^#/);
		chomp $line;
		my ($groupname, $omecluster, $SMclass) = split/\t/, $line;
		$omecluster =~ m/^(.+)_[^_]+$/; #capture everything up to last _
		my $ome = $1;
		if ($SMclass =~ m/n\.d\./) {
			$ome2ndcount->{$ome}++;
		} else {
			$ome2newcount->{$ome}++;
		}
	}
	return($ome2newcount, $ome2ndcount);
}

sub SMURF_parser {
	my ($clusterfile, $ome2smurf, $ome2smurfsize) = @_;
	my ($omecode) = fileparse($clusterfile, ".smclusters");
	open(my $clusterin, '<', $clusterfile) or usage("Error: cannot open $clusterfile\n");
	while (my $line = <$clusterin>) {
		next if ($line =~ m/^#/);
		chomp $line;
		if ($line =~ m/Cluster_size: (\d+)/) {
			my $clustersize = $1;
			$ome2smurf->{$omecode}++ if ($clustersize > 1);
			push @{$ome2smurfsize->{$omecode}}, $clustersize;
		}
	}
	return($ome2smurf, $ome2smurfsize);
}

sub Print_graph3 {
	my ($groupstats, $singlestats, $out) = @_;
	print $out "SMclass\tnumber_of_clustergroups\ttype\n";
	foreach my $SMtype (nsort keys %{$groupstats}) {
		print $out "$SMtype\t".scalar(keys(%{$groupstats->{$SMtype}}))."\tclustergroup\n";
	}
	foreach my $SMtype (nsort keys %{$singlestats}) {
		print $out "$SMtype\t".scalar(keys(%{$singlestats->{$SMtype}}))."\tsingleton\n";
	}
}

sub Parse_cluster_file_graph3 {
	my ($groupsfile, $singletonsfile, $MINSIZE) = @_;
	my (%groupstats, %singlestats);
	open(my $groupin, '<', $groupsfile) or usage("Error: could not open $groupsfile for reading\n");
	while (my $line = <$groupin>) {
		next if ($line =~ m/^#/);
		chomp $line;
		my ($groupname, $cluster, $SMclass) = split/\t/, $line;
		my ($numbers, $size) = split/_/, $groupname;
		$size =~ s/s//;
		$groupstats{$SMclass}{$groupname} = 1 if ($size >= $MINSIZE);
	}
	open(my $singlein, '<', $singletonsfile) or usage("Error: could not open $singletonsfile for reading\n");
	while (my $line = <$singlein>) {
		next if ($line =~ m/^#/);
		chomp $line;
		my ($groupname, $cluster, $SMclass) = split/\t/, $line;
		my ($numbers, $size) = split/_/, $groupname;
		$size =~ s/s//;
		$singlestats{$SMclass}{$groupname} = 1 if ($size >= $MINSIZE);
	}
	return(\%groupstats, \%singlestats);		
}
sub Print_graph2 {
	my ($memberstats, $graph2out) = @_;
	print $graph2out "clustergroup_size\tnumber_of_members\n";
	foreach my $size (nsort keys %{$memberstats}) {
		foreach my $membersize (@{$memberstats->{$size}}) {
			print $graph2out "$size\t$membersize\n";
		}
	}
}

sub Parse_cluster_file_graph2 {
	my ($groupsfile, $MINSIZE) = @_;
	my (%groupstats);
	open(my $groupin, '<', $groupsfile) or usage("Error: could not open $groupsfile for reading\n");
	my %group2membercount;
	while (my $line = <$groupin>) {
		next if ($line =~ m/^#/);
		chomp $line;
		my ($groupname) = split/\t/, $line;
		my ($numbers, $size) = split/_/, $groupname;
		$size =~ s/s//;
		$group2membercount{$groupname}++ if ($size >= $MINSIZE);
	}
	foreach my $group (keys %group2membercount) {
		my ($groupname, $size) = split/_/, $group;
		$size =~ s/s//;
		push @{$groupstats{$size}}, $group2membercount{$group};
	}
	return(\%groupstats);
}
		
sub Print_graph1 {
	my ($groupstats, $singlestats, $out) = @_;
	print $out "clustergroup_size\tclustergroup_count\ttype\n";
	foreach my $size (nsort keys %{$groupstats}) {
		print $out "$size\t".scalar(keys(%{$groupstats->{$size}}))."\tclustergroup\n";
	}
	foreach my $size (nsort keys %{$singlestats}) {
		print $out "$size\t".scalar(keys(%{$singlestats->{$size}}))."\tsingleton\n";
	}
}

sub Parse_cluster_file_graph1 {
	my ($groupsfile, $singletonsfile, $MINSIZE) = @_;
	my (%groupstats, %singlestats);
	open(my $groupin, '<', $groupsfile) or usage("Error: could not open $groupsfile for reading\n");
	while (my $line = <$groupin>) {
		next if ($line =~ m/^#/);
		chomp $line;
		my ($groupname) = split/\t/, $line;
		my ($numbers, $size) = split/_/, $groupname;
		$size =~ s/s//;
		$groupstats{$size}{$groupname} = 1 if ($size >= $MINSIZE);
	}
	open(my $singlein, '<', $singletonsfile) or usage("Error: could not open $singletonsfile for reading\n");
	while (my $line = <$singlein>) {
		next if ($line =~ m/^#/);
		chomp $line;
		my ($groupname) = split/\t/, $line;
		my ($numbers, $size) = split/_/, $groupname;
		$size =~ s/s//;
		$singlestats{$size}{$groupname} = 1 if ($size >= $MINSIZE);
	}
	return(\%groupstats, \%singlestats);
}


sub Opts_check {
	my ($opts) = @_;
	usage("Error: please provide a cluster groups file\n") if (not defined $opts->{'g'});
	usage("Error: please provide a cluster singleton file\n") if (not defined $opts->{'s'});
	usage("Error: please provide an output directory\n") if (not defined $opts->{'o'});
	usage("Error: please provide a minimum cluster size\n") if (not defined $opts->{'m'});
	usage("Error: cannot find the cluster groups file\n") if (! -f $opts->{'g'});
	usage("Error: cannot find the cluster singleton file\n") if (! -f $opts->{'s'});
	usage("Error: cannot find the output directory\n") if (! -d $opts->{'o'});
}
