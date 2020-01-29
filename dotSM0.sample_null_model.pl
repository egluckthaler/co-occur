#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use File::Basename;
use Data::Dumper;
use FileHandle;
use Sort::Naturally;
use Math::Random::Secure qw(irand);
use JSON::XS qw(encode_json decode_json);
use File::Slurp qw(read_file write_file);
use Ebinf::Utils qw(Coo3_parse Dependency_check dim_1_hash dim_0_hash);
use Ebinf::Cluster_utils qw(Clusterfy);
$| = 1;

#ome is defined by everything up to the last _
sub usage {
	my $message = shift;
	$message = '' if (not defined $message);
	my $usage = qq/
usage: $0 [OPTIONS]
This script randomly samples genomes of interest for gene pairs with 1-6 intervening
genes and calculates some stats about the frequency of their pairing (by MCL) in all genomes. 
If previous output is provided, then MCL pairs that have already been sampled are taken into
account and are not resampled\n\n/;
	die($usage, $message);
}
 
main: {
	#my %opts;
	#getopt('i:m:o', \%opts);
	#Opts_check(\%opts);
	my $bigcoo3file = "/fs/project/PAS1046/projects/dotSM/coo3/all_dothideo.coo3";
	my $bigprotidfile = "/fs/project/PAS1046/projects/dotSM/ProtDB/all_dothideo_prots.ids";
	my $rawmclfile = "/fs/project/PAS1046/projects/dotSM/OrthoMCL/orthoMCL.output";
	my $omesoifile = "/fs/project/PAS1046/projects/dotSM/Metadata/ome2order";
	Dependency_check($bigcoo3file, $bigprotidfile, $rawmclfile, $omesoifile);
	my $sampledmclfile = "/fs/project/PAS1046/projects/dotSM/Results/null_distribution_cooccurrences.txt";
	my $maxiterations = 500000;
	open(my $nullout, '>>', $sampledmclfile) or usage("Error: cannot open $sampledmclfile for writing\n");
	$nullout->autoflush(1);
#1. Grab all omes of interest
	my ($omesoi) = dim_1_hash($omesoifile, "\t", "0:1");

#2. Grab all proteins from all dothideomes
	my ($prothash) = Hash_all_prots($bigprotidfile); #temporary sub to conform with coo3_parse specifications
	my @protdb = keys %{$prothash};
	
#3. Store all 7 prots upstream and 7 prots downstream of all proteins. Corresponds to a max of 6 intervening genes. structured{prot}{'up'} = [..]
#	Takes approx 3 min to run
	my $prothoods;
	my $prothoodjsonfile = "/fs/project/PAS1046/projects/dotSM/Results/json/prothoods.7.json";
	system("mkdir /fs/project/PAS1046/projects/dotSM/Results/json/") if (! -d "/fs/project/PAS1046/projects/dotSM/Results/json/");
	if (! -f $prothoodjsonfile) {
		print "Precomputing 7 prots upstream and 7 prots downstream of all proteins..\n";
		my ($bigcoo3) = Coo3_parse($bigcoo3file, $prothash);	
		($prothoods) = Prot_hoods($bigcoo3, 7);
		my $prothoodjson = encode_json($prothoods); write_file($prothoodjsonfile, { binmode => ':raw' }, $prothoodjson);
	} else {
		print "Loading previous data of prothoods from json file..\n";
		my $prothoodjson = read_file($prothoodjsonfile, { binmode => ':raw' }); $prothoods = decode_json($prothoodjson);
	}
	
#3. Read previous sampled MCLs into mem, if provided
#	sampledmcl is structured {mcl1-mcl2} = 1, where mcl1 always comes before mcl2 by nsort
	my $sampledmcl = {};
	my $iteration = 0;
	if (-f $sampledmclfile) {
		($sampledmcl) = Read_sampled_MCL($sampledmclfile) if (defined $sampledmclfile);
		$iteration = scalar keys %{$sampledmcl}; #keeps track of how many non-redundant MCL pairs have been sampled
		print "Reading in previously sampled MCL cooccurrence file.. starting at iteration $iteration\n";
	}
#4. Read ortholog groups file into a hash where {protid => OG} and another where {OG => protid => 1}
#	name clusters starting at MCL000001
	my ($prot2og, $og2prot) = OMCL_raw_parser($rawmclfile);

#5. Build null model
	ITER: while ($iteration <= $maxiterations) {
		#a. Randomly pick a MCL-MCL pair separated by 1-6 intervening genes that has not been observed before
		my ($mclpair) = Pick_new_MCL_pair($prothoods, $prot2og, \@protdb, $sampledmcl);
		goto ITER if (exists $sampledmcl->{$mclpair}); #ignore instances of mcls that have been sampled already
		
		#b. Retrieve all prots associated with each mcl
		my ($mcl1, $mcl2) = split/-/, $mclpair;
		my %protids1 = map {$_ => 1} keys %{$og2prot->{$mcl1}};
		my %protids2 = map {$_ => 1} keys %{$og2prot->{$mcl2}};
		
		#c. Find all instances of minimal clustering of mcl1 and mcl2, keeping track of which protids have already been counted as a co-occurring pair
		my ($cooccurrences, $unknown1, $unknown2) = (0, 0, 0);
		my $mcl1count = scalar keys %protids1;
		my $mcl2count = scalar keys %protids2;
		goto ITER if ($mcl1count == 1 || $mcl2count == 1); #ignore instances where mcl is a singleton
		my %countedprots;
		foreach my $protid1 (keys %protids1) { #iterate through all prots assigned to MCL1
			#neighbor structure {abs_dist_away}{protid} = 1;
			my ($neighbors, $nearend) = Retrieve_all_neighbor_prots($protid1, $prothoods); 
			my $cooccurfound = 0;
			COOCCUR: foreach my $dist (nsort keys %{$neighbors}) {#pick the closest available match
				foreach my $prot (keys %{$neighbors->{$dist}}) {
					next if (exists $countedprots{$prot}); #don't double count proteins that have already been counted as part of a cooccurrence
					if (exists $protids2{$prot}) { #is this prot assigned to MCL2?
						$cooccurrences++;
						$countedprots{$prot} = 1; #no need to keep track of protid1, since it will only ever be checked once for clustering with prots from MCL2
						$cooccurfound = 1; 
						last COOCCUR;
					}
				}
			}
			#if all dists away from protid1 have been examined, but no cooccurrence found, check if protid1 is within 7 genes of the end of the contig
			if ($cooccurfound == 0) {
				$unknown1++ if ($nearend == 1); #is this prot within 7 genes of end of a contig?
			}
		}
		#now, check those protid2s that were not found to be co-occurring with protid1 to see if they fall within contig end
		foreach my $protid2 (keys %protids2) {
			next if (exists $countedprots{$protid2}); #don't look at protids that are part of co-occurrences
			my ($neighbors, $nearend) = Retrieve_all_neighbor_prots($protid2, $prothoods);
			$unknown2++ if ($nearend == 1); #is this prot within 7 genes of end of a contig?
		}
		#assess frequency of cooccurence of 1->2 and 2->1, discounting unknowns
		my ($freq1, $freq2) = ("n.d.", "n.d.");
		$freq1 = sprintf("%.3f", $cooccurrences/($mcl1count-$unknown1)) if (($mcl1count-$unknown1) != 0);
		$freq2 = sprintf("%.3f", $cooccurrences/($mcl2count-$unknown2)) if (($mcl2count-$unknown2) != 0);
		print $nullout "$mcl1\t$mcl2\t$cooccurrences\t$mcl1count\t$mcl2count\t$unknown1\t$unknown2\t$freq1\t$freq2\n";
		$sampledmcl->{$mclpair} = 1;
		$iteration++;
	}
}

sub Retrieve_all_neighbor_prots {
	my ($query, $prothoods) = @_;
	my (%neighbors);
	my $nearend = 0;
	foreach my $direction (keys %{$prothoods->{$query}}) {
		for(my $i = 0; $i < scalar @{$prothoods->{$query}->{$direction}}; $i++) {
			$neighbors{$i}{$prothoods->{$query}->{$direction}[$i]} = 1; #where i is the abs distance from query
		}
		#test to see if protein is near the end of a contig; this will happen when there are not 7 prots either upstream or downstream of query
		if (scalar @{$prothoods->{$query}->{$direction}} < 7) {
			$nearend = 1;
		}
	}
	return(\%neighbors, $nearend);
}


sub Pick_new_MCL_pair {
	usage("Error: incorrect number of args to Pick_new_MCL_pair\n") if (scalar @_ != 4);
	my ($prothoods, $prot2og, $protdb, $sampledmcl) = @_;
	my $mclpair;
	while (not defined $mclpair) {
		#a. First, randomly pick a pair of proteins separated by 1-6 genes from set of all proteins. Keep picking new pairs if no pair is found
		my ($ranprot, $neighborprot) ;
		while (not defined $neighborprot) {
			$ranprot = Pick_random_element($protdb);
			#b. Randomly pick a distance between 0-6 (represents up to 6 intervening genes)
			my $randist = irand(7);
			#print "checking for neighbors within $randist of $ranprot\n";
			#c. Retrieve protein within random distance of random protein. Try upstream first, then downstream
			if ((exists $prothoods->{$ranprot}->{'up'}) && ($randist + 1 <= scalar @{$prothoods->{$ranprot}->{'up'}})) {
				$neighborprot = $prothoods->{$ranprot}->{'up'}[$randist];
				#print "found $neighborprot\n";
			} elsif ((exists $prothoods->{$ranprot}->{'down'}) && ($randist + 1 <= scalar @{$prothoods->{$ranprot}->{'down'}})) {
				$neighborprot = $prothoods->{$ranprot}->{'down'}[$randist];
				#print "found $neighborprot\n";
			}
		}
		#b. Second, check that the mcl of the protpair has not been observed before
		my @mcls = ($prot2og->{$ranprot}, $prot2og->{$neighborprot});
		my $mclputative = join("-", nsort @mcls);
		if ($prot2og->{$ranprot} ne $prot2og->{$neighborprot}) {
			if (not exists $sampledmcl->{$mclputative}) {
				$mclpair = $mclputative;
			}
		}			
	}
	return($mclpair);
}



sub Hash_all_prots {
	my ($protfile) = @_;
	open(my $protin, '<', $protfile);
	my %protids;
	while (my $line = <$protin>) {
		chomp $line;
		push @{$protids{$line}}, "null";
	}
	return(\%protids);
}

sub Prot_hoods {
	my $usage = qq/
(<\%prothoods>) = Prot_hoods(<big_coo3_hash>, <upstream+downstream_max>)
sub Prot_hood grabs all prots that fall within max upstream and downstream
of each prot in a coo3 hash. If a prot is within max of end of contig, then the array of
upstream or downstream prots will not be max size\n/;
	die($usage, "Error: incorrect number of args\n") if scalar (@_ != 2);
	my ($bigcoo3, $maxprotaway) = @_;
	my %prothoods;
	foreach my $ome (keys %{$bigcoo3}) {
		foreach my $contig (keys %{$bigcoo3->{$ome}}) {
			my %contigpositions = map {$_ => 1} sort keys %{$bigcoo3->{$ome}->{$contig}}; #copy all positions in the contig, to protect against autovivificaiton
			foreach my $position (sort keys %{$bigcoo3->{$ome}->{$contig}}) {
				foreach my $prot (sort keys %{$bigcoo3->{$ome}->{$contig}->{$position}}) {
					#find the range of possible upstream positions, then check if proteins exist at those positions
					my ($upstream_putative, $downstream_putative) = Find_putative_positions($position, $maxprotaway);
					my ($upstream_true, $downstream_true) = Find_true_positions($upstream_putative, $downstream_putative, \%contigpositions);
					#Extract all proteins at true positions, sorted from closest to farthest away
					foreach my $trueup (sort @{$upstream_true}) {
						my @prot = keys %{$bigcoo3->{$ome}->{$contig}->{$trueup}}; #there will only ever be 1 prot per position
						push @{$prothoods{$prot}{'up'}}, @prot;
					}
					foreach my $truedown (sort @{$downstream_true}) {
						my @prot = keys %{$bigcoo3->{$ome}->{$contig}->{$truedown}}; #there will only ever be 1 prot per position
						push @{$prothoods{$prot}{'down'}}, @prot;
					}
				}
			}
		}
	}
	return(\%prothoods);
}


sub Find_true_positions {
	my ($upstream_putative, $downstream_putative, $real_positions) = @_;
	my (@upstream_true, @downstream_true);
	foreach my $putative (sort @{$upstream_putative}) {
		if (exists $real_positions->{$putative}) {
			push @upstream_true, $putative;
		}
	}
	foreach my $putative (sort @{$downstream_putative}) {
		if (exists $real_positions->{$putative}) {
			push @downstream_true, $putative;
		}
	}
	return(\@upstream_true, \@downstream_true);
}


sub Find_putative_positions {
	my ($currentposition, $maxdist) = @_;
	my (@dists) = (1..$maxdist);
	my (@uppositions, @downpositions);
	foreach my $dist (sort @dists) {
		push @uppositions, $currentposition + $dist;
		push @downpositions, $currentposition - $dist;
		
	}
	return(\@uppositions, \@downpositions);
}

sub Pick_random_element {
	my ($dbref) = (@_);
	my @db = @{$dbref};
	my $elementnum = scalar @db;
	my $int = irand($elementnum); #generate a random integer between 0 and the total number of omes - 1
	my $randelement = $db[$int];
	return($randelement);
}



sub OMCL_raw_parser {
	my $usage = qq/
my (<\%prot2og) = OMCL_raw_parser(<path_to_OMCL_output>)
sub OMCL_raw_parser reads in a raw OMCL output file and returns a hash where 
{proteinID => OG}. OG names start at MCL000001, and increment as the file is 
processed line by line.\n\n/;
	die($usage, "Error: incorrect number of args\n") if (scalar @_ != 1);
	my ($omclfile) = @_;
	open(my $omclin, '<', $omclfile) or die($usage, "Error: cannot open $omclfile\n");
	my (%prot2og, %og2prot, $mclcount);
	while(my $line = <$omclin>) {
		$mclcount++;
		my ($formattedOG) = sprintf("%06d", $mclcount);
		$formattedOG = "MCL$formattedOG";
		chomp $line;
		my (@protids) = split/\s/, $line;
		foreach my $protid (@protids) {
			$protid =~ s/\|/_/;
			$prot2og{$protid} = $formattedOG;
			$og2prot{$formattedOG}{$protid} = 1;
		}
	}
	return(\%prot2og, \%og2prot);
}
	

sub Read_sampled_MCL {
	my ($sampledmclfile) = @_;
	my %sampledmcl;
	open(my $sampledin, '<', $sampledmclfile) or usage("Error: cannot open $sampledmclfile for reading\n");
	while (my $line = <$sampledin>) {
		chomp $line;
		my ($mcl1, $mcl2) = split/\t/, $line;
		my @mcls = ($mcl1, $mcl2);
		@mcls = nsort @mcls;
		$sampledmcl{"${mcl1}-${mcl2}"} = 1;
	}
	return(\%sampledmcl);
}

