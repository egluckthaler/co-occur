#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use File::Basename;
use Data::Dumper;
use FileHandle;
use Sort::Naturally;
use Ebinf::Utils qw(dim_1_hash Dependency_check Open_FH);
use JSON::XS qw(encode_json decode_json);
use File::Slurp qw(read_file write_file);
$| = 1;

#ome is defined by everything up to the last _
sub usage {
	my $message = shift;
	$message = '' if (not defined $message);
	my $usage = qq/
usage: $0 [OPTIONS]
-o <output_directory>
This script retrieves all proteins associated with MCL groups from 
SMURF cluster predictions (irregular format supplied by JGI; can be customized in future 
to accept other formats). The set of MCLs associated with this set of proteins is retrieved
and then all unique MCL cooccurrences are counted. The frequency of cooccurrence is then
be compared against the null model of MCL cooccurrences (created by dotSM.sample_null_model.pl)
Any pairs with 5 or less cooccurrences are not deemed significant.
The frequency of each cooccurrence in newSMURF clusters (produced by dotSM1.cluster_all_SM_prots.pl)
is also computed. \n\n/;
	die($usage, $message);
}
 
main: {

	my %opts;
	getopt('o:h', \%opts);
	Opts_check(\%opts);
	my ($outdir) = ($opts{'o'});
	my $clusterdir = "/fs/project/PAS1046/projects/dotSM/SecondaryMetabolism";
	my $orthofile = "/fs/project/PAS1046/projects/dotSM/OrthoMCL/orthoMCL.output";
	my $prothoodjsonfile = "/fs/project/PAS1046/projects/dotSM/Results/json/prothoods.7.json"; #output by dotSM0.sample_null_model.pl
	my $ome2orderfile = "/fs/project/PAS1046/projects/dotSM/Metadata/ome2order";
	my $nullmodelfile = "/fs/project/PAS1046/projects/dotSM/Results/null_distribution_cooccurrences.500001.txt"; #output by dotSM0.sample_null_model.pl
	my $mclannfile = "/fs/project/PAS1046/projects/dotSM/Results/dotSM_OG_annotation.txt"; #manually output by Emapper analysis
	my $combinedfile = "/fs/project/PAS1046/projects/dotSM/Results/dotSM_combined.txt"; #output by dotSM3.<SM_algorithm>2combined.pl
	Dependency_check($ome2orderfile, $prothoodjsonfile, $nullmodelfile, $mclannfile, $orthofile, $combinedfile);
	
	#usage("${analysisname}_all_SM_cooccurrences.txt already exists in $outdir\n") if (-f "$outdir/${analysisname}_all_SM_cooccurrences.txt");
	my $header = "#output by dotSM3.test_all_SM_cooccurrences.pl\n";
	my $pvaluefile = "$outdir/dotSM3.all_cooccurrence_pvalues.txt";
	my ($pvalueout) = Open_FH($pvaluefile, $header);
	my $allpvaluefile = "$outdir/dotSM3.truly_all_cooccurrence_pvalues.txt";
	my ($allpvalueout) = Open_FH($allpvaluefile, $header);
	my $sigpvaluefile = "$outdir/dotSM3.all_cooccurrence_pvalues_signetwork.txt";
	my ($sigpvalueout) = Open_FH($sigpvaluefile, $header);
	my $cytofile = "$outdir/dotSM3.all_cooccurrence_pvalues_significant.network.txt";
	my ($cytoout) = Open_FH($cytofile, $header);
	my $smclassfile = "$outdir/dotSM3.all_cooccurrence_pvalues_smclass.txt";
	my ($smout) = Open_FH($smclassfile, $header);
	my $sigmcloutfile = "$outdir/dotSM3.all_cooccurrence_pvalues_significantMCLs.txt";
	my ($sigogout) = Open_FH($sigmcloutfile, $header);

	#grab list of omes of interest
	my ($omesoi) = dim_1_hash($ome2orderfile, "\t", "0:1");

#1. Read ortholog groups file into a hash where {protid => OG} and {OG=> protid => 1}
#	name clusters starting at MCL000001
	my ($prot2og, $og2prot) = OMCL_raw_parser($orthofile);

#5. For each non-redundant MCL pair separated by 6 genes or less, find all cooccurrences separated by no more than 6 intervening genes. Print out file for recordkeeping and return each line in file
#	structure: {mcl1-mcl2} = line from coocurrence file
#	has to be submitted as job to run.. take ~ an hour
	my $allmclpairs;
	my $cooccurrencefile = "/fs/project/PAS1046/projects/dotSM/Results/dotSM3.all_SM_cooccurrences.sorted.txt";
	if (! -f $cooccurrencefile) {
		#	Read in combined file with all new SMURF bootstrapped cluster predictions
		#	structured $cluster2prot->{$clusterid}->{$protid} = clusterclass
		#	spit out another hash: prot2clusterclass, $prot2cluster
		#	'unknown' as SM class are not stored
		my ($prot2clusterclass, $prot2cluster, $cluster2prot, $clusterogs) = Combined_parser($combinedfile, $omesoi);
		#	Store all 7 prots upstream and 7 prots downstream of all proteins. Corresponds to a max of 6 intervening genes. structured{prot}{'up'} = [..]
		#	Takes approx 3 min to run
		my $prothoods;
		print "Loading previous data of prothoods from json file..\n";
		my $prothoodjson = read_file($prothoodjsonfile, { binmode => ':raw' }); $prothoods = decode_json($prothoodjson);
		#	Iterate through all proteins in all newSMURF clusters and find set of MCL cooccurrences that are separated by no more than 6 intervening genes
		#	structure: {mcl1-mcl2} = 1
		print "Finding set of MCL cooccurrences that are separated by no more than 6 intervening genes in at least 1 genome\n";
		my ($sampledmclpairs) = Find_all_MCL_pairs($cluster2prot, $prot2og, $prothoods);
		open(my $coout, '>', $cooccurrencefile) or usage("Error: cannot open $outdir/dotSM3_all_SM_cooccurrences.sorted.txt for printing\n");
		$coout->autoflush(1);
		print "Finding all occurrences of MCL cooccurrences and printing them out\n";
		($allmclpairs) = Find_and_print_all_cooccurrences($sampledmclpairs, $og2prot, $prothoods, $prot2cluster, $prot2clusterclass, $coout);
	} else {
		print "$cooccurrencefile already exists.. reading it in..\n";
		($allmclpairs) = Cooccurrence_parser($cooccurrencefile);
	}
		
#6. Read in null model. For each OG size range in increments of 25, gather # of cooccurrences
#	structured {smaller OG}{larger OG} = [distribution of cooccurrences]
	print "Determining null model distributions based on pairwise MCL family size in increments of 25..\n";
	my ($nullmodel) = Null_model_distributions($nullmodelfile);

#7. Test each MCL pair for significance vs. the null model
#	Significance is determined by (# cooccurrences | sizes of MCL groups in pair)
#	Return all mclpairs where p<0.05 

#	Retrieve the most frequent MCL annotation
	my ($mcl2ann) = Retrieve_MCL_annotation($mclannfile);

#	Retrieve all MCLs associated with signature SM genes, as identified by SMURF
#	structured {signatureMCLS}->{signatureclass} = count
	print "Parsing signature gene data from SMURF predictions in $clusterdir\n";
	my @files = glob("$clusterdir/*.smclusters");
	my ($sigogs);
	foreach my $clusterfile (@files) {
		($sigogs) = SMURF_parser($clusterfile, $prot2og, $sigogs); #update sigogs with each iteration
	}
	my ($sig2ann) = Retrieve_signature_annotation($sigogs);

	print "Determining significance of all mcl pairs by comparing against cooccurrence distributions in $nullmodelfile\n";
	
	#unhash next section to assess co-occurrences occurring in entire set (global)
	#print "Estimating pvalue using cooccurrences that occur globally\n"; 
	#my ($sigpairs, $sigmcls, $sigSMcount, $sigpairclustercount, $sigcooccs) = Significance_testing($allmclpairs, $nullmodel, $pvalueout);
	
	#unhash next section to assess co-occurrences occurring in predicted newSMURF clusters only
	print "Estimating pvalue using cooccurrences that occur only in predicted clusters\n";
	my ($sigpairs, $sigmcls, $sigSMcount, $sigpairclustercount, $sigcooccs) = Significance_testing($allmclpairs, $nullmodel, $pvalueout, $allpvalueout);
	
	#unhash next section to assess co-occurrences occurring in predicted newSMURF clusters only AND that are found in a network of significant co-occurrences with a signature SM gene
	#print "Estimating pvalue using cooccurrences that occur only in predicted clusters\n";
	#my ($sigpairs, $sigmcls, $sigSMcount, $sigpairclustercount, $sigcooccs) = Significance_testing($allmclpairs, $nullmodel, $pvalueout);
	#print "Extracting unexpected mclpairs that occur in a network with a signature gene..\n";
	#my ($signetworkpairs) = Check_network_for_signature($sigpairs, $sig2ann); #Check_network_for_signature works as intended, manually checked output
	#no real need to reassess significance at this point, but this sub conveniently has a bunch of useful functions in it, and all signetwork pairs will be significant anyways
	#($sigpairs, $sigmcls, $sigSMcount, $sigpairclustercount, $sigcooccs) = Significance_testing($signetworkpairs, $nullmodel, $sigpvalueout);

	Print_SMcount($sigSMcount, $sigpairclustercount, $sigcooccs, $mcl2ann, $sig2ann, $smout);
	Print_network($sigpairclustercount, $mcl2ann, $sig2ann, $cytoout);
	Print_MCLs($sigmcls, $sigogout);
	
	#only needs to be run once to print out null model
	my $nullmodelout = "$outdir/null_model_distribution_binsizes.txt";
	if (! -f $nullmodelout) {
		open(my $nullout, '>', $nullmodelout) or usage("Error: cannae open $nullmodelout for printing\n");
		Print_null_model($nullmodel, $nullout);
	}
	
}

sub Print_null_model {
	my ($nullmodel, $nullout) = @_;
	foreach my $binsize1 (sort {$a<=>$b} keys %{$nullmodel}) {
		foreach my $binsize2 (sort {$a<=>$b} keys %{$nullmodel->{$binsize1}}) {
			my $obs_number = scalar @{$nullmodel->{$binsize1}->{$binsize2}};
			print $nullout "$binsize1\t$binsize2\t$obs_number\t".join(",",nsort @{$nullmodel->{$binsize1}->{$binsize2}})."\n";
		}
	}
}

sub Check_network_for_signature {
	my ($sigpairs, $signaturemcls) = @_;
	
	#split each pair into a 2d hash
	my %nodes;
	foreach my $pair (nsort keys %{$sigpairs}) {
		my ($mcl1, $mcl2) = split/-/, $pair;
		$nodes{$mcl1}{$mcl2} = 1;
		$nodes{$mcl2}{$mcl1} = 1;
	}
	#define self-contained networks; iterate through all nodes in each network to find signature mcls
	my (%observednodes, %networknodes);
	my $networkcount = 0;
	foreach my $startingnode (keys %nodes) { #start at random place in network
		next if (exists $observednodes{$startingnode}); #skip testing nodes that have already been tested
		$networkcount++;
		$networknodes{"network$networkcount"}{$startingnode} = 1;
		foreach my $neighbornode (keys %{$nodes{$startingnode}}) {
			$networknodes{"network$networkcount"}{$neighbornode} = 1;
		}
		my $neighborfound = 0;
		$neighborfound = 1 if (scalar keys %{$networknodes{"network$networkcount"}} >= 1); #this should always be true if initiating with mcl pairs, as each mcl group is necessarily found in a pair
		while ($neighborfound == 1) {
			$neighborfound = 0;
			foreach my $node (keys %{$networknodes{"network$networkcount"}}) {
				$observednodes{$node} = 1;
				foreach my $neighbornode (keys %{$nodes{$node}}) {
					if (not exists $networknodes{"network$networkcount"}{$neighbornode}) { #as soon as no new neighboring nodes are found, then the no new nodes can be found in network
						$networknodes{"network$networkcount"}{$neighbornode} = 1;
						$observednodes{$neighbornode} = 1;
						$neighborfound = 1;
					}
				}
			}
		}
	}
	#now search find networks that contain signature MCLs
	my %signaturenetworkpairs;
	foreach my $network (keys %networknodes) {
		my $signaturefound = 0;
		foreach my $mcl (keys %{$networknodes{$network}}) {
			if (exists $signaturemcls->{$mcl}) { #is this mcl assigned to a SMURF signature gene?
				$signaturefound = 1;
				last;
			}
		}
		if ($signaturefound == 1) {
			my %mcls = 	map {$_ => 1} (keys %{$networknodes{$network}});
			foreach my $pair (nsort keys %{$sigpairs}) {
				my ($mcl1, $mcl2) = split/-/, $pair;
				if (exists $mcls{$mcl1} || exists $mcls{$mcl2}) { #is one of the mcls in this network in the pair?
					$signaturenetworkpairs{$pair} = $sigpairs->{$pair};
				}
			}
		}
	}
	
	return(\%signaturenetworkpairs);	
}




sub Cooccurrence_parser {
	my ($cooccurrencefile) = @_;
	my (%coi);#cooccurrences of interest
	open(my $coocin, '<', $cooccurrencefile) or usage("Error: cannae open $cooccurrencefile\n");
	while (my $line = <$coocin>) {
		chomp $line;
		my ($mcl1, $mcl2, $cooccs, $mcl1count, $mcl2count, $mcl1unknown, $mcl2unknown, $mcl1freq, $mcl2freq, $SMstring, $clusterstring) = split/\t/, $line;
		my $mclpair = "$mcl1-$mcl2";
		$coi{$mclpair} = $line;
	}
	return(\%coi);
}



sub Print_MCLs {
	my ($sigmcls, $sigogout) = @_;
	foreach my $mcl (nsort keys %{$sigmcls}) {
		print $sigogout "$mcl\n";
	}
}



sub Retrieve_signature_annotation {
	my ($sigogs) = @_;
	my (%sig2ann);
	foreach my $mcl (keys %{$sigogs}) {
		my @sortedclasses =  (sort {$sigogs->{$mcl}->{$b} <=> $sigogs->{$mcl}->{$a} } keys %{$sigogs->{$mcl}});
		my $SMclass = shift @sortedclasses;
		$sig2ann{$mcl} = $SMclass;
	}
	return(\%sig2ann);
}

sub SMURF_parser {
	my $usage = qq/
sub SMURF_parser reads in a an SMURF cluster prediction file suuplied by JGI and 
returns a hash where {signatureMCLS}->{signatureclass} = count
Cluster IDs are formatted: omecode_C<number>\n\n/;
	die($usage, "Error: incorrect number of args\n") if (scalar @_ !=3);
	my ($clusterfile, $prot2og, $sigogs) = @_;
	my ($omecode) = fileparse($clusterfile, ".smclusters");
	open(my $clusterin, '<', $clusterfile) or die($usage, "Error: cannot open $clusterfile\n");
	while (my $line = <$clusterin>) {
		chomp $line;
		if ($line =~ m/^\s/) {
			$line =~ s/^\s+//;
			my ($geneid, $protid, $signatureclass) = split/\s/, $line;
			if ($signatureclass =~ m/\w/) { #only signature genes will have text
				$protid = "${omecode}_$protid";
				$sigogs->{$prot2og->{$protid}}->{$signatureclass}++;
			}
		}
	}
	return($sigogs);
}

sub Retrieve_MCL_annotation {
	my ($mclannfile) = @_;
	open(my $in, "<", $mclannfile) or usage("Error: could not open $mclannfile for reading\n");
	my (%mcl2ann);
	while(my $line = <$in>) {
		chomp $line;
		my ($mcl, $annstring) = split/\t/, $line;
		my @anns = split/,/, $annstring;
		my $topann = shift (@anns);
		$mcl2ann{$mcl} = $topann;
	}
	return(\%mcl2ann);
}


sub Null_model_distributions {
#Each mcl in a pair gets its own datapoint, depending on the OG sizes within the pair
#800 is the maximum bin size; using binsizes of 25, we have 496 unique combinations of bins (32 choose 2)
#output hash is structured {smaller OG}{larger OG} = # of cooccurrences
	my ($nullmodelfile) = @_;
	my (%nullmodel);
	my $binsize = 25;
	open(my $nullin, '<', $nullmodelfile) or usage("Error: cannot open $nullmodelfile\n");
	while (my $line = <$nullin>) {
		chomp $line;
		my ($mcl1, $mcl2, $cooccs, $mcl1count, $mcl2count, $mcl1unknown, $mcl2unknown, $mcl1freq, $mcl2freq) = split/\t/, $line;
		my ($bin1) = Determine_binsize($mcl1count, $mcl1unknown, $binsize);
		my ($bin2) = Determine_binsize($mcl2count, $mcl2unknown, $binsize);
		my @bins = ($bin1, $bin2);
		@bins = sort {$a<=>$b} @bins;
		push @{$nullmodel{$bins[0]}{$bins[1]}},$cooccs;
	}
	return(\%nullmodel);
}

sub Determine_binsize {
	#Determine the bin that an MCL group should fall into based on its size and the number of members that fall near a contig (unknowns)
	my ($mclcount, $mclunknown, $binsize) = @_;
	my $mclcorrected = $mclcount - $mclunknown;
	my ($quotient, $remainder) = quotient_remainder($mclcorrected, $binsize);
	my $bin = $quotient * $binsize; #determine which bin of OG sizes this one way mcl MCL falls into. Bins are determined by the highest value in the bin (e.g., 10, 20, 30 etc)
	if ($remainder != 0) {
		$bin = ($quotient + 1) * $binsize;
	}
	if ($bin > 800) {
		$bin = 800; #put all gene family sizes above 800 into their own bin
	}
	return($bin);
}

#from http://www.perlmonks.org/index.pl?node_id=981240
sub quotient_remainder {
    use integer;
    my( $dividend, $divisor ) = @_;
    my $quotient = $dividend / $divisor;
    my $remainder = $dividend % $divisor;
    return ( $quotient, $remainder );
}



sub Print_network {
	my ($sigpairs, $mcl2ann, $sig2ann, $cyto1out) = @_;
	print $cyto1out "source\ttarget\teweight\n";
	foreach my $mclpair (nsort keys %{$sigpairs}) {
		my ($mcl1, $mcl2) = split/-/, $mclpair;
		my ($mcl1annotation, $mcl2annotation) = Retrieve_mclpair_annotation($mcl1, $mcl2, $mcl2ann, $sig2ann);
		print $cyto1out "$mcl1 $mcl1annotation\t$mcl2 $mcl2annotation\t$sigpairs->{$mclpair}\n";
	}
}

sub Retrieve_mclpair_annotation {
	my ($mcl1, $mcl2, $mcl2ann, $sig2ann) = @_;
	my $mcl1annotation;
	if (exists $mcl2ann->{$mcl1}) {
		$mcl1annotation = $mcl2ann->{$mcl1};
	} else {
		$mcl1annotation = "unknown";
	}
	my $mcl2annotation;
	if (exists $mcl2ann->{$mcl2}) {
		$mcl2annotation = $mcl2ann->{$mcl2};
	} else {
		$mcl2annotation = "unknown";
	}
	$mcl1annotation = $sig2ann->{$mcl1} if (exists $sig2ann->{$mcl1});
	$mcl2annotation = $sig2ann->{$mcl2} if (exists $sig2ann->{$mcl2});
	return($mcl1annotation, $mcl2annotation);
}

sub Print_SMcount {
	my ($sigSMcount1, $sigSMmclpair, $sigcooccs, $mcl2ann, $sig2ann, $smout) = @_;
	print $smout "SMclass\tmcl1\tmcl2\tcooccs_in_SMclass\tcooccs_in_clusters\ttotal_cooccs\tmcl1_ann\tmcl2_ann\n";
	foreach my $SMclass (nsort keys %{$sigSMcount1}) {
		foreach my $mclpair (sort {$sigSMcount1->{$SMclass}->{$b} <=> $sigSMcount1->{$SMclass}->{$a}} keys %{$sigSMcount1->{$SMclass}}) {
			my ($mcl1, $mcl2) = split/-/, $mclpair;
			my ($mcl1annotation, $mcl2annotation) = Retrieve_mclpair_annotation($mcl1, $mcl2, $mcl2ann, $sig2ann);
			print $smout "$SMclass\t$mcl1\t$mcl2\t$sigSMcount1->{$SMclass}->{$mclpair}\t$sigSMmclpair->{$mclpair}\t$sigcooccs->{$mclpair}\t$mcl1annotation\t$mcl2annotation\t\n";
		}
	}
}

sub Significance_testing {
	my ($mclpairs, $nullmodel, $out, $allout) = @_;
	my (%sigpairs, %sigmcls, %sigSMcount, %sigpairclustercount, %sigcooccs);
	my $minimum_cooccurrences = 5;
	foreach my $mclpair (keys %{$mclpairs}) {
		my ($mcl1, $mcl2, $cooccs, $mcl1count, $mcl2count, $mcl1unknown, $mcl2unknown, $mcl1freq, $mcl2freq, $SMstring, $clusterstring) = split/\t/, $mclpairs->{$mclpair};
		my @clusters = split/,/, $clusterstring;
		my $clustercooccs = scalar @clusters; #how many times does this cooccurrence occur within a predicted cluster?

		my $mcl1corrected = $mcl1count - $mcl1unknown; #correct mcl size based on number of members found near contig ends
		my $mcl2corrected = $mcl2count - $mcl2unknown;
		next if ($cooccs <= $minimum_cooccurrences); #don't count any pairs with less cooccurrences as significant
		
		#print "Estimating pvalue using cooccurrences that occur globally\n";
		#my ($mclpvalue, $nulltotalobs) = Estimate_pvalue_from_null($mclpair, $cooccs, $mcl1corrected, $mcl2corrected, $nullmodel);
		
		#print "Estimating pvalue using cooccurrences that occur only in predicted clusters\n";
		my ($mclpvalue, $nulltotalobs) = Estimate_pvalue_from_null($mclpair, $clustercooccs, $mcl1corrected, $mcl2corrected, $nullmodel);
		
		my $sigfound = 0;
		if ($mclpvalue =~ m/\d/) { #estimations with warnings return 'n.d.'
			if ($mclpvalue <= 0.05) {
				$sigfound = 1;
			}
		}
		print $allout "$mcl1-$mcl2\t$cooccs\t$clustercooccs\t$mclpvalue\n";
		if ($sigfound == 1) {
			print $out "$mclpairs->{$mclpair}\t$clustercooccs\t$mclpvalue\t$nulltotalobs\n";
			$sigmcls{$mcl1} = 1;
			$sigmcls{$mcl2} = 1;
			$sigpairs{$mclpair} = $mclpairs->{$mclpair};
			$sigpairclustercount{$mclpair} = $clustercooccs;
			$sigcooccs{$mclpair} = $cooccs;
			my (@SMcounts) = split/,/, $SMstring;
			foreach my $SMcount (@SMcounts) {
				my ($SMclassstring, $SMcoocc) = split/=/, $SMcount;
				$sigSMcount{$SMclassstring}{$mclpair} = $SMcoocc;
			}
		}
	}
	return(\%sigpairs, \%sigmcls, \%sigSMcount, \%sigpairclustercount, \%sigcooccs);
}

sub Estimate_pvalue_from_null {
	#Sorting logic: sort mcl group size smallest to largest, and sort the binsizes in null model smallest to largest
	my ($mclpair, $cooccs, $mcl1count, $mcl2count, $nullmodel) = @_;
	my $pvalue;
	my $minimum_nullmodel_size = 10;
	my $maximum_mclgroup_size = 800;
	my @mclpairsize = ($mcl1count, $mcl2count);
	@mclpairsize = sort {$a<=>$b} @mclpairsize;
	#find the bin that the smallest MCL group falls into
	my $bin1;
	my $bin1found = 0;
	foreach my $querybin (sort {$a<=>$b} keys %{$nullmodel}) { #iterate through bin sizes from smallest to largest
		if ($bin1found == 0) { #once appropriate bin is found, don't search any other bins
			if ($mclpairsize[0] <= $querybin) { #works because we are iterating through bin sizes from smallest to largest
				$bin1found = 1;
				$bin1 = $querybin;
			}
			if ($mclpairsize[0] > $maximum_mclgroup_size) {
				$bin1found = 1;
				$bin1 = $maximum_mclgroup_size;
			}
		}
	}
	if ($bin1found == 0) {
		print "WARNING: could not find binsize1 for $mclpair mcl1=$mclpairsize[0], mcl2=$mclpairsize[1], coocc=$cooccs..skipping pvalue estimation\n";
		return("n.d.");
	}

	#Now using the smallest binsize, find the bin that the largest MCL group falls into
	my $bin2;
	my $bin2found = 0;
	foreach my $querybin (sort {$a<=>$b} keys %{$nullmodel->{$bin1}}) { #iterate through all bins found with bin1
		if ($bin2found == 0) { #once appropriate bin is found, don't search any other bins
			if ($mclpairsize[1] <= $querybin) { #works because we are iterating through bin sizes from smallest to largest
				$bin2found = 1;
				$bin2 = $querybin;
			}
			if ($mclpairsize[1] > $maximum_mclgroup_size) {
				$bin2found = 1;
				$bin2 = $maximum_mclgroup_size;
			}
		}
	}
	if ($bin2found == 0) { #this will only happen if bin2size was never sampled with bin1size
		print "WARNING: there are no observations to test $mclpair for null model bin1=$bin1, mcl2=$mclpairsize[1]..skipping pvalue estimation\n";
		return("n.d.");
	}
	#using the retrieved binsizes, find the number of cooccurrences in the null model greater than or equal to the observed number of cooccurrences
	if (not exists $nullmodel->{$bin1}->{$bin2}) {
		print "WARNING: there are no observations to test $mclpair for null model bin1=$bin1, bin2=$bin2..skipping pvalue estimation\n";
		return("n.d.");
	}
	my $totalobs = scalar @{$nullmodel->{$bin1}->{$bin2}};
	if ($totalobs <= $minimum_nullmodel_size) {
		print "WARNING: there are <= $minimum_nullmodel_size observations to test $mclpair for null model bin1=$bin1, bin2=$bin2..skipping pvalue estimation\n";
		return("n.d.");
	}
	my $gtcount = 0; #counts number of obs in null model that are greater than or equal to query value
	foreach my $obs (@{$nullmodel->{$bin1}->{$bin2}}) {
		if ($obs >= $cooccs) {
			$gtcount++;
		}
	}
	$pvalue = sprintf("%.5f", $gtcount / $totalobs);
	return($pvalue, $totalobs);
}


sub Find_all_MCL_pairs {
	my ($cluster2prot, $prot2og, $prothoods) = @_;
	my (%sampledmclpairs);
	foreach my $cluster (keys %{$cluster2prot}) {
		my %sampledprotpairs;
		foreach my $prot (keys %{$cluster2prot->{$cluster}}) {
			my $mcl1 = $prot2og->{$prot};
			my ($neighbors, $nearend) = Retrieve_all_neighbor_prots($prot, $prothoods); 
			foreach my $dist (nsort keys %{$neighbors}) {#find all matches
				foreach my $neighbor (keys %{$neighbors->{$dist}}) {
					my $mcl2 = $prot2og->{$neighbor};
					my @mclpair = ($mcl1, $mcl2);
					my $mclstring = join("-", nsort @mclpair);
					next if (exists $sampledmclpairs{$mclstring}); #since the purpose of this sub is to just identify mcl pairs within 6 of each other, skip pairs that have already been observed
					next if (not exists $cluster2prot->{$cluster}->{$neighbor}); #only look at prots that are part of cluster
					my @pair = ($prot, $neighbor);
					my $protpair = join("-", nsort @pair); #always sort prots before making them a string
					next if (exists $sampledprotpairs{$protpair}); #skip pairs of prots that have already been counted
					next if ($mcl1 eq $mcl2); #skip self cooccurrences
					$sampledmclpairs{$mclstring} = 1;
				}
			}
		}
	}
	print scalar(keys(%sampledmclpairs))." unique mcl-pairs have been retrieved from newSMURF clusters!\n";
	return(\%sampledmclpairs);
}

sub Find_and_print_all_cooccurrences {
	my ($sampledmclpairs, $og2prot, $prothoods, $prot2cluster, $prot2clusterclass, $coout) = @_;
	my %sampledcooccurrences; 
	foreach my $mclpair (keys %{$sampledmclpairs}) { #all mclpairs are unique
		#c. Find all instances of minimal clustering of mcl1 and mcl2, keeping track of which protids have already been counted as a co-occurring pair
		my ($cooccurrences, $unknown1, $unknown2) = (0, 0, 0);
		my ($mcl1, $mcl2) = split/-/, $mclpair;
		my %protids1 = map {$_ => 1} keys %{$og2prot->{$mcl1}};
		my %protids2 = map {$_ => 1} keys %{$og2prot->{$mcl2}};
		my $mcl1count = scalar keys %protids1;
		my $mcl2count = scalar keys %protids2;
		next if ($mcl1count == 1 || $mcl2count == 1); #don't count cooccurrences of singletons
		my (%clustercocount, %smcocount);
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
						if (exists $prot2clusterclass->{$prot}) { #is this prot and its pair in cluster that has had (at least partially) a SM class initially assigned from SMURF?
							$clustercocount{$prot2cluster->{$prot}}++; #only need to check one protein, because protid1 is necessarily in the same cluster, since it cooccurs within 6 intervening genes
							$smcocount{$prot2clusterclass->{$prot}}++;
						}
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
		my $datastring = "$mcl1\t$mcl2\t$cooccurrences\t$mcl1count\t$mcl2count\t$unknown1\t$unknown2\t$freq1\t$freq2\t";
		print $coout $datastring;
		foreach my $smclass (sort keys %smcocount) {
			print $coout "$smclass=$smcocount{$smclass},";
			$datastring .= "$smclass=$smcocount{$smclass},";
		}
		print $coout "\t";
		$datastring .= "\t";
		foreach my $smcluster (nsort keys %clustercocount) {
			print $coout "$smcluster=$clustercocount{$smcluster},";
			$datastring .= "$smcluster=$clustercocount{$smcluster},";
		}
		print $coout "\n";
		$sampledcooccurrences{$mclpair} = $datastring;
	}
	return(\%sampledcooccurrences);
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


sub Combined_parser {
	my ($combinedfile, $omesoi) = @_;
	my (%cluster2clusterclass); #internal use only
	my (%prot2clusterclass, %prot2cluster, %clusterogs, %cluster2prot);
	open(my $combin, '<', $combinedfile) or usage("Error: could not open $combinedfile for reading\n");
	while(my $line = <$combin>) {
		chomp $line;
		my ($ome, $clusterid, $protid, $SMclass, $contig, $coords, $strand, $mcl) = split/\t/, $line;
		if (not exists $omesoi->{$ome}){
			warn("$ome is not in the restricted ome file\n");
			next;
		}	
		$cluster2clusterclass{$clusterid}{$SMclass} = 1 if ($SMclass ne 'unknown'); #some clusters might have multiple SM classes assigned to them if smurf predicted clusters were joined in the newSMURF prediction
		$prot2cluster{$protid} = $clusterid;
		$clusterogs{$mcl} = 1;
		$cluster2prot{$clusterid}{$protid} = 1;
	}
	foreach my $cluster (keys %cluster2clusterclass) {
		my $SMstring = join("+", nsort keys %{$cluster2clusterclass{$cluster}});#some clusters might have multiple SM classes assigned to them if smurf predicted clusters were joined in the newSMURF prediction
		foreach my $prot (keys %{$cluster2prot{$cluster}}) {
			$prot2clusterclass{$prot} = $SMstring;
		}
	}

	return(\%prot2clusterclass, \%prot2cluster, \%cluster2prot, \%clusterogs);
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

sub Opts_check {
	my ($opts) = @_;
	usage() if (exists $opts->{'h'});
	usage("Error: please specify an output directory\n") if (not defined $opts->{'o'});
	usage("Error: cannot find the output directory.. is it in the right place?\n") if (! -d $opts->{'o'});
}
