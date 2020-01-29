#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use File::Basename;
use Data::Dumper;
use FileHandle;
use Sort::Naturally;
use Ebinf::Utils qw(Blastdbp Coo3_parse dim_1_hash Dependency_check);
use Ebinf::Cluster_utils qw(Clusterfy);
$| = 1;

#ome is defined by everything up to the last _
sub usage {
	my $message = shift;
	$message = '' if (not defined $message);
	my $usage = qq/
usage: $0 [OPTIONS]
-m <raw_orthoMCL_output>
-s <directory_with_SMURF_cluster_predictions>
-c <directory_with_all_coo3_files>
-o <output_directory>
This script retrieves all proteins associated with MCL groups from 
SMURF cluster predictions (irregular format supplied by JGI; can be customized in future 
to accept other formats). All proteins are then re-clustered to obtain a 'boot-strapped'
set of clusters that is coherent across all SMURF cluster predictions. Also has the
benefit of including degraded clusters that might not have been picked up by SMURF but
are nonetheless important for counting co-occurrences
Outputs a combined table of all relevant info to be used in subsequent steps\n
Main output: dotSM_combined.txt\n\n/;
	die($usage, $message);
}
 
main: {

	my %opts;
	getopt('m:c:s:o:h', \%opts);
	Opts_check(\%opts);
	my ($orthofile, $clusterdir, $coodir, $outdir) = ($opts{'m'}, $opts{'s'}, $opts{'c'}, $opts{'o'});
	my $analysisname = "dotSM";
	#usage("${analysisname}_combined.txt already exists in $outdir\n") if (-f "$outdir/${analysisname}_combined.txt");
	open(my $combout, '>', "$outdir/${analysisname}_combined.txt") or usage("Error: cannot open $outdir/${analysisname}_combined.txt for printing\n");
	
	#grab list of omes of interest
	my $ome2orderfile = "/fs/project/PAS1046/projects/dotSM/Metadata/ome2order";
	my ($omesoi) = dim_1_hash($ome2orderfile, "\t", "0:1");

#1. Read ortholog groups file into a hash where {protid => OG}
#	name clusters starting at MCL000001
	my ($prot2og) = OMCL_raw_parser($orthofile);

#2. Read SMURF cluster prediction files supplied by JGI for each genome
#	structured $cluster2prot->{$clusterid}->{$protid} = clusterclass
#	spit out another hash: prot2clusterclass
#	This step can be swapped out for other cluster file parsers in the future
	my @files = glob("$clusterdir/*.smclusters");
	my ($cluster2prot, $prot2clusterclass);
	foreach my $clusterfile (@files) {
		($cluster2prot, $prot2clusterclass) = SMURF_parser($clusterfile, $cluster2prot, $prot2clusterclass); #update cluster2prot with each iteration
	}

#3. Cross-reference prot2og and cluster2prot to find all MCL groups in clusters. structured: $clusterogs{OG} = 1, $singletons{$ome}{$protid} = singleton
	my ($clusterogs, $singletons) = Retrieve_OG_from_SMURF($prot2og, $cluster2prot);

#4. Retrieve all proteins associated with cluster OGs. Structured {ome}{prot} = 1;
	my ($protsoi) = Retrieve_prot_from_OG($prot2og, $clusterogs);

#5. Add singletons to protsoi (any protein without a MCL group necessarily does not have a homolog present in another cluster.. ignore)
	#($protsoi) = Add_singletons_to_protsoi($singletons, $protsoi);

#6. Find clusters made of genes separated by no more than 6 intervening genes using all proteins in clustered OG
	foreach my $ome (sort keys %{$protsoi}) {
		if (not exists $omesoi->{$ome}) {
			print "$ome is not an ome of interest.. skipping..\n";
			next;
		}
		my $coo_file = "$coodir/$ome.coo3";
		my %protids = %{$protsoi->{$ome}}; #copy {protid => og/singleton} 
		my ($hit_info) = Coo3_parse($coo_file, \%protids);
		my ($clusters, $loners) = Clusterfy($hit_info); #structured clusters{$cluster_name}{$protid} = coo3_info; updated each iteration
		print "working on finding clusters of SMURF OGs in $ome.. ".scalar(keys(%{$clusters}))." clusters found!\n";
		Print_combined_table($clusters, $prot2clusterclass, $combout);
	}
	
	
	
}

sub Add_singletons_to_protsoi {
	my ($singletons, $protsoi) = @_;
	foreach my $ome (keys %{$protsoi}) {
		foreach my $prot (keys %{$singletons->{$ome}}) {
			push @{$protsoi->{$ome}->{$prot}}, 'singleton'; #for compliance with sub Coo3_parse
		}
	}
	return($protsoi);
}

sub Print_combined_table {
	my $usage = qq/
Print_combined_table(<\%clusters>, <\%prot2clusterclass>)
Print out the combined table where:
ome	ome_cluster	protein	protein_SMURF_class	scaffold coords strand MCL\n\n/;
	my ($clusters, $prot2clusterclass, $combout) = @_;
	foreach my $clusterid (nsort keys %{$clusters}) {
		foreach my $protid (nsort keys %{$clusters->{$clusterid}}) {
			$protid =~ m/^(.+?)_[^_]+$/; #grab everything up to last _
			my $ome = $1;
			my $protclass;
			if (defined $prot2clusterclass->{$protid}) {
				$protclass = $prot2clusterclass->{$protid};
			} else {
				$protclass = 'unknown';
			}
			print $combout "$ome\t$clusterid\t$protid\t$protclass\t".join("\t", @{$clusters->{$clusterid}->{$protid}})."\n";
		}
	}
}

sub Retrieve_prot_from_OG {
	my $usage = qq/
(<\%prots_of_interest) = Retrieve_OG_from_SMURF(<\%prot2og>, <\%clusterogs>)
sub Retrieve_prot_from_OG retrieves a list of non-redundant proteins that are part of an
OG present in SMURF-predicted clusters\n\n/;
	my ($prot2og, $clusterogs) = @_;
	my %protsoi;
	foreach my $prot (keys %{$prot2og}) {
		if (exists $clusterogs->{$prot2og->{$prot}}) { #is this prots OG found in a SMURF cluster?
			$prot =~ m/^(.+?)_[^_]+$/; #grab everything up to last _
			my $ome = $1;
			push @{$protsoi{$ome}{$prot}}, $prot2og->{$prot}; #for compliance with sub Coo3_parse
		}
	}
	return(\%protsoi);
}

sub Retrieve_OG_from_SMURF {
	my $usage = qq/
(<\%clusterogs>, <\%singletons>) = Retrieve_OG_from_SMURF(<\%prot2og>, <\%cluster2prot>)
sub Retrieve_OG_from_SMURF retrieves a list of non-redundant ortholog groups whose 
members are present in SMURF-predicted clusters. Also returns list of protids without OG\n\n/;
	my ($prot2og, $cluster2prot) = @_;
	my (%clusterogs, %singletons);
	foreach my $cluster (keys %{$cluster2prot}) {
		foreach my $prot (keys %{$cluster2prot->{$cluster}}) {
			if (exists $prot2og->{$prot}) { #is this prot part of an OG? if not, then singleton
				$clusterogs{$prot2og->{$prot}} = 1;
			} else {
				$prot =~ m/^(.+?)_[^_]+$/; #grab everything up to last _
				my $ome = $1;
				$singletons{$ome}{$prot} = 'singleton';
			}
		}
	}
	return(\%clusterogs, \%singletons);
}

sub SMURF_parser {
	my $usage = qq/
(<\%cluster2protid, \%prot2clusterclass) = SMURF_parser(<path_to_SMURF_output>, <\%existing_cluster2protid>, <\%existing_prot2clusterclass>)
sub SMURF_parser reads in a an SMURF cluster prediction file suuplied by JGI and 
returns a hash where {clusterID => protID => SM class} and a hash where {protID => SM class}
Cluster IDs are formatted: omecode_C<number>\n\n/;

	die($usage, "Error: incorrect number of args\n") if (scalar @_ != 3);
	my ($clusterfile, $cluster2prot, $prot2clusterclass) = @_;
	my ($omecode) = fileparse($clusterfile, ".smclusters");
	open(my $clusterin, '<', $clusterfile) or die($usage, "Error: cannot open $clusterfile\n");
	my ($clusterid, $SMclass);
	while (my $line = <$clusterin>) {
		chomp $line;
		if ($line =~ m/^Cluster_id: (\d+)/) {
			$clusterid = $1; #will be updated each time a new clusterid is passed
			$line =~ m/BackboneGene_class: (.+)$/;
			$SMclass = $1; #Reset SM class when new clusterID is observed
		} elsif ($line =~ m/^\s/) {
			$line =~ s/^\s+//;
			my ($geneid, $protid, $thirdcol) = split/\s/, $line;
			$protid = "${omecode}_$protid";
			$cluster2prot->{"${omecode}_c$clusterid"}->{$protid} = $SMclass;
			$prot2clusterclass->{$protid} = $SMclass;
		}
	}
	return($cluster2prot, $prot2clusterclass);
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
	my (%prot2og, $mclcount);
	while(my $line = <$omclin>) {
		$mclcount++;
		my ($formattedOG) = sprintf("%06d", $mclcount);
		$formattedOG = "MCL$formattedOG";
		chomp $line;
		my (@protids) = split/\s/, $line;
		foreach my $protid (@protids) {
			$protid =~ s/\|/_/;
			$prot2og{$protid} = $formattedOG;
		}
	}
	return(\%prot2og);
}


sub Opts_check {
	my ($opts) = @_;
	usage() if (exists $opts->{'h'});
	usage("Error: please provide an orthoMCL file\n") if (not defined $opts->{'m'});
	usage("Error: please provide a directory with cluster prediction files\n") if (not defined $opts->{'s'});
	usage("Error: please specify an output directory\n") if (not defined $opts->{'o'});
	usage("Error: please provide the coo3 directory\n") if (not defined $opts->{'c'});
	usage("Error: cannot find the orthoMCL file.. is it in the right place?\n") if (! -f $opts->{'m'});
	usage("Error: cannot find the directory with cluster prediction files.. is it in the right place?\n") if (! -d $opts->{'s'});
	usage("Error: cannot find the output directory.. is it in the right place?\n") if (! -d $opts->{'o'});
	usage("Error: cannot find the coo3 directory\n") if (! -d $opts->{'c'});	
}
