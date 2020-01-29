#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use File::Basename;
use Data::Dumper;
use FileHandle;
use Sort::Naturally;
use Ebinf::Utils qw(Fasta_hash Coo3_parse_by_prot Dependency_check dim_1_hash dim_0_hash Ecofy_omes Open_FH);
$| = 1;

sub usage {
	my $message = shift;
	$message = '' if (not defined $message);
	my $usage = qq/
usage: $0 [OPTIONS]
-o <output_directory>
This script parses output from SMURF supplied by JGI and prints a combined 
file to be input to dotSM4.group_clusters_by_content.pl for diversity analyses, or to
dotSM5.cross_reference_combined.pl for methods validation.
NOTE: only clusters with > 1 gene will be printed out
 \n\n/;
	die($usage, $message);
}
 
main: {
	my %opts;
	getopt('o:h', \%opts);
	Opts_check(\%opts);
	my ($outdir) = ($opts{'o'});
	my $clusterdir = "/fs/project/PAS1046/projects/dotSM/SecondaryMetabolism";
	my $orthofile = "/fs/project/PAS1046/projects/dotSM/OrthoMCL/orthoMCL.output";
	my $ome2orderfile = "/fs/project/PAS1046/projects/dotSM/Metadata/ome2order";
	my $coodir = "/fs/project/PAS1046/projects/dotSM/coo3/";
	my $header = "#output by dotSM3.smurf2combined.pl\n";
	my $combinedoutfile = "$outdir/smurf2combined.txt";
	my ($combout) = Open_FH($combinedoutfile, $header);
	
	#grab list of omes of interest
	my ($ome2order) = dim_1_hash($ome2orderfile, "\t", "0:1");

#1. Read ortholog groups file into a hash where {protid => OG} and {OG=> protid => 1}
#	name clusters starting at MCL000001
	my ($prot2og, $og2prot) = OMCL_raw_parser($orthofile);

#2.	Retrieve all MCLs associated with signature SM genes, as identified by SMURF
#	structured {ome}->{clusterid}->{protid}->{mcl} = annotation
	print "Parsing cluster data from SMURF predictions in $clusterdir\n";
	my @files = glob("$clusterdir/*.smclusters");
	my ($smurfs, $sigprots);
	foreach my $clusterfile (@files) {
		($smurfs, $sigprots) = whole_SMURF_parser($clusterfile, $prot2og, $smurfs, $sigprots); #update smurfs with each iteration
	}

#3.	Retrieve coordinates for all proteins in clusters; print out combined file
	#ome ome_cluster protein protein_type scaffold coords strand MCL
	foreach my $ome (sort keys %{$smurfs}){
		if (not exists $ome2order->{$ome}) {
			print "$ome is not an ome of interest.. skipping..\n";
			next;
		}
		my %protids;
		foreach my $cluster (keys %{$smurfs->{$ome}}) {
			foreach my $protid (keys %{$smurfs->{$ome}->{$cluster}}) {
				foreach my $mcl (keys %{$smurfs->{$ome}->{$cluster}->{$protid}}) {
					push@{$protids{$protid}}, $mcl, $smurfs->{$ome}->{$cluster}->{$protid}->{$mcl};
				}
			}
		}
		my $coo_file = "$coodir/$ome.coo3";
		my ($hit_info) = Coo3_parse_by_prot($coo_file, \%protids);
		foreach my $cluster (nsort keys %{$smurfs->{$ome}}) {
			foreach my $protid (nsort keys %{$smurfs->{$ome}->{$cluster}}) {
				next if (scalar keys %{$smurfs->{$ome}->{$cluster}} == 1); #skip clusters with only 1 protein
				if (exists $sigprots->{$protid}) {
					print $combout "$ome\t$cluster\t$protid\t$sigprots->{$protid}\t".join("\t", @{$hit_info->{$protid}})."\n";
				} else {
					print $combout "$ome\t$cluster\t$protid\tunknown\t".join("\t", @{$hit_info->{$protid}})."\n";
				}
			}
		}
	}

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

sub whole_SMURF_parser {
	my $usage = qq/
sub SMURF_parser reads in a an SMURF cluster prediction file suuplied by JGI and 
returns a hash where {signatureMCLS}->{signatureclass} = count
Cluster IDs are formatted: omecode_C<number>\n\n/;
	die($usage, "Error: incorrect number of args\n") if (scalar @_ !=4);
	my ($clusterfile, $prot2og, $smurfs, $sigprot) = @_;
	my ($omecode) = fileparse($clusterfile, ".smclusters");
	my $clusterid;
	open(my $clusterin, '<', $clusterfile) or die($usage, "Error: cannot open $clusterfile\n");
	while (my $line = <$clusterin>) {
		chomp $line;
		if ($line =~ m/^Cluster_id: (\d+)/) {
			my $clustercount = $1;
			$clusterid = "${omecode}_smurf$clustercount";
		} elsif ($line =~ m/^\s/) {
			$line =~ s/^\s+//;
			$line =~ s/\s+/ /g;
			my ($geneid, $protid, $signatureclass, $domains) = split/\s/, $line;
			$protid = "${omecode}_$protid";
			if ($signatureclass =~ m/\w/) {
				$smurfs->{$omecode}->{$clusterid}->{$protid}->{$prot2og->{$protid}} = $signatureclass;
			} 
			if (defined $domains) { #only signature genes will have text here
				$sigprot->{$protid} = $signatureclass;
			}
		} 
	}
	return($smurfs, $sigprot);
}


sub Opts_check {
	my ($opts) = @_;
	usage() if (exists $opts->{'h'});
	usage("Error: please specify an output directory\n") if (not defined $opts->{'o'});
	usage("Error: cannot find the output directory.. is it in the right place?\n") if (! -d $opts->{'o'});
}
