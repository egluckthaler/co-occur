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
This script parses output from antiSMASH v4 and prints a combined 
file to be input to dotSM4.group_clusters_by_content.pl for diversity analyses, or to
dotSM5.cross_reference_combined.pl for methods validation.
NOTE: only clusters with > 1 gene will be printed out\n\n/;	
die($usage, $message);
}
 
main: {
	my %opts;
	getopt('o:h', \%opts);
	Opts_check(\%opts);
	my ($outdir) = ($opts{'o'});
	my $clusterdir = "/fs/project/PAS1046/projects/dotSM/antiSMASH";
	my $orthofile = "/fs/project/PAS1046/projects/dotSM/OrthoMCL/orthoMCL.output";
	my $ome2orderfile = "/fs/project/PAS1046/projects/dotSM/Metadata/ome2order";
	my $coodir = "/fs/project/PAS1046/projects/dotSM/coo3/";
	my $header = "#output by dotSM3.smash2combined.pl\n";
	my $combinedoutfile = "$outdir/smash2combined.txt";
	my ($combout) = Open_FH($combinedoutfile, $header);
	
	#grab list of omes of interest
	my ($ome2order) = dim_1_hash($ome2orderfile, "\t", "0:1");

#1. Read ortholog groups file into a hash where {protid => OG} and {OG=> protid => 1}
#	name clusters starting at MCL000001
	my ($prot2og, $og2prot) = OMCL_raw_parser($orthofile);

#2.	Retrieve all MCLs associated with signature SM genes, as identified by antiSMASH
#	structured {ome}->{clusterid}->{protid}->{mcl} = annotation
	print "Parsing cluster data from antiSMASH predictions in $clusterdir\n";
	my @dirs = glob("$clusterdir/*");
	my ($smashes, $sigprots);
	foreach my $omedir (@dirs) {
		my ($ome) = basename($omedir);
		if (not exists $ome2order->{$ome}) {
			print "$ome is not an ome of interest.. skipping..\n";
			next;
		}
		($smashes, $sigprots) = antiSMASH_parser($omedir, $prot2og, $smashes, $sigprots); #update smashes with each iteration
	}

#3.	Retrieve coordinates for all proteins in clusters; print out combined file
	#ome ome_cluster protein protein_type scaffold coords strand MCL
	foreach my $ome (sort keys %{$smashes}){
		if (not exists $ome2order->{$ome}) {
			print "$ome is not an ome of interest.. skipping..\n";
			next;
		}
		my %protids;
		foreach my $cluster (keys %{$smashes->{$ome}}) {
			foreach my $protid (keys %{$smashes->{$ome}->{$cluster}}) {
				foreach my $mcl (keys %{$smashes->{$ome}->{$cluster}->{$protid}}) {
					push@{$protids{$protid}}, $mcl, $smashes->{$ome}->{$cluster}->{$protid}->{$mcl};
				}
			}
		}
		my $coo_file = "$coodir/$ome.coo3";
		my ($hit_info) = Coo3_parse_by_prot($coo_file, \%protids);
		foreach my $cluster (nsort keys %{$smashes->{$ome}}) {
			foreach my $protid (nsort keys %{$smashes->{$ome}->{$cluster}}) {
				next if (scalar keys %{$smashes->{$ome}->{$cluster}} == 1); #skip clusters with only 1 protein
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

sub antiSMASH_parser {
	my $usage = qq/
sub antiSMASH_parser reads in antiSMASH cluster prediction file and 
returns a hash where {signatureMCLS}->{signatureclass} = count
Cluster IDs are formatted: omecode_C<number>\n\n/;
	die($usage, "Error: incorrect number of args\n") if (scalar @_ !=4);
	my ($clusterdir, $prot2og, $smashes, $sigprots) = @_;
	my ($omecode) = basename($clusterdir);
	my ($clusterid, $protid, $desc, $type, $clustertype, @sigproteins, %temp_prots);
	my ($stop, $clusterinfo) = (0,0);
	open(my $clusterin, '<', "$clusterdir/geneclusters.js") or die($usage, "Error: cannot open $clusterdir/geneclusters.js\n");
	while (my $line = <$clusterin>) {
		if ($stop == 0) {
			chomp $line;
			if ($line =~ m/cluster-(\d+)/) {
				$clusterinfo = 0;
				my $clustercount = $1;
				$clusterid = "${omecode}_smash$clustercount";
			} elsif ($line =~ m/"locus_tag": "([^\n]+)",/) {
				$protid = $1;
				if ($protid !~ m/jgi.p_/) {
					warn("WARNING: ${protid}'s name is not formatted as expected (line 134).. skipping..\n");
				} else {
					$protid =~ s/^jgi.p_//;
				}
			} elsif ($line =~ m/smCOG:\s+([^<]+)/) {
				$desc = $1;
			} elsif ($line =~ m/"type": "([^\n]+)",/ && $clusterinfo == 0) {
				$type = $1;
			} elsif ($line =~ m/"strand": /) { #end of protein info section, load up info
				if ($type !~ m/other/) { #control which proteins will be considered part of the cluster. In this case, do not load up proteins that were included just by existing within the window around signature gene
					$temp_prots{$protid} ="$type $desc";
					if ($type =~ m/^biosynthetic$/) { #keep track of signature biosynthetic genes
						push @sigproteins, $protid;
					}
				}
			} elsif ($line =~ m/"knowncluster": "/) {
				$clusterinfo = 1;
			} elsif ($line =~ m/"type": "([^\n]+)",/ && $clusterinfo == 1) { #this will be true once we have entered cluster info section
				$clustertype = $1;
				if ($clustertype !~ m/^cf_/) { #control types of clusters accepted. Skip clusterfinder clusters without any signature biosynthetic genes
					foreach my $clusterprot (keys %temp_prots) {
						$smashes->{$omecode}->{$clusterid}->{$clusterprot}->{$prot2og->{$clusterprot}} = $temp_prots{$clusterprot};
					}
					foreach my $sigprot (@sigproteins) {
						$sigprots->{$sigprot} = $clustertype;#once we know the cluster type, we can append it to sigprot hash
					}
				}
				%temp_prots = ();
				@sigproteins = ();
				$clusterinfo = 0;
			} elsif ($line =~ m/^var details_data/) { #do not parse file beyond this point
				$stop = 1;
			} 
		}
	}
	return($smashes, $sigprots);
}


sub Opts_check {
	my ($opts) = @_;
	usage() if (exists $opts->{'h'});
	usage("Error: please specify an output directory\n") if (not defined $opts->{'o'});
	usage("Error: cannot find the output directory.. is it in the right place?\n") if (! -d $opts->{'o'});
}
