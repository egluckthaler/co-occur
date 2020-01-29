#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use File::Basename;
use Data::Dumper;
use FileHandle;
use Sort::Naturally;
use Statistics::Descriptive;
use Ebinf::Utils qw(Dependency_check dim_1_hash dim_0_hash Open_FH);
use JSON::XS qw(encode_json decode_json);
use File::Slurp qw(read_file write_file);
$| = 1;

sub usage {
	my $message = shift;
	$message = '' if (not defined $message);
	my $usage = qq/
usage: $0 [OPTIONS]
-o: <output_directory>

This script cross-references hard-coded combined cluster files
(output by previous dotSM scripts) to assess a number of summary statistics, including
the frequency by size of clusters found in each file, a venn diagram for gene cluster
content according to signature loci present in each combined file\n\n/;
	die($usage, $message);
}

main: {
	my %opts;
	getopt('o', \%opts);
	Opts_check(\%opts);
	my ($outdir) = ($opts{'o'});
	my @combinedfiles = ("/fs/project/PAS1046/projects/dotSM/Results/nonsigSM_cluster_cooc_500001/dotSM4.clusters_with_significant_MCLs.combined.txt","/fs/project/PAS1046/projects/dotSM/Results/smurf/smurf2combined.txt","/fs/project/PAS1046/projects/dotSM/Results/antiSMASH/smash2combined.txt");
	my $filecount = 0;
	my $header = "#output by dotSM5.cross_reference_combined.pl using ".join(",", @combinedfiles)."\n";
	foreach my $combfile (@combinedfiles) {
		$filecount++;
		print "file$filecount is $combfile\n";
		$header .= "#file$filecount is $combfile\n";
		Dependency_check($combfile);
	}
	my $outfile = "$outdir/cross_reference_combined.txt";
	my $overout = Open_FH($outfile, $header);
	my $vennoutfile = "$outdir/venn_overlap_summary.txt";
	my $vennout = Open_FH($vennoutfile, $header);
	my $freqfile = "$outdir/cluster_size_frequency.txt";
	my $freqout = Open_FH($freqfile, $header);
	my $emapperfile = "/fs/project/PAS1046/projects/dotSM/Annotation/dotSM_combined.annotations";
	Dependency_check($emapperfile);
	
	#Extract all relevant info
	my ($clusters, $cluster2prots, $cluster2file, $prot2clusters, $prot2mcl) = Extract_cluster_info(\@combinedfiles);

	#Count frequency x cluster size
	Print_size_frequency($clusters, $cluster2file, $freqout);
	
	#Compare all clusters from all files to find sets of overlapping proteins
	#structured {set#}{protein} = 1;
	my $setsjsonfile = "$outdir/sets_of_overlapping_proteins.json";
	my $sets;
	if (! -f $setsjsonfile) {
		($sets) = Find_all_sets_of_proteins($prot2clusters, $cluster2prots);
		my $setjson = encode_json($sets); write_file($setsjsonfile, { binmode => ':raw' }, $setjson);
	} else {
		print "all sets of overlapping proteins have already been determined.. reading in $setsjsonfile..\n";
		my $setjson = read_file($setsjsonfile, { binmode => ':raw' }); $sets = decode_json($setjson);
	}

	my ($venn, $vennsets) = Print_cross_ref_file($sets, $prot2clusters, $cluster2file, $overout);
	#Cross-check all gene recovery: global venn diagram
	#Retrieve each prot's annotation and process(es)
	my ($prot2ann, $prot2process) = Retrieve_prot_annotation($emapperfile);
	Print_venn_files($venn, $vennsets, $prot2ann, $prot2mcl, $vennout);
}

sub Print_size_frequency {
	my ($clusters, $cluster2file, $freqout) = @_;
	my %freqs;
	foreach my $cluster (keys %{$clusters}) {
		my $size = scalar @{$clusters->{$cluster}};
		$freqs{$size}{$cluster2file->{$cluster}}++;
	}
	foreach my $size (nsort keys %freqs) {
		foreach my $file (nsort keys %{$freqs{$size}}) {
			print $freqout "$size\t$file\t$freqs{$size}{$file}\n";
		}
	}	
}

sub Print_venn_files {
	my ($venn, $vennsets, $prot2ann, $prot2mcl, $vennout) = @_;
	#first print summary stats
	print $vennout "#protein venn\n";
	foreach my $section (nsort keys %{$venn}) {
		print $vennout "#$section\t".scalar(keys %{$venn->{$section}})."\n";
	}
	print $vennout "#set (metacluster) venn\n";
	foreach my $section (nsort keys %{$vennsets}) {
		print $vennout "#$section\t".scalar(keys %{$vennsets->{$section}})."\n";
	}
	#now print out detailed report
	foreach my $section (nsort keys %{$venn}) {
		foreach my $protein (nsort keys %{$venn->{$section}}) {
			my $protannstring = "-";
			$protannstring = join(",", @{$prot2ann->{$protein}}) if (exists $prot2ann->{$protein});
			my $protmclstring = "-";
			$protmclstring = $prot2mcl->{$protein} if (exists $prot2mcl->{$protein});
			print $vennout "$section\t$protein\t$protannstring\t$protmclstring\n";
		}
	}
}

sub Print_cross_ref_file {
	my ($loci, $prot2cluster, $cluster2file, $out) = @_;
	print $out "set_ids\tfile_ids\tcluster_ids\tprotein_ids\n";
	my (%venndiagram, %vennsets);
	foreach my $set (nsort keys %{$loci}) {
		my %lines;
		foreach my $prot (keys %{$loci->{$set}}) {
			my $clusterstring = join(",", nsort keys %{$prot2cluster->{$prot}});
			my %files;
			foreach my $cluster (keys %{$prot2cluster->{$prot}}) {
				$files{"$cluster2file->{$cluster}"} = 1;
			}
			my $filestring = join(",", nsort keys %files);
			$venndiagram{$filestring}{$prot} = 1;
			$lines{$filestring}{$clusterstring}{$prot} = 1;
		}
		foreach my $filestring (nsort keys %lines) {
			$vennsets{$filestring}{$set} = 1;
			foreach my $clusterstring (keys %{$lines{$filestring}}) { #should only ever be one key here
				print $out "$set\t$filestring\t$clusterstring\t".join(",", nsort keys %{$lines{$filestring}{$clusterstring}})."\n";
			}
		}
	}
	return(\%venndiagram, \%vennsets);
}

sub Find_all_sets_of_proteins {
	#works as intended; checked output manually
	my ($prot2clusters, $cluster2prots) = @_;
	my (%sigsets, %examined_clusters);
	my $setcount = 0;
	foreach my $starting_cluster (keys %{$cluster2prots}) {
		#Skip recorded clusters HERE
		next if (exists $examined_clusters{$starting_cluster});
		
		#get your starting point 
		my %tempclusters;
		$tempclusters{$starting_cluster} = 1;
		my $last_cluster_count = 0;
		my $current_cluster_count = scalar keys %tempclusters;
		while ($last_cluster_count < $current_cluster_count) {
			#find all clusters with the proteins that are part of tempclusters
			foreach my $tempclust1 (keys %tempclusters) {
				foreach my $tempprot (keys %{$cluster2prots->{$tempclust1}}) { 
					foreach my $tempclust2 (keys %{$cluster2prots}) { #this step is very memory intensive
						next if ($tempclust1 eq $tempclust2); #skip self comparisons
						if (exists $cluster2prots->{$tempclust2}->{$tempprot}) {
							$tempclusters{$tempclust2} = 1;
						}
					}
				}
			}
			$last_cluster_count = scalar keys %tempclusters;
		
			#now find all prots in those clusters 
			my %tempprots;
			foreach my $tempclust1 (keys %tempclusters) {
				foreach my $prot (keys %{$cluster2prots->{$tempclust1}}) {
					$tempprots{$prot} = 1;
				}
			}
		
			#now find all clusters with those new prots, and add to tempclusters hash
			foreach my $tempprot (keys %tempprots) { 
				foreach my $tempclust1 (keys %{$prot2clusters->{$tempprot}}) {
					$tempclusters{$tempclust1} = 1;
				}
			}
			$current_cluster_count = scalar keys %tempclusters;
		}
		$setcount++;
		foreach my $cluster (keys %tempclusters) {
			$examined_clusters{$cluster} = 1;
			foreach my $prot (keys %{$cluster2prots->{$cluster}}) {
				$sigsets{"set$setcount"}{$prot} = 1;
			}
		}
	}
	return(\%sigsets);
}

sub Extract_cluster_info {
	my ($combinedfiles) = @_;
	my (%clusters, %cluster2prots, %cluster2file, %prot2clusters, %prot2mcl);
	my $filecount = 0;
	foreach my $file (@{$combinedfiles}) {
		$filecount++;
		my $tag = "file$filecount";
		open(my $filein, '<', $file) or usage("Error: cannae open $file\n");
		while(my $line = <$filein>) {
			next if ($line =~ m/^#/);
			chomp $line;
			my ($ome, $ome_cluster, $protein, $protein_type, $scaffold, $coords, $strand, $MCL) = split/\t/, $line;
			push @{$clusters{$ome_cluster}}, $protein;
			$cluster2prots{$ome_cluster}{$protein} = 1;
			$cluster2file{$ome_cluster} = "file$filecount";
			$prot2clusters{$protein}{$ome_cluster} = 1;
			$prot2mcl{$protein} = $MCL; #this will not change across combined files (or shouldnt for the analysis for which this script was created)
		}
	}
	return(\%clusters, \%cluster2prots, \%cluster2file, \%prot2clusters, \%prot2mcl);
}

sub Retrieve_prot_annotation {
	my ($emapperfile) = @_;
	my (%prot2ann, %prot2process);
	open(my $emapin, '<', $emapperfile) or usage("Error: could not open $emapperfile for reading\n");
	while(my $line = <$emapin>) {
		next if ($line=~ m/^#/);
		chomp $line;
		my @fields = split/\t/, $line;
		my ($prot, $process, $annotation) = ($fields[0], $fields[10], $fields[11]);
		push @{$prot2process{$prot}}, $process;
		push @{$prot2ann{$prot}}, $annotation;
	}
	return(\%prot2ann, \%prot2process);
}

sub Opts_check {
	my ($opts) = @_;
	usage() if (defined $opts->{'h'});
	usage("Error: please provide an output directory\n") if (not defined $opts->{'o'});
	usage("Error: cannot find the output directory\n") if (! -d $opts->{'o'});
}
