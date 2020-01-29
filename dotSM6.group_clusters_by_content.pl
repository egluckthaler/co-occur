#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use File::Basename;
use Data::Dumper;
use FileHandle;
use Sort::Naturally;
use Statistics::Descriptive;
use Ebinf::Utils qw(Fasta_hash Coo3_parse Dependency_check dim_1_hash dim_0_hash Ecofy_omes Open_FH);
use Ebinf::Cluster_utils qw(Clusterfy);
use JSON::XS qw(encode_json decode_json);
use File::Slurp qw(read_file write_file);
$| = 1;
sub usage {
	my $message = shift;
	$message = '' if (not defined $message);
	my $usage = qq/
usage: $0 [OPTIONS]
-s: <minimum_content_similarity_cutoff>
-c: <min_cluster_size_to_print_out_for_matrix>
-g: <min_group_size_to_print_out_for_matrix> number of clusters in group
-o: <output_directory>
-m: <mcl_file> (optional, if not -f) Restricts clusters to contain these mcl groups (for printing combined file)
-f: <combined_file> (optional, if not -m) output by some dotSM3 scripts

This script groups clusters containing mcls from the 
restricted MCL file based on a user-supplied minimum cutoff of gene content
similarity (I suggest 90% minimum content conservation). Number of minimum allowed matches
is rounded down, such that for 90% concent conservation, clusters sized 1-10, a minimum of 
1 gene is allowed to be missing; for clusters sized 11-20, a minimum of 2, etc. 
For every cluster size, algorithm does a preliminary assignment to groups and then uses
the 'best' scoring cluster (the one most similar to every other cluster) as the reference
cluster in the next iteration of cluster group assignment. Assumes bigger clusters are the 
better reference, so starts assessing similarities from largest clusters to smallest.
The restricted mcl file is output by dotSM4.test_all_co-occurrences.pl and should contain
all MCLs that unexpectedly cluster with another mcl group
This script is meant to follow dotSM4.test_all_co-occurrences.pl\n\n/;
	die($usage, $message);
}

#currently using dotSM3_all_cooccurrence_pvalues_significantMCLs.txt for mcl_file
#clustersize is determined by the number of unique mcls
#upon second iteration of sub Group_clusters_by_content, some reference clusters might be smaller than the size they are assigned to.. nevertheless, they are categorized as such to allow them to pick up clusters that are less than or equal to the 'theoretical' size

main: {
	my %opts;
	getopt('s:c:g:m:f:o', \%opts);
	Opts_check(\%opts);
	my ($MINSIM, $MINCLUSTERSIZE, $MINGROUPSIZE, $mclfile, $outdir, $newcombinedfile) = ($opts{'s'}, $opts{'c'}, $opts{'g'}, $opts{'m'}, $opts{'o'}, $opts{'f'});
	my $ome2orderfile = "/fs/project/PAS1046/projects/dotSM/Metadata/ome2order";
	my $ome2ecofile = "/fs/project/PAS1046/projects/dotSM/Metadata/ome2eco.simplified";
	my $orthofile = "/fs/project/PAS1046/projects/dotSM/OrthoMCL/orthoMCL.output";
	my $signaturemclfile = "/fs/project/PAS1046/projects/dotSM/Results/signatureMCLs.txt"; #not output by a standard script in pipeline, double check
	my $mclannfile = "/fs/project/PAS1046/projects/dotSM/Results/dotSM_OG_annotation.txt";
	my $raupcrickscript = "/users/PAS1046/osu9029/analyses/scripts/dotSM/dotSM_clustergroup_modeling.R";
	Dependency_check($ome2orderfile, $orthofile, $signaturemclfile, $ome2ecofile, $mclannfile, $raupcrickscript);
	my $coodir = "/fs/project/PAS1046/projects/dotSM/coo3/";
	
	my $fileheader = "#output by dotSM4.group_clusters_by_content.pl\n";
	my $clustergroupfile = "$outdir/dotSM4.cluster_groups_minsim${MINSIM}_mincl${MINCLUSTERSIZE}_mingr$MINGROUPSIZE.txt";
	my ($groupout) = Open_FH($clustergroupfile, $fileheader);
	my $raupmatrixfile = "$outdir/dotSM4.cluster_groups_minsim${MINSIM}_mincl${MINCLUSTERSIZE}_mingr${MINGROUPSIZE}_matrix_for_raupcrick.txt";
	my ($raupout) = Open_FH($raupmatrixfile, $fileheader);
	my $singletonfile = "$outdir/dotSM4.cluster_singletons_minsim${MINSIM}_mincl${MINCLUSTERSIZE}_mingr${MINGROUPSIZE}.txt";
	my ($singleout) = Open_FH($singletonfile, $fileheader);
	my $matrixfile = "$outdir/dotSM4.groups_minsim${MINSIM}_mincl${MINCLUSTERSIZE}_mingr${MINGROUPSIZE}.matrix";
	my ($matout) = Open_FH($matrixfile, $fileheader);
	my $meltedfile = "$outdir/dotSM4.groups_minsim${MINSIM}_mincl${MINCLUSTERSIZE}_mingr${MINGROUPSIZE}.melted.table";
	my ($meltout) = Open_FH($meltedfile);
	my $metafile = "$outdir/dotSM4.metagroups_minsim${MINSIM}_mincl${MINCLUSTERSIZE}_mingr${MINGROUPSIZE}.table";
	my ($metafileout) = Open_FH($metafile, $fileheader);
	my $metamatrixfile = "$outdir/dotSM4.metagroups_minsim${MINSIM}_mincl${MINCLUSTERSIZE}_mingr${MINGROUPSIZE}.matrix";
	my ($metamatout) = Open_FH($metamatrixfile, $fileheader);
	my $ecomatfile = "$outdir/dotSM4.ecology.matrix";
	my ($ecomatout) = Open_FH($ecomatfile, $fileheader);
	my $newnetworkfile = "$outdir/dotSM4.clusters_with_significant_MCLs.network.txt";
	my ($networkout) = Open_FH($newnetworkfile, $fileheader);

	my ($ome2order) = dim_1_hash($ome2orderfile, "\t", "0:1");

#0. Retrieve the best SM class associated with signature MCLs
	my ($signaturemcl2class) = Retrieve_best_SM_class($signaturemclfile);

#0. Retrieve the most frequent MCL annotation
	my ($mcl2ann) = Retrieve_MCL_annotation($mclannfile);

#steps 1-2 taken from dotSM_cluster_all_SM_prots.pl, with minor changes to Print_combined_table
#1. Read ortholog groups file into a hash where {protid => OG}
#	name clusters starting at MCL000001
#2. Find clusters made of genes separated by no more than 6 intervening genes using all proteins in clustered OG
	if (! defined $newcombinedfile && defined $mclfile) { #if no combined file supplied, print one out
		#0. Retrieve all restricted mcls
		my ($mcls) = dim_0_hash($mclfile, "\t", "0");
		my ($prot2og) = OMCL_raw_parser($orthofile);
		my ($protsoi) = Retrieve_prot_from_OG($prot2og, $mcls);
		$newcombinedfile = "$outdir/dotSM4.clusters_with_significant_MCLs.combined.txt";
		my ($combout) = Open_FH($newcombinedfile, $fileheader);
		foreach my $ome (sort keys %{$protsoi}) {
			if (not exists $ome2order->{$ome}) {
				print "$ome is not an ome of interest.. skipping..\n";
				next;
			}
			my $coo_file = "$coodir/$ome.coo3";
			my %protids = %{$protsoi->{$ome}}; #copy {protid => og/singleton} 
			my ($hit_info) = Coo3_parse($coo_file, \%protids);
			my ($clusters, $loners) = Clusterfy($hit_info); #structured clusters{$cluster_name}{$protid} = coo3_info; updated each iteration
			print "working on finding clusters of signficant MCLs in $ome.. ";
			Print_combined_table($clusters, $signaturemcl2class, $combout);
		}
	} else {
		print "combined file $newcombinedfile supplied.. skipping printing combined file\n";
	}
	
	
#3. Retrieve all clusters from combined.txt file (steps 3-5 from original group_clusters_by_content.pl)
#	structure: {size}{clusterid}{mcl} = mclcount
#	in this version of script, size is dictated by the number of UNIQUE mcl groups present in a cluster
	print "Retrieving all clusters from $newcombinedfile\n";
	my ($cluster2clusterclass, $clusters) = Combined_parser($newcombinedfile, $ome2order);


#4. Do a preliminary grouping of the clusters by similarity and retrieve the best 'reference' cluster for each group
#	The second passing of $clusters is to allow any cluster within a given size to be the reference cluster
	print "Retrieving references for cluster groups based on minsim of $MINSIM for ortholog group content..\n";
	my ($prelimgroups, $grouprefs) = Group_clusters_by_content($clusters, $clusters, $MINSIM);

#4. Group clusters based on user supplied cutoff for ortholog group content similarity
#	structure: {group}{cluster}{mcl}++
	print "Grouping clusters by comparing against references based on minsim of $MINSIM for ortholog group content..\n";
	my ($clustergroups, $finalgrouprefs) = Group_clusters_by_content($clusters, $grouprefs, $MINSIM);


#5. Print out singleton groups and groups with >1 member, as well as the putative annotation of
#	structure groupsummar= {smstring_groupname}{mcl}++
	my ($omeclustergroups, $allgroups, $groupsummary, $group2size) = Print_cluster_groups($clustergroups, $cluster2clusterclass, $mcl2ann, $groupout, $singleout);

#6. Print out network of newly predicted clusters
#	Find all ortholmcl groups in all clusters. output structured as {cluster}{mcl} = 1 (does not count multiple instances of same MCL in cluster)	
	my ($cluster2mcl) = Retrieve_MCL_from_cluster($newcombinedfile);
	
#	Count all MCL co-occurences. output structured as {mcl1}{mcl2}++
	my ($cooccur) = Count_MCL_cooccurrences($cluster2mcl);
	
#	Print out network
	Print_cooccurrences($cooccur, $mcl2ann, $networkout);

#6. Print out hierarchical clustering tree of clusters
#	Add option to retain groups that have been selected a priori (e.g., those groups that contain loci that hit to MIBIG clusters)
#	%mibiggroups was determined by manual inspection
	#my %mibiggroups = ("group41_s6" => 1,"group1_s16" => 1,"group4_s12" => 1,"group6_s10" => 1,"group7_s10" => 1,"group10_s9" => 1,"group35_s7" => 1,"group61_s5" => 1,"group101_s4" => 1,"group180_s3" => 1,"group42_s6" => 1,"group82_s5" => 1,"group99_s4" => 1,"group164_s3" => 1,"group121_s4" => 1,"group180_s3" => 1,"group180_s3" => 1,"group18_s8" => 1,"group89_s4" => 1,"group244_s2" => 1,"group101_s4" => 1,"group161_s3" => 1,"group217_s3" => 1);
	my %mibiggroups;
	Print_clustergroup_matrix($groupsummary, $group2size, $MINCLUSTERSIZE, $MINGROUPSIZE, \%mibiggroups, $raupout);
	my ($failcheck) = system("Rscript --vanilla $raupcrickscript $raupmatrixfile");
	usage("Error: could not execute $raupcrickscript\n") if ($failcheck != 0);
	my $labelorderfile = "$raupmatrixfile.labelorder";
	my ($clusterlabelorder) = Parse_label_order($labelorderfile); #parse order of labels on cluster group tree in order to print out presence/absence matrix in the same order

#8. Print out a ome to eco matrix to juxtapose next to species tree
	my ($ome2eco) = Ecofy_omes($ome2ecofile);
	Print_eco_matrix($ome2eco, $ecomatout);

	
#7. Create count matrix for juxtaposing against species tree; sort order of clusters in clustergrouptree
	print "Now grouping omes by cluster group repertoire using raup-crick\n";
	Print_group_matrix($omeclustergroups, $clusterlabelorder, $ome2order, $matout);
#8. Print out raup-crick tree for ome by cluster group
	($failcheck) = system("Rscript --vanilla $raupcrickscript $matrixfile");
	usage("Error: could not execute $raupcrickscript\n") if ($failcheck != 0);
#9.	Print out file for visualizing multi-cluster profiles
	my $omelabelorderfile = "$matrixfile.labelorder";
	my ($omelabelorder) = Parse_label_order($omelabelorderfile); #parse order of labels on metaprofile tree in order to print out presence/absence matrix in the same order
	my ($metaprofiletable) = "$matrixfile.table";
	my ($ome2metaprofile) = Parse_raup_table($metaprofiletable);
	Print_melted_multicluster_file($ome2order, $ome2eco, $ome2metaprofile, $omeclustergroups, $omelabelorder, $clusterlabelorder, $meltout);

#9. Print out file to juxtapose taxonomy vs. metaprofile tree


#9. Print metagroups matrix (rows = omes; columns = metagroups, as determined by raup-crick)
	my $metagroupsfile = "$raupmatrixfile.table"; #output by line 138, $raupcrickscript
	my ($group2meta) = Parse_raup_table($metagroupsfile); #renames metagroups from largest to smallest
	Print_meta_file($group2meta, $metafileout);
	Print_meta_matrix($omeclustergroups, $group2meta, $ome2order, $metamatout);

#10. Print out raup-crick tree for ome by metagroup
	print "Now grouping omes by cluster metagroup repertoire using raup-crick\n";
	($failcheck) = system("Rscript --vanilla $raupcrickscript $metamatrixfile");
	usage("Error: could not execute $raupcrickscript\n") if ($failcheck != 0);

}

sub Print_melted_multicluster_file {
	my ($ome2order, $ome2eco, $ome2metaprofile, $omeclustergroups, $omelabelorder, $clusterlabelorder, $meltedout) = @_;
	#grab all ecologies present among genomes
	my %ecologies;
	foreach my $ome (keys %{$ome2eco}) {
		foreach my $ecology (keys %{$ome2eco->{$ome}}) {
			$ecologies{$ecology} = 1;
		}
	}
	#grab all taxonomic orders present among genomes
	my %orders;
	foreach my $ome (keys %{$ome2order}) {
		$orders{$ome2order->{$ome}} = 1;
	}
	
	print $meltedout "omecode\tprofile\tcluster_group\tcount\n";
	foreach my $ome (reverse @{$omelabelorder}) {
		$ome =~ m/^(.+)_[^_]+$/; #some of the earlier figure printed out taxonomic order alongside omecode
		my $omecode = $1;
		foreach my $clustergroup (@{$clusterlabelorder}) {
			print $meltedout "$omecode\t$ome2metaprofile->{$ome}\t$clustergroup\t";
			if (exists $omeclustergroups->{$omecode}->{$clustergroup}) {
				print $meltedout "$omeclustergroups->{$omecode}->{$clustergroup}\n";
			} else {
				print $meltedout "0\n";
			}
		}
		foreach my $order (nsort keys %orders) {
			print $meltedout "$omecode\t$ome2metaprofile->{$ome}\t$order\t";
			if ($ome2order->{$omecode} eq $order) {
				print $meltedout "1\n";
			} else {
				print $meltedout "0\n";
			}
		}
		foreach my $ecology (nsort keys %ecologies) {
			print $meltedout "$omecode\t$ome2metaprofile->{$ome}\t$ecology\t";
			if (exists $ome2eco->{$omecode}->{$ecology}) {
				print $meltedout "1\n";
			} else {
				print $meltedout "0\n";
			}
		}
			
			
	}
	#print out ome2order in melted style
}

sub Print_meta_matrix {
	my ($omeclustergroups, $group2meta, $ome2order, $matout) = @_;
	#derive the ome2metagroups hash
	my (%omemetagroups, %metagroups);
	foreach my $ome (keys %{$omeclustergroups}) {
		foreach my $group (keys %{$omeclustergroups->{$ome}}) {
			if (exists $group2meta->{$group}) {
				$omemetagroups{$ome}{$group2meta->{$group}}++;
				$metagroups{$group2meta->{$group}} = 1;
			}
		}
	}
	#print all headers
	print $matout "Omecodes";
	foreach my $metagroup (nsort keys %metagroups) {
		print $matout "\t$metagroup";
	}
	print $matout "\n";
	#print out counts
	foreach my $ome (nsort keys %omemetagroups) {
		print $matout "${ome}_$ome2order->{$ome}";
		foreach my $metagroup (nsort keys %metagroups) { #iterate through all metagroups present in group2meta hash
			if (exists $omemetagroups{$ome}{$metagroup}) {
				print $matout "\t$omemetagroups{$ome}{$metagroup}";
			} else {
				print $matout "\t0";
			}
		}
		print $matout "\n";
	}
}

sub Print_meta_file {
	my ($group2meta, $metafileout) = @_;
	foreach my $group (nsort keys %{$group2meta}) {
		print $metafileout "$group2meta->{$group}\t$group\n";
	}
}

sub Parse_raup_table {
	my ($metagroupsfile) = @_;
	my (%group2meta, %meta2count);
	open(my $meta_in, '<', $metagroupsfile) or usage("Error: could not open $metagroupsfile for reading\n");
	my $header = <$meta_in>;
	while (my $line = <$meta_in>) {
		chomp $line;
		$line =~ s/\"//g;
		my ($group, $meta) = split/\t/, $line;
		$group2meta{$group} = $meta;
		$meta2count{$meta}++;
	}
	#assign new metanames based on number of groups in the metagroup
	my %meta2newmeta;
	my $newmeta = 1;
	foreach my $meta (sort {$meta2count{$b} <=> $meta2count{$a}} keys %meta2count) {
		$meta2newmeta{$meta} = "metagroup_$newmeta";
		$newmeta++;
	}
	#rename each meta according to the order of its prevalence
	foreach my $group (keys %group2meta) {
		my $newmeta = $meta2newmeta{$group2meta{$group}};
		$group2meta{$group} = $newmeta;
	}
	return(\%group2meta);
}


sub Print_cooccurrences {
	my ($cooccur, $mcl2ann, $cytoout) = @_;
	print $cytoout "source\ttarget\teweight\n";
	my (%sorted);
	foreach my $mcl1 (nsort keys %{$cooccur}) {
		my $mcl1annotation;
		if (exists $mcl2ann->{$mcl1}) {
			$mcl1annotation = $mcl2ann->{$mcl1};
		} else {
			$mcl1annotation = "unknown";
		}
		foreach my $mcl2 (nsort keys %{$cooccur->{$mcl1}}) {
			my $mcl2annotation;
			if (exists $mcl2ann->{$mcl2}) {
				$mcl2annotation = $mcl2ann->{$mcl2};
			} else {
				$mcl2annotation = "unknown";
			}
			$sorted{"$mcl1 $mcl1annotation\t$mcl2 $mcl2annotation"} = $cooccur->{$mcl1}->{$mcl2};
		}
	}
	foreach my $pair (sort {$sorted{$b} <=> $sorted{$a}} keys %sorted) { #sort by cooccurrence frequency
		print $cytoout "$pair\t$sorted{$pair}\n";
	}
}

sub Count_MCL_cooccurrences {
	my ($cluster2mcl) = @_;
	my (%cooccur);
	foreach my $cluster (keys %{$cluster2mcl}) {
		my @mcls = nsort keys %{$cluster2mcl->{$cluster}};
		while (scalar @mcls > 1) { #only count each cooccurrence once: by sorting mcls smallest to largest, this is ensured
			my $lowermcl = shift @mcls;
			foreach my $highermcl (@mcls) {
				$cooccur{$lowermcl}{$highermcl}++;
			}
		}
	}
	return(\%cooccur);
}

sub Retrieve_MCL_from_cluster {
	my ($combinedfile) = @_;
	my (%cluster2mcl);
	open(my $combin, '<', $combinedfile) or usage("Error: could not open $combinedfile for reading\n");
	while(my $line = <$combin>) {
		next if ($line =~ m/^#/);
		chomp $line;
		my ($ome, $clusterid, $protid, $SMclass, $contig, $coords, $strand, $mcl) = split/\t/, $line;
		$cluster2mcl{$clusterid}{$mcl} = 1;
	}
	return(\%cluster2mcl);
}

sub Parse_label_order{
	#top to bottom in this file corresponds to left to right in presence absence matrix
	my ($labelorderfile) = @_;
	my @order;
	open(my $labin, '<', $labelorderfile) or usage("Error: cannae open $labelorderfile for reading\n");
	my $header = <$labin>;
	while(my $line = <$labin>) {
		chomp $line;
		$line =~ s/"//g;
		my ($order, $groupname) = split/\t/, $line;
		push @order, $groupname;
	}
	return(\@order);
}

sub Retrieve_group_name {
	my ($group) = @_;
	my @fields = split/_/, $group;
	my $size = pop @fields;
	my $groupnum = pop @fields;
	return("${groupnum}_$size");
}

sub Print_clustergroup_matrix {
	my ($groupsummary, $group2size, $MINCLUSTERSIZE, $MINGROUPSIZE, $mibiggroups, $clustergroupout) = @_;
	#retrieve set of mcls present in groups
	my %mcls;
	foreach my $group (keys %{$groupsummary}) {
		my ($clustersize) = Retrieve_cluster_size($group);
		my ($groupname) = Retrieve_group_name($group);
		if ($clustersize >= $MINCLUSTERSIZE || exists $mibiggroups->{$groupname}) {
			if ($group2size->{$group} >= $MINGROUPSIZE || exists $mibiggroups->{$groupname}) {
				foreach my $mcl (keys %{$groupsummary->{$group}}) {
					$mcls{$mcl} = 1;
				}
			}	
		}
	}
	#print header
	print $clustergroupout "Clustergroups";
	foreach my $mcl (nsort keys %mcls) {
		print $clustergroupout "\t$mcl";
	}
	print $clustergroupout "\n";
	#print out each cluster groups profile (don't print out actual mcl counts, just 1 for presence)
	foreach my $group (nsort keys %{$groupsummary}) {
		my ($clustersize) = Retrieve_cluster_size($group);
		my ($groupname) = Retrieve_group_name($group);
		if ($clustersize >= $MINCLUSTERSIZE || exists $mibiggroups->{$groupname}) {
			if ($group2size->{$group} >= $MINGROUPSIZE || exists $mibiggroups->{$groupname}) {
				print $clustergroupout "$group";
				foreach my $mcl (nsort keys %mcls) {
					if (exists $groupsummary->{$group}->{$mcl}) {
						print $clustergroupout "\t1";
					} else {
						print $clustergroupout "\t0";
					}
				}
				print $clustergroupout "\n";
			}
		}
	}
}

sub Print_eco_matrix {
	my ($ome2eco, $ecomatout) = @_;
	print $ecomatout "#Names";
	#grab headers
	my %allecos;
	foreach my $ome (keys %{$ome2eco}) {
		foreach my $eco (keys %{$ome2eco->{$ome}}) {
			$allecos{$eco} = 1;
		}
	}
	foreach my $eco (nsort keys %allecos) {
		print $ecomatout "\t$eco";
	}
	print $ecomatout "\n";
	#print out values;
	foreach my $ome (nsort keys %{$ome2eco}) {
		print $ecomatout "$ome";
		foreach my $eco (nsort keys %allecos) {
			if (exists $ome2eco->{$ome}->{$eco}) {
				print $ecomatout "\t1";
			} else {
				print $ecomatout "\t0";
			}
		}
		print $ecomatout "\n";
	}
}

sub Print_combined_table {
	my $usage = qq/
Print_combined_table(<\%clusters>, <\%prot2clusterclass>, <out_file_handle>)
Print out the combined table where:
ome	ome_cluster	protein	mcl_SMURF_class	scaffold coords strand MCL\n\n/;
	my ($clusters, $signaturemcl2class, $combout) = @_;
	my $clustercount = 0;
	foreach my $clusterid (nsort keys %{$clusters}) {
		my %mcls;
		foreach my $protid (nsort keys %{$clusters->{$clusterid}}) {
			#check that > 1 type of mcl is present in each cluster
			my $protmcl = ${$clusters->{$clusterid}->{$protid}}[-1]; #only 1 mcl will ever be assigned per protein
			$mcls{$protmcl} = 1;
		}
		if (scalar keys %mcls > 1) {
			$clustercount++;
			foreach my $protid (nsort keys %{$clusters->{$clusterid}}) {
				$protid =~ m/^(.+?)_[^_]+$/; #grab everything up to last _
				my $ome = $1;
				my $protmcl = pop @{$clusters->{$clusterid}->{$protid}}; #only 1 mcl will ever be assigned per protein
				my $protclass;
				if (defined $signaturemcl2class->{$protmcl}) {
					$protclass = $signaturemcl2class->{$protmcl};
				} else {
					$protclass = 'unknown';
				}
				print $combout "$ome\t$clusterid\t$protid\t$protclass\t".join("\t", @{$clusters->{$clusterid}->{$protid}})."\t$protmcl\n";
			}
		}
	}
	print "found $clustercount clusters containing >1 type of mcl!\n";
}

sub Retrieve_best_SM_class {
	my ($signaturefile) = @_;
	open(my $sigin, '<', $signaturefile) or usage("Error: cannae open $signaturefile..\n");
	my %signaturemcl2class;
	my %temp;
	while(my $line = <$sigin>) {
		chomp $line;
		my ($mcl, $ann, $count) = split/\t/, $line;
		$temp{$mcl}{$ann} = $count;
	}
	foreach my $mcl (keys %temp) {
		foreach my $ann (sort {$temp{$mcl}{$b} <=> $temp{$mcl}{$a}} keys %{$temp{$mcl}}) { #sort by annotation count
			$signaturemcl2class{$mcl} = $ann;
			last;
		}
	}
	return(\%signaturemcl2class);
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


sub Print_group_matrix {
	my ($omeclustergroups, $labelorder, $ome2order, $matout) = @_;
	#print all headers
	print $matout "Omecodes";
	foreach my $group (@{$labelorder}) {
		print $matout "\t$group";
	}
	print $matout "\n";
	#retrieve set of omes with at LEAST 1 cluster (don't want any rows of zeroes)
	my %goodomes;
	foreach my $ome (nsort keys %{$omeclustergroups}) {
		foreach my $group (@{$labelorder}) {
			if (exists $omeclustergroups->{$ome}->{$group}) {
				$goodomes{$ome}++;
			}
		}
	}
	#print out counts of omes with at least 1 cluster
	foreach my $ome (nsort keys %{$omeclustergroups}) {
		if (exists $goodomes{$ome}) {
			print $matout "${ome}_$ome2order->{$ome}";
			foreach my $group (@{$labelorder}) {
				if (exists $omeclustergroups->{$ome}->{$group}) {
					print $matout "\t$omeclustergroups->{$ome}->{$group}";
				} else {
					print $matout "\t0";
				}
			}
			print $matout "\n";
		} else {
			print "$ome does not have any clusters that made the cutoff\n";
		}
	}
}
 
sub Retrieve_cluster_size {
	my ($clustergroup) = @_;
	my @terms = split/_/, $clustergroup;
	my $clustersize = pop @terms;
	$clustersize =~ s/s//;
	return($clustersize);
}

sub Print_cluster_groups {
	my ($clustergroups, $cluster2clusterclass, $mcl2ann, $groupout, $singleout) = @_;
	my ($groupcount, $singlecount) = 0;
	my (%omeclustergroups, %allgroups, %groupsummary, %group2size);
	foreach my $group (nsort keys %{$clustergroups}) {
		my (%mclanns, %mclsummary, %smclasssummary); #does not reflect mcl counts
		my ($originalgroupnumber, $clustersize) = split/_/, $group;
		if (scalar keys %{$clustergroups->{$group}} > 1) {
			$groupcount++;
			my $groupname = "group${groupcount}_$clustersize";
			foreach my $cluster (nsort keys %{$clustergroups->{$group}}) {
				print $groupout "$groupname\t$cluster\t";
				my $smclass = "n.d.";
				$smclass = join("+", sort keys %{$cluster2clusterclass->{$cluster}}) if (exists $cluster2clusterclass->{$cluster});
				print $groupout "$smclass\t";
				$smclasssummary{$smclass} = 1; 
				foreach my $mcl (nsort keys %{$clustergroups->{$group}->{$cluster}}) {
					my $mclann = "$mcl unknown";
					$mclann = "$mcl $mcl2ann->{$mcl}" if (exists $mcl2ann->{$mcl});
					$mclanns{$mclann}++;
					$mclsummary{$mcl}++;
					my $mclcount = $clustergroups->{$group}->{$cluster}->{$mcl};
					print $groupout "$mcl," x $mclcount;
				}
				print $groupout "\n";
			}
			my %smreduced; #reduce redundancy in composite sm terms
			foreach my $class (keys %smclasssummary) {
				my (@composites) = split/\+/, $class;
				foreach my $term (@composites) {
					$smreduced{$term} = 1;
				}
			}
			my $smstring = join("+", nsort keys %smreduced);
			foreach my $mcl (keys %mclsummary) {
				$groupsummary{"${smstring}_$groupname"}{$mcl} = $mclsummary{$mcl}; #copy mcl counts over to group summary once the sm string has been found
			}
			$allgroups{$smstring}{$groupname}++;
			#assign each cluster to the new groupname
			foreach my $cluster (nsort keys %{$clustergroups->{$group}}) {
				$cluster =~ m/^(.+?)_[^_]+$/; #ome is defined by up to last _
				my $ome = $1;
				$omeclustergroups{$ome}{"${smstring}_$groupname"}++;
			}
			print $groupout "#$groupname\tsummary\t$smstring\t".join(",", sort keys %mclanns)."\n";
			$group2size{"${smstring}_$groupname"} = scalar(keys(%{$clustergroups->{$group}}));
		} else {
			$singlecount++;
			foreach my $cluster (nsort keys %{$clustergroups->{$group}}) {
				$cluster =~ m/^(.+?)_[^_]+$/; #ome is defined by up to last _
				my $ome = $1;
				print $singleout "singleton${singlecount}_$clustersize\t$cluster\t";
				my $smclass = "n.d.";
				$smclass = join("+", sort keys %{$cluster2clusterclass->{$cluster}}) if (exists $cluster2clusterclass->{$cluster});
				print $singleout "$smclass\t";
				foreach my $mcl (nsort keys %{$clustergroups->{$group}->{$cluster}}) {
					my $mclcount = $clustergroups->{$group}->{$cluster}->{$mcl};
					print $singleout "$mcl," x $mclcount;
				}
				print $singleout "\n";
				#unhash next 9 lines to include singletons in the groupsummary hash
				$group2size{"${smclass}_singleton${singlecount}_$clustersize"} = 1; #singletons only have 1 cluster in them
				foreach my $mcl (nsort keys %{$clustergroups->{$group}->{$cluster}}) {
					my $mclann = "$mcl unknown";
					$mclann = "$mcl $mcl2ann->{$mcl}" if (exists $mcl2ann->{$mcl});
					$mclsummary{$mcl}++;
				}
				foreach my $mcl (keys %mclsummary) {
					$groupsummary{"${smclass}_singleton${singlecount}_$clustersize"}{$mcl} = $mclsummary{$mcl}; #copy mcl counts over to group summary once the sm string has been found
				}
				$omeclustergroups{$ome}{"${smclass}_singleton${singlecount}_$clustersize"}++;
			}
		}
	}
	return(\%omeclustergroups, \%allgroups, \%groupsummary, \%group2size);
}

sub Group_clusters_by_content {
	my ($clusters, $referenceclusters, $MINSIM) = @_;
	my (%groups, %sizerefs, %countedclusters, %allclusters);
	my $groupcount = 0;
	foreach my $querysize (sort {$b <=> $a} keys %{$clusters}) { #start from largest clusters to smallest
		foreach my $querycluster (keys %{$referenceclusters->{$querysize}}) { #notice the iteration through reference clusters
			next if (exists $countedclusters{$querycluster});
			print "*" if (scalar(keys(%countedclusters)) % 100 == 0);
 			$allclusters{$querycluster} = 1; #keep track of ALL clusters in analysis
			$countedclusters{$querycluster} = 1;
			$groupcount++;
			my $groupname = "group${groupcount}_s$querysize";
			my (%querymcls);
			foreach my $querymcl (keys %{$clusters->{$querysize}->{$querycluster}}) {
				$groups{$groupname}{$querycluster}{$querymcl} = $clusters->{$querysize}->{$querycluster}->{$querymcl};
				#$querymcls{$querymcl} = $clusters->{$querysize}->{$querycluster}->{$querymcl}; #each mcl appears at least 1 number of times.. CHANGE HERE to allow multiple copies of a given mcl 
				$querymcls{$querymcl} = 1; 
			}
			foreach my $targetsize (keys %{$clusters}) {
				my $minsize = int($MINSIM * $querysize);
				my $maxsize = $querysize + ($querysize - $minsize); 
				if ($targetsize <= $maxsize && ($targetsize >= $minsize)) { #if the number of matches is greater than or equal to X% of the size of query cluster, rounded down
 					foreach my $targetcluster (keys %{$clusters->{$targetsize}}) {
 						$allclusters{$targetcluster} = 1; #keep track of ALL clusters in analysis
 						my $matches = 0;
 						foreach my $targetmcl (keys %{$clusters->{$targetsize}->{$targetcluster}}) {
							if (exists $querymcls{$targetmcl}) {
								$matches++;
								#always increment matches by the maximum number of possible matches
								#CHANGE HERE to allow multiple copies of a given mcl
# 								if ($querymcls{$targetmcl} >= $clusters->{$targetsize}->{$targetcluster}->{$targetmcl}) {
# 									$matches += $clusters->{$targetsize}->{$targetcluster}->{$targetmcl};
# 								} elsif ($querymcls{$targetmcl} < $clusters->{$targetsize}->{$targetcluster}->{$targetmcl}) {
# 									$matches += $querymcls{$targetmcl};
# 								}
							}
						}
						if (($matches >= int($MINSIM * $querysize))) { #if the number of matches is greater than or equal to X% of the size of query cluster, rounded down, which allows for a maximum of 1 gene missing for clusters size 2-10; 2 genes missing for clusters size 11-20, 3 genes missing for clusters size 21-30 etc
							next if (exists $countedclusters{$targetcluster});
							$groups{$groupname}{$targetcluster} = $clusters->{$targetsize}->{$targetcluster}; #copy whole hash over, including mcls and their counts
 							print "*" if (scalar(keys(%countedclusters)) % 100 == 0);
							$countedclusters{$targetcluster} = 1;
 						}
 					}
 				}
 			}
 			#Score all clusters based on their pairwise similarity to all other clusters in group.. Retrieve the reference cluster
 			#allow self comparisons to avoid errors when groups have a single cluster
 			my %clusterscores;
 			foreach my $cluster1 (keys %{$groups{$groupname}}) {
 				my %mcls1 = %{$groups{$groupname}{$cluster1}};
 				my $size1 = scalar keys %mcls1;
 				#my $size1 = Get_total_MCL_count(\%mcls1); #CHANGE HERE to allow multiple copies of a given mcl
 				my $denominator = $size1;
 				foreach my $cluster2 (keys %{$groups{$groupname}}) {
					my %mcls2 = %{$groups{$groupname}{$cluster2}};
 					my $size2 = scalar keys %{$groups{$groupname}{$cluster2}};
 					#my $size2 = Get_total_MCL_count(\%mcls2); #CHANGE HERE to allow multiple copies of a given mcl
 					$denominator = $size2 if ($size2>$size1);
 					#my ($matches) = Get_all_MCL_matches(\%mcls1, \%mcls2); #CHANGE HERE to allow multiple copies of a given mcl
 					my ($matches) = Get_all_unique_MCL_matches(\%mcls1, \%mcls2);
 					my $similarity = sprintf("%.3f", $matches / $denominator);
 					$clusterscores{$cluster1} += $similarity; #add all the similarity scores together
 				}
 			}
 			#make the highest scoring cluster in this group a reference for clusters of this size
 			foreach my $cluster1 (sort {$clusterscores{$b} <=> $clusterscores{$a}} keys %clusterscores) {
 				$sizerefs{$querysize}{$cluster1} = 1;
 				last;
 			}
 					
 		}
 	}
	#if a cluster present in the combined file was not assigned to a group, assign it to its own group (this is only a problem with clusters size 2 (ie only 2 unique mcls), sometimes they get skipped during group assignment)
 	foreach my $combinedcluster (keys %allclusters) {
 		if (not exists $countedclusters{$combinedcluster}) { 
 			$groupcount++;
 			my $querysize = 2;
 			my $groupname = "group${groupcount}_s$querysize";
			my (%querymcls);
			foreach my $querymcl (keys %{$clusters->{$querysize}->{$combinedcluster}}) {
				$groups{$groupname}{$combinedcluster}{$querymcl} = $clusters->{$querysize}->{$combinedcluster}->{$querymcl};
			}
		}
	}				
 	print "\n";
	return(\%groups, \%sizerefs);
}

sub Get_all_unique_MCL_matches { 
	my ($mcls1, $mcls2) = @_;
	my $matches = 0;
	foreach my $mcl1 (keys %{$mcls1}) {
		foreach my $mcl2 (keys %{$mcls2}) {
			if (exists $mcls1->{$mcl2}) {
				$matches++;
			} 
		}
	}
	return($matches);
}

sub Get_all_MCL_matches {
	my ($mcls1, $mcls2) = @_;
	my $matches = 0;
	foreach my $mcl1 (keys %{$mcls1}) {
		foreach my $mcl2 (keys %{$mcls2}) {
			if (exists $mcls1->{$mcl2}) {
				if ($mcls1->{$mcl1} >= $mcls2->{$mcl2}) { #add the maximum possible number of matches (dictated by the lowest mcl count)
					$matches += $mcls2->{$mcl2};
				} else {
					$matches += $mcls1->{$mcl1};
				}			
			} 
		}
	}
	return($matches);
}


sub Get_total_MCL_count {
	my ($mcls) = @_;
	my $size = 0;
	foreach my $mcl (keys %{$mcls}) {
		$size += $mcls->{$mcl};
	}
	return($size);
}

sub Combined_parser {
	#my ($newcombinedfile, $restrictedmcls, $omesoi) = @_; #implementation with restrictedmcls, redundant since clusters were made with restricted mcls
	my ($newcombinedfile, $omesoi) = @_;
	my (%cluster2clusterclass, %clusters);
	my (%clusteredmcls, %tempclusters); #internal use only
	#first find all mcls present in >1 cluster
	open(my $combin, '<', $newcombinedfile) or usage("Error: could not open $newcombinedfile for reading\n");
	while(my $line = <$combin>) {
		next if ($line =~ m/^#/);
		chomp $line;
		my ($ome, $clusterid, $protid, $SMclass, $contig, $coords, $strand, $mcl) = split/\t/, $line;
		if (not exists $omesoi->{$ome}){
			warn("$ome is not in the restricted ome file\n");
			next;
		}
		#next if (not exists $restrictedmcls->{$mcl}); #restrict clustered mcls to those in restricted hash, redundant, since these clusters were made using only mcls in restricted hash 
		$cluster2clusterclass{$clusterid}{$SMclass} = 1 if ($SMclass ne 'unknown'); #some clusters might have multiple SM classes assigned to them if smurf predicted clusters were joined in the newSMURF prediction
		$clusteredmcls{$mcl}++;
	}
	#now assign mcl content to clusters	
	open(my $combin2, '<', $newcombinedfile) or usage("Error: could not open $newcombinedfile for reading\n");
	while(my $line = <$combin2>) {
		next if ($line =~ m/^#/);
		chomp $line;
		my ($ome, $clusterid, $protid, $SMclass, $contig, $coords, $strand, $mcl) = split/\t/, $line;
		if (not exists $omesoi->{$ome}){
			warn("$ome is not in the restricted ome file\n");
			next;
		}	
		#next if (not exists $restrictedmcls->{$mcl}); #restrict clustered mcls to those in restricted hash 
		if ($clusteredmcls{$mcl} > 1) { 
			$tempclusters{$clusterid}{$mcl}++;
		}
	}
	#now sort clusters by size
	foreach my $cluster (keys %tempclusters) {
		my $clustersize = 0;
		foreach my $mcl (keys %{$tempclusters{$cluster}}) {
			#$clustersize += $tempclusters{$cluster}{$mcl}; #CHANGE HERE to allow multiple copies of a given mcl
			$clustersize++;
		}
		next if ($clustersize == 1); #don't include clusters of size 1
		foreach my $mcl (keys %{$tempclusters{$cluster}}) {
			$clusters{$clustersize}{$cluster}{$mcl} = $tempclusters{$cluster}{$mcl};
		}
	}
	print scalar(keys(%tempclusters))." clusters retrieved from $newcombinedfile for grouping analysis..\n";
	return(\%cluster2clusterclass, \%clusters);
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


sub Opts_check {
	my ($opts) = @_;
	usage("Error: please provide a minimum content similarity cutoff\n") if (not defined $opts->{'s'});
	usage("Error: please provide a minimum cluster size to print for cluster matrix\n") if (not defined $opts->{'c'});
	usage("Error: please provide a minimum clustergroup size to print for cluster matrix\n") if (not defined $opts->{'g'});
	usage("Error: minimum content similarity cutoff must be between 0 and 1\n") if ($opts->{'s'} < 0 || $opts->{'s'} > 1);
	usage("Error: please specify either a restricted mcl file or a combined file\n") if (! defined $opts->{'m'} && ! defined $opts->{'f'});
	usage("Error: please provide an output directory\n") if (not defined $opts->{'o'});
	usage("Error: cannot find the output directory\n") if (! -d $opts->{'o'});
	if (defined $opts->{'m'}) {
		usage("Error: cannot find the restricted mcl file\n") if (! -f $opts->{'m'});
	}
	if (defined $opts->{'f'}) {
		usage("Error: cannot find the combined file\n") if (! -f $opts->{'f'});
	}
}
