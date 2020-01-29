#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use File::Basename;
use Data::Dumper;
use FileHandle;
use Sort::Naturally;
use Ebinf::Utils qw(dim_1_hash);
$| = 1;

#ome is defined by everything up to the last _
#backbone genes should be identified and highlighted in this analysis
#also, combined.txt will eventually have a column with eggnog annotations
sub usage {
	my $message = shift;
	$message = '' if (not defined $message);
	my $usage = qq/
usage: $0 [OPTIONS]
-c <combined_file>
-a <dotSM_emapper_annotation.txt>
-o <output_directory>
This script calculates the frequencies of emapper annotations assigned to 
proteins in each homolog group (ortholog group: OG) found in a SMURF cluster\n
Main output: dotSM_OG_annotation.txt\n\n/;
	die($usage, $message);
}
 
main: {
	my %opts;
	getopt('a:c:o:h', \%opts);
	Opts_check(\%opts);
	my ($emapperfile, $combinedfile, $outdir) = ($opts{'a'}, $opts{'c'}, $opts{'o'});
	open(my $annout, '>', "$outdir/dotSM_OG_annotation.txt") or usage("Error: could not open $outdir/dotSM_OG_annotation.txt for writing\n");
	$annout->autoflush(1);
	open(my $processout, '>', "$outdir/dotSM_OG_process.txt") or usage("Error: could not open $outdir/dotSM_OG_process.txt for writing\n");
	$processout->autoflush(1);
	
#1. Find each prot's OG from combined file
	my ($prot2og) = Prot2OG($combinedfile);
	
#2. Retrieve each prot's annotation and process(es)
	my ($prot2ann, $prot2process) = Retrieve_prot_annotation($emapperfile);
	
#3. Assess frequency of annotations within each OG. output structued {og}{ann} = freq
	my ($OGprocessfreq) = Assess_OG_annotation_freq($prot2process, $prot2og);
	my ($OGannfreq) = Assess_OG_annotation_freq($prot2ann, $prot2og);
#4. Print out report	
	Print_annotation_freq($OGprocessfreq, $processout);
	Print_annotation_freq($OGannfreq, $annout);
}

sub Print_annotation_freq {
	my ($annfreq, $out) = @_;
	foreach my $og (nsort keys %{$annfreq}) {
		print $out "$og\t";
		foreach my $ann (sort { $annfreq->{$og}->{$b} <=> $annfreq->{$og}->{$a} } keys %{$annfreq->{$og}}) {
			print $out "$ann ($annfreq->{$og}->{$ann}),";
		}
		print $out "\n";
	}
}


sub Assess_OG_annotation_freq {
	my ($prot2ann, $prot2og) = @_;
	my (%OGfreq, %protfreq);
	foreach my $prot (keys %{$prot2ann}) {
		$protfreq{$prot2og->{$prot}}{$prot} = 1;
		foreach my $ann (@{$prot2ann->{$prot}}) {
			my @fields = split/,/, $ann; #someitmes there can be multiple comma separated annotations, (e.g., for processes)
			foreach my $field (@fields) {
				$field = lc $field; #some annotations are annoyingly not controlled by case (e.g., you can have chitin and Chitin)
				$field =~ s/^\s+//;
				$OGfreq{$prot2og->{$prot}}{$field}++;
			}
		}
	}
	my (%annfreq);
	foreach my $og (keys %OGfreq) {
		my $protcount = scalar keys %{$protfreq{$og}}; #how many prots in this OG?
		foreach my $ann (keys %{$OGfreq{$og}}) {
			my $anncount = $OGfreq{$og}{$ann};
			my $freq = sprintf("%.3f", $anncount/$protcount); #divide the number of annotations by the number of prots in OG
			$annfreq{$og}{$ann} = $freq;
		}
	}
	return(\%annfreq);
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



sub Prot2OG {
	my ($combinedfile) = @_;
	my (%prot2og);
	open(my $combin, '<', $combinedfile) or usage("Error: could not open $combinedfile for reading\n");
	while(my $line = <$combin>) {
		chomp $line;
		my ($ome, $clusterid, $protid, $SMclass, $contig, $coords, $strand, $og) = split/\t/, $line;
		$prot2og{$protid} = $og;
	}
	return(\%prot2og);
}

sub Opts_check {
	my ($opts) = @_;
	usage() if (exists $opts->{'h'});
	usage("Error: please provide the combined.txt file\n") if (not defined $opts->{'c'});
	usage("Error: please provide the emapper annotation file \n") if (not defined $opts->{'a'});
	usage("Error: please specify an output directory\n") if (not defined $opts->{'o'});
	usage("Error: cannot find the combined.txt file.. is it in the right place?\n") if (! -f $opts->{'c'});
	usage("Error: cannot find the emapper annotation file.. is it in the right place?\n") if (! -f $opts->{'a'});
	usage("Error: cannot find the output directory.. is it in the right place?\n") if (! -d $opts->{'o'});
}

