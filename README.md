# Identifying gene clusters using gene co-occurrences
This repository contains all of the necessary scripts for implementing the CO-OCCUR gene cluster detection pipeline, as described in:

Gluck-Thaler E, Haridas S, Binder M, Grigoriev I, Crous PW, Spatafora JW, Bushley K, Slot JC. The architecture of metabolism maximizes biosynthetic diversity in the largest class of fungi. Submitted.

# Pipeline implementation
Scripts are sequentially numbered according to the order in which they are executed. Specific descriptions are found within each script, but briefly:

| Script name | Description |
| --- | --- |
| dotSM0.sample_null_model.pl | Samples a random distribution of homolog group co-occurrences |
| dotSM1.cluster_all_SM_prots.pl | Extends predicted gene clusters to include all homolog groups found across entire cluster set |
| dotSM2.HG_annotation.pl | Calculates frequency of predicted annotations across members of each homolog group |
| dotSM3.smash2combined.pl | Converts antiSMASH v4 output files to a custom "combined" file format |
| dotSM3.smurf2combined.pl | Converts SMURF output files to a custom "combined" file format |
| dotSM4.test_all_co-occurrences.pl | Empirically estimates probability of observing co-occurrences of interest using null model |
| dotSM5.group_clusters_by_content.pl | Groups predicted clusters together based on minimum threshold of gene content similarity |
| dotSM6.cross_reference_clusters.pl | Utility script for comparing the output of various cluster detection algorithms |
| dotSM7.summary_statistics.pl | Utility script for extracting summary statistics of interest for graphical plotting |
