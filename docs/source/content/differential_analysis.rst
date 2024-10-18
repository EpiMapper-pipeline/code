=====================
differential_analysis
=====================

.. contents::
    :local:

Preforms differntial analysis  using the student´s t-test on enriched regions/peaks before annotating the stastitically significant changes to spesific genomic regions and visulizing the results. 

For ATAC-seq, fold enrichment is utilized since ATAC-seq experiments usually do not include control samples. Thus, the read counts from alignment are affected by sequencing depth, hence, the samples should not be compared directly. 

The annotation step is derived from another epigenetic Python Package: `HMST-seq-Analyzer <https://hmst-seq.github.io/hmst/>`_ and allows for annotation to gene regions, transcription start sites (TSS), transcription end sites (TES), and end 5´ distance regions. Further, annotation to enhancer regions may be done, however, a separate enchancer BED must provided in the ``-e, --enhancer`` parameter. These files must be processed to have this format:

.. code-block:: bash

    chr1	50218	50564	chr1:50218:50564||enhancer:2||346
    chr1	111057	111395	chr1:111057:111395||enhancer:1||338
    chr1	173031	174331	chr1:173031:174331||enhancer:3||1300
    chr1	174502	176296	chr1:174502:176296||enhancer:6||1794

Arguments
===========

**Required arguments:**

- ``-p, --peaks``: Input file folder of sorted BED files from peak calling.
- ``-cs, --chromosome_sizes``: Input file of sorted chromosome sizes information.
- ``-bl, --blacklist``: Input file with genomic reagion that should be excluded.
- ``-r, --refrence_refFlat``: RefFlat text (.txt) file to which the peaks will be annotated to. The reference file will be used for the division of genomic regions (i.e., Transcription start site/promotor, intrageneic etc.)".
- ``-la, --list_a``: List of sample names that will be used to do differential analysis.
- ``-lb, --list_b``: List of sample names that will be used to do differential analysis.

**Required if not fold enrichment is used**

- ``bg, --bedgraph``: Input file folder of BedGraph fragment files for differntial analysis. This parameter has no default value, and is required if the ``-fold, --fold_enrichment`` parameter is not used or "False". 


**Optional arguments:**

- ``-tm, --test_methods``: Four hypothesis testing methods are available for selection: t-test, Kolmogorov-Smirnov test (K-S test), Mann-Whitney U test, and Wilcoxon rank-sum test. The corresponding values are "ttest", "kstest", "mannwhitneyu", and "ranksumtest". Default=ttest.
- ``-an, --annotate``: Wheither or not to annotate the data. The annotation step might take a considerable amount of time, especially with large samples (i.e geneom wide ATAC-seq). Default=False 
- ``-fold, --fold_enrichment``: Input must be either "True" or "False". The function will use fold enrichemnt as well as normalization before differential analysis. Default=False
- ``-cut, --p_value_cutoff``:  Cut-off p-value for differential analysis, default =0.05
- ``-n , --normalize``: Whether or not to quantile normalize the data. Might be beneficial if spike in calibration is not preformed. Default = False.
- ``-X``: The number of upstream basepaires from TSS, TES, gene, when creating genomic region files. Default = 1000.
- ``-Y``: The number of downstream basepairs from TSS, TES, gene, when creating genomic region files. Default = 1000.
- ``-M``: The number of basepairs from gene start site, 5dist, when creating genomic region files. Default=10000.
- ``-N``: The number of basepairs from gene start site, 5dits, when creating genomic region files. Default=1000000.
- ``-l , --minIntergenicLen``: Minimim intergenic region distance. Default = 2000.
- ``-xL, --maxIntergenicLen``: Maximum intergentic region distance. Default = 10000000.
- ``-i, --intergenicBTGenes``:  Wheither intergenic regions is considrered between gene body regions (True), or betweeen TSS and TES (False). Default=True.
- ``-e, --enhancer``: Sorted BED file with defined enhancer regions for annotation.



Example Usage
================


Example usage default, compareing normal vs cancer samples, no annotation:

.. code-block:: bash

    epimapper differential_analysis -p /Users/me/results/Epimapper/peakCalling/seacr/control -cs /Users/me/in_folder/hg38_chromosome_sizes.txt -bg /Users/me/results/Epimapper/alignment/bedgraph -bl /Users/me/in_folder/hg38_blacklist.bed -r  /Users/me/in_folder/hg38_refFlat.txt -la H3K4me3-cancer_rep1 H3K4me3-cancer_rep2 -lb H3K4me3-healthy_rep1 H3K4me3-healthy_rep2


Example usage default, compareing normal vs cancer samples, with annotation:

.. code-block:: bash

    epimapper differential_analysis  -an True  -e /Users/me/in_folder/enhancer.bed -p /Users/me/results/Epimapper/peakCalling/seacr/control -cs /Users/me/in_folder/hg38_chromosome_sizes.txt -bg /Users/me/results/Epimapper/alignment/bedgraph -bl /Users/me/in_folder/hg38_blacklist.bed -r  /Users/me/in_folder/hg38_refFlat.txt -la H3K4me3-cancer_rep1 H3K4me3-cancer_rep2 -lb H3K4me3-healthy_rep1 H3K4me3-healthy_rep2


Example usage ATAC-seq: 

.. code-block:: bash

    epimapper differential_analysis  -an True  -e /Users/me/in_folder/enhancer.bed -p /Users/me/results/Epimapper/peakCalling/macs2/top_peaks -cs /Users/me/in_folder/hg38_chromosome_sizes.txt -bl /Users/me/in_folder/hg38_blacklist.bed -r  /Users/me/in_folder/hg38_refFlat.txt -la H3K4me3-cancer_rep1 H3K4me3-cancer_rep2 -lb H3K4me3-healthy_rep1 H3K4me3-healthy_rep2



Output
=======

Like all the other functions in EpiMapper Python package, the function will create a main ``Epimapper`` output directiry, if it is not already present in the chosen output directory.
Further, this function will create multiple directories and sub-folders to store all output files, seen below. However, the output directories look different if annotation as not done. 

Firstly a "differential_analysis" folder is created where blacklist 100 bp bin files of every samples  as well as three sub-folders: DAR, out_combined_file, region_files are as stored.

The "out_combined_files" subfolder contains several files and a "genome" folder. The "combined_peaks.bed" contains all peaks in all samples, while "combined_peaks_merged.bed" contains the same peaks merged based on overlapping intervals.
The "mapped_peaks_100bp_merged.bed.gz" contains peaks mapped to the blacklist 100 bp bin files. 
In the "genome" folder, files containing annoted differntially enriched peaks are stored, and "combined_signals_100b.bed.gz" contians every samples signal, it is these two files that are used as input for the t-test. "combined_signals_100b.head" is the column names of "combined_signals_100b.bed.gz".

In the "DAR" folder, the output from the differntial analysis t-test is stored, with "combined_peaks_merged_pval.csv" containing all peaks and their p-values, "combined_peaks_merged_pvals0.01.csv" containing only the peaks with a p-value < 0.01, and "combined_peaks_merged_pvals0.01_genome.csv" containing the annotated p-vaule < 0.01 peaks.

The "region_files" folder contains the region files which are used for the annotation. 

Lastly, this function also creates several plots: total_AR_pie.pdf - containing a pie plot of where the differntial peaks are annotated, "groupA-vs_groupb"_pca.pdf - containing a 3D PCA plot of the differnt samples, and Epimapper_DAR_ttest_pval_0.01.pdf containing a bar plot explaining up and down enrichment as well as where in the genome the peaks are annotated.

.. code-block:: bash

    Epimapper
    |- differential_analysis
    |   |- "sample-name"_100b.windows.BlackListFiltered.bed.gz
    |   |- out_combined_files
    |   |   |    |- combined_peaks.bed 
    |   |   |    |- combined_peaks_merged.bed 
    |   |   |    |- combined_peaks_merged_pval_min_0.01.bed 
    |   |   |    |- combined_signals_100b.head
    |   |   |    |- mapped_peaks_100bp_merged.bed.gz
    |   |   |    |- combined_signals_100b.bed.gz 
    |   |   |    |- genome
    |   |   |    |    |- combined_peaks_merged_pval_min_0.01_gene_Up1000_Down1000removedShort_overlap1e-09.bed
    |   |   |    |    |- combined_peaks_merged_pval_min_0.01_intergenic_uniqueSorted_betweenGenes_minLen2000_overlap1e-09.bed
    |   |   |    |    |- combined_peaks_merged_pval_min_0.01_noGenes_5dist_Down1000000_Down10000removedShort_overlap1e-09.bed
    |   |   |    |    |- combined_peaks_merged_pval_min_0.01_noGenes_5dist_Up1000000_Up10000removedShort_overlap1e-09.bed
    |   |   |    |    |- combined_peaks_merged_pval_min_0.01_TES_Up1000_Down1000removedShort_overlap1e-09.bed
    |   |   |    |    |- combined_peaks_merged_pval_min_0.01_TSS_Up1000_Down1000_removedShort_overlap1e-09.bed
    |   |   |    |    |- Epimapper_DAR_ttest_pval_0.01.csv              
    |   |- DAR
    |   |   |- combined_peaks_merged_pval.csv
    |   |   |- combined_peaks_merged_pvals0.01.csv
    |   |   |- combined_peaks_merged_pvals0.01_genome.csv
    |   |- region_files
            |- 5dist_Down1000000_Down10000removedShort.bed
            |- hg19.refFlat_clean_sorted.txt
            |- TES_Up1000_Down1000removedShort.bed
            |- 5dist_Up1000000_Up10000removedShort.bed
            |- intergenic_uniqueSorted_betweenGenes_minLen2000.bed
            |- TSS_Up1000_Down1000_removedShort.bed
            |- gene_Up1000_Down1000removedShort.bed
            |- list_region_files.txt
            |- refFlat_clean_sorted.bed
            |- removed_regions_all_TSS_TES_5dist_geneBodyLessThan0.bed
    |- summary_tables
    |   |- total_AR_pie.pdf
    |   |- "groupA-vs_groupb"_pca.pdf 
    |   |- Epimapper_DAR_ttest_pval_0.01.pdf


