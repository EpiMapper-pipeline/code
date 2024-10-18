=================
filtering
=================

.. contents::
    :local:

Performs data filtering for mapped reads based on their alignment quality, and file format conversion before high-level data analysis,before  visulizing reproducibility among biological replicates.

Arguments
==========

**Required arguments:**

- ``-s, --sam``: Input file folder of SAM files exported from alignment for filtering.
- ``-cs, --chromosome_sizes``: Input file of sorted chromosome sizes information.
- ``-bl, --blacklist``: Input file in BED format with genomic region that should be excluded.

**Optional arguments:**

- ``-o, out_dir``: Output directory, default = current working directory.
- ``-sn, --spike_in_norm``  Input must be either "True" or "False". If the files are being spike-in normalized after filtering or not. If set to "True" the function will not create BedGraph files, as these are created by the ``spike_in_calibration`` function. Default=False
- ``-atac, --atac_seq_shit`` Input must be either "True" or "False". Shift start point in ATAC-seq data +4 on the forward (+) strand and -5 on the backwards (-) strand. Default=False

Example Usage
==============

Example usage for samples that should not be shifted:

.. code-block:: bash

    $ epimapper filtering -s /Users/me/results/Epimapper/alignment/removeDuplicate/sam_duplicates_removed -cs /Users/me/in_folder/hg38_chromosome_sizes.txt -bl /Users/me/in_folder/hg38_blacklist.bed -o /Users/me/results


Example usage for ATAC-seq data that should be shiftet:

.. code-block:: bash

    $ epimapper filtering -s /Users/me/results/Epimapper/alignment/removeDuplicate/sam_duplicates_removed -atac True -cs /Users/me/in_folder/hg38_chromosome_sizes.txt -bl /Users/me/in_folder/hg38_blacklist.bed -o /Users/me/results

Example usage for data that should be spike-in normalized after filtering:

.. code-block:: bash

    $ epimapper filtering -s /Users/me/results/Epimapper/alignment/removeDuplicate/sam_duplicates_removed -sn True -cs /Users/me/in_folder/hg38_chromosome_sizes.txt -bl /Users/me/in_folder/hg38_blacklist.bed -o /Users/me/results



Output
=======

Like all the other functions in EpiMapper Python package, the function will create a main ``Epimapper`` output directiry, if it is not already present in the chosen output directory. Further, this function will create multiple folders for each file conversion.
Natrually, the "bam" folder will contain BAM files, both unedited as well as a mapped and sorted version. The "bed" folder will contain BED files, which are sorted and filtered on alignment quality score framgent length (< 1000 bp) and blacklist removed. Here, the blacklist removed files are in 500 pb bin form, meaning that the BED file has this form [chromosome    start   end     count], where each "start" and "end" are a 500bp intervall and "counts" contains the number of fragments which were present in this intervall.
The "bedgraph" folder contains the BedGraph version of the sorted bed files. Lastly, the function will also create two JPG files contraining correlation heatmap plots. These plots are based on the log-transformed counts in blacklist removed 500 bin files, and the "corrcoef_heatmap4logCount_all.jpg" includes every bin entry, while "corrcoef_heatmap4logCount_filtered_gt8.jpg" only accepts bin entries with a total count of 8 across all samples.

.. code-block:: bash

    Epimapper
    |- alignment
    |   |- bam
    |   |   |- "sample-name".bam
    |   |   |- "sample-name".mapped_sorted.bam
    |   |- bed
    |   |   |- "sample-name".fragments.bed  
    |   |   |- "sample-name".fragments_sorted.bed
    |   |   |- "sample-name".fragments_sorted.500b.windows.BlackListFiltered.bed
    |   |- bedgraph
    |   |   |- "sample-name".fragments.bedgraph 
    |- summary_tables
    |   |- corrcoef_heatmap4logCount_all.jpg 
    |   |- corrcoef_heatmap4logCount_filtered_gt8.jpg

