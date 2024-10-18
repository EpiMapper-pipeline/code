=================
peak_calling
=================

.. contents::
    :local:

Finds enriched regions/calls for peaks from chromatin profiling data with `SEACR <https://github.com/FredHutch/SEACR>`_ or `MACS2 <https://github.com/macs3-project/MACS>`_, then visulizes results.

Both of these peak calling software have an alternative where control samples are used to create a background noise level before finding enriched reagions, peaks. However, not every experiment does include control samples, therefore this a optional feature in the EpiMapper pipeline. If your experiment includes control samples and you wish it perform peak calling based on these, you may utilize the ``-c, --control_index`` parameter to insert the names of the control sampels.
The peak_calling function will separate all the samples containing the ``-c, --control_index``, and use these at input for background-noise. 

Arguments 
===========


**Required arguments:**

- ``-soft, --software``: input has to be either "macs2" or "seacr". Here you decide which peak calling software you would like to use. Generally, MACS2 is more used and prefered for samples with higher background noise ( ie. ATAC-seq and CHiP-seq), while SEACR may be used for samples with lower noise.

- ``-f, --fragments``: Input file folder of filterd BED files for peak calling.

*If peak calling with SEACR software:*

- ``-bg, --bedgraph``: Input file folder of BedGraph files for peak calling with SEACR.
- ``-norm, --seacr_norm``: Seacr peakcalling normalization option (norm, non), default= non for Seacr peakcalling", default="non".

*If peak calling with MACS2 software:*

- ``-b --bam``: Input folder containing BAM files for MACS2 peak calling.
- ``-gs, --genome_size``: Relative genome size of organism being studied.  About 90% or 70% of the genome size. (i.e, 2.7e9 for humans, 1.87e9 for mice, 9e7 for *Caenorhabditis elegans*, or 1.2e8 for fruitfly).

*If peak calling with with control samples:*

- ``-la, --list_a``: A list of sample names.
- ``-lb, --list_b``: A list of control sample names.

.. note::

    The sample names and control sample names must be entered and they should be in the correct paired order in the list.

- ``-c, --control_index``: Indexes of control files (i.e, 'control' 'IgG' ect.) that should be used as a background-noise reference during peak calling. If not provided the peak calling softwares will select peaks based on a cut-off value (i.e p-value for SEACR). Default = False. An example would be if you have the samples: [H3K4me3_rep1, H3K4me3_rep2, H3K27me3_rep1, H3K27me3_rep2, IgG_rep1, IgG_rep2], you would input "IgG" into the ``-c, --control_index``. If you have the samples [healthy_rep1, healthy_rep2, cancer_rep1, cancer_rep2, control_rep1, control_rep2] you would input "control" into the ``-c, --control_index``. 



**Optional arguments:**

- ``-p, --percentage``: Cut-off percentage for peak calling software, default= 0.01.
- ``-qval, --macs2_qvalue``: Macs2 callpeak Q-value, default is None, if qvalue is used then Percentage or P-value will not be considered, default= None.
- ``-eB, --export_bdg``: Macs2 Whether or not to save extended fragment pileup, defualt=False for not export, use True or exporting", default= False.
- ``-tbl, --fragment_table``: Input CSV file with the following columns = ["Sample", "Replication", "SequencingDepth", "MappedFragments", "AlignmentRate"] with corresponding sample information, default = “bowtie2_alignment_ref.csv” exported by this pipeline function: ``bowtie2_alignment``.
- ``-o, out_dir``: Output directory, default = current working directory.

.. note::

    If you first perform peak calling without control samples and then want to perform peak calling with control samples, you need to use a new output folder.


Example Usage:
===============

Example usage SEACR peak calling with control samples and specified cut-off percentage:

.. code-block:: bash
    
    $ epimapper peak_calling -soft seacr -bg /Users/me/results/Epimapper/alignment/bedgraph -f /Users/me/results/Epimapper/alignment/bed -la H3K27me3_rep1 H3K27me3_rep2 H3K4me3_rep1 H3K4me3_rep2 -lb IgG_rep1 IgG_rep2 IgG_rep1 IgG_rep2 -c IgG -p 0.05 -o /Users/me/results



Example usage SEACR peak calling without control and default cut-off precentage:

.. code-block:: bash
    
    $ epimapper peak_calling -soft seacr -bg /Users/me/results/Epimapper/alignment/bedgraph -f /Users/me/results/Epimapper/alignment/bed  -o /Users/me/results 


Example usage with MACS2 and control samples and normalized data:

.. code-block:: bash

    $ epimapper peak_calling -soft macs2 -b /Users/me/results/Epimapper/alignment/bam -gs 1.87e9 -norm norm -f /Users/me/results/Epimapper/alignment/bed -la ZNF143-48h_rep1 ZNF143-48h_rep2 ZNF143-72h_rep1 ZNF143-72h_rep2 -lb ZNF143-Control-48h_rep1 ZNF143-Control-48h_rep2 ZNF143-Control-72h_rep1 ZNF143-Control-72h_rep2 -c Control -o /Users/me/results 


Example usage with MACS2, no control samples:

.. code-block:: bash
    
    $ epimapper peak_calling -soft macs2 -b /Users/me/results/Epimapper/alignment/bam -gs 1.87e9 -f /Users/me/results/Epimapper/alignment/bed -o /Users/me/results 



Output
========

Like all the other functions in EpiMapper Python package, the function will create a main ``Epimapper`` output directiry, if it is not already present in the chosen output directory. Further, this function will create multiple folders for each file conversion.
Depending on your prefered peak calling software as well as if control samples are used, the output directory will differ. For detailed information about all output files created by `SEACR <https://github.com/FredHutch/SEACR>`_ and `MACS2 <https://github.com/macs3-project/MACS>`_ please visit their websites. 
Further, this function will create a summary table (peak_summary.csv) as well as several plots of various peak statisitcs: peak width violin plot, peak number, peak reproducibility and fragment proposition in peaks (FRiPS). 


.. code-block:: bash

    Epimapper
    |- peakCalling
    |   |- macs2
    |   |   |- top_peaks
    |   |   |    |- "sample-name"_macs2_top._peaks.narrowPeak 
    |   |   |    |- "sample-name"_macs2_top._peaks_sorted.bed   
    |   |   |    |- "sample-name"_macs2_top._peaks.xls 
    |   |   |    |- "sample-name"_macs2_top._summits.bed
    |   |   |- control_peaks
    |   |   |    |- "sample-name"_macs2_control_peaks.stringent.bed
    |   |- seacr
    |   |   |- top_0.01 
    |   |   |    |- "sample-name"_seacr_top.0.01_peaks.stringent.bed
    |   |   |- control_peaks
    |   |   |    |- "sample-name"_seacr_control_peaks.stringent.bed
    |- summary_tables
    |   |- peak_summary.csv
    |   |- Peak_width.png
    |   |- Peak_numbers.png
    |   |- frips.png
    |   |- peaks_reproducibility_rate.png




