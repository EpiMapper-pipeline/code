=====================
spike_in_calibration
=====================

.. contents::
    :local:

Removes experimental bias by normalizing fragment counts based on sequencing depth to a spike-in genome and visulizes results.

Arguments
==========

**Required arguments:**

- ``-b, --bed``: Input file folder of filterd BED files for normalization
- ``-ss, --sam_spike_in``: Input file folder of SAM files exported from alignment to a spike in genome.
- ``-cs, --chromosome_sizes``: Input file of sorted chromosome sizes information.


**Optional arguments:**

- ``-tbl, --fragment_table ALIGNMENT SUMMARY TABLE``:Input CSV file containing the following columns = ["Sample",	"Replication", "SequencingDepth", "MappedFragments", "AlignmentRate", "MappedFragments_SpikeIn",	"AlignmentRate_SpikeIn"] with corresponding sample information , default = “bowtie2_alignment_ref_and_spike_in.csv” exported by this pipeline function: ``bowtie2_alignment``.
- ``-o, --out_dir``: Output directory, default = current working directory.

Example usage
==============

The function will assume that the “bowtie2_alignment_ref_and_spike_in.csv" file is present $out_dir/"Epimapper/summary_tables". Therefore, it is important to use the same output directory "-o/--out_dir" as the one you utilized for the spike-in alignment. This will make sure that the fucntion will find the table.

.. code-block:: bash
    
    $ epimapper spike_in_calibration -b /Users/me/results/Epimapper/alignment/bed -ss /Users/me/results/Epimapper/alignment/sam_spike_in -cs /Users/me/in_folder/hg38_chromosome_sizes.txt -o /Users/me/results

If you want a differnet output directory you may choose to input the path to the table manually:

.. code-block:: bash
    
    $ epimapper spike_in_calibration -b /Users/me/results/Epimapper/alignment/bed -tbl /Users/me/results/Epimapper/summary_tables/bowtie2_alignment_ref_and_spike_in.csv -ss /Users/me/results/Epimapper/alignment/sam_spike_in -cs /Users/me/in_folder/hg38_chromosome_sizes.txt -o /Users/me/results 


If you have not used this pipelines ``bowtie2_alignment`` to preform the reference genome and spike-in alignment, you must manually create a summary table containing the following columns:["Sample", "Replication", "SequencingDepth", "MappedFragments", "AlignmentRate", "MappedFragments_SpikeIn",	"AlignmentRate_SpikeIn"] with corresponding infromation for each sample. 
Therefore, it is recommended to use the pipeline as a whole to avoid any manual labor.


.. code-block:: bash
    
    $ epimapper spike_in_calibration -b /Users/me/results/Epimapper/alignment/bed -tbl /Users/me/results/my_table.csv -ss /Users/me/results/sam_spike_in -cs /Users/me/in_folder/hg38_chromosome_sizes.txt -o /Users/me/results 


Output
=======

Like all the other functions in EpiMapper Python package, the function will create a main ``Epimapper`` output directiry, if it is not already present in the chosen output directory. Further, this function will create a "bedgraph" folder to store the spike-in normalized files. Further, this function will create a summay table and a PNG figure with boxplots of spike-in scaling factors and normalized fragment count.

.. code-block:: bash

    Epimapper
    |- alignment
    |   |- bedgraph
    |   |   |- "sample-name".fragments.normalized.bedgraph
    |- summary_tables
    |   |- spike_in_calibration_summary.csv 
    |   |- spike_in_calibration.png



