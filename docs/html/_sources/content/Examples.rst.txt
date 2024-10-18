================================
Examples - Ready-made Pipeline
================================

.. contents::
    :local:


For novice command-line users we have made this page with examples on usage, as well as a ready-made file for useing the complete pipeline.



Ready-made Files:
==================

Here is a shell script running the entire pipeline :download:`full pipeline example script <_static/example_full_pipeline.sh>`

Here is the configuration script needed to run the shell script of the full pipeline :download:`configuration script<_static/config.txt>`



Set Up Configuration File:
==========================

Before using the ready-made file for useing the complete pipeline, you need to set up the configuration file with all the paths to your files.

We recommend to create a new folder containing everything you need to run the pipeline, so everything is in one place and is easier to find on a later occasion. 

Setting up the configuration file itself is quite easy as the file itself contains a lot of comments explaining where you should provide the various information.

Here is an example setup for an input folder with the directory structure showcased below:


.. code-block:: bash

    in_folder_cut_tag
    |- fastq/
    |   |- diabetic-1_rep1.fastq 
    |   |- diabetic-1_rep2.fastq 
    |   |- diabetic-2_rep1.fastq 
    |   |- diabetic-2_rep2.fastq 
    |   |- healthy-1_rep1.fastq 
    |   |- healthy-1_rep2.fastq 
    |   |- healthy-2_rep1.fastq 
    |   |- healthy-2_rep2.fastq 
    |   |- healthy-3_rep1.fastq 
    |   |- healthy-3_rep2.fastq 
    |   |- IgG-1_rep1.fastq 
    |   |- IgG-1_rep2.fastq 
    |   |- IgG-2_rep1.fastq 
    |   |- IgG-2_rep2.fastq 
    |- hg19_bowtie2_index/
    |- ecoli_bowtie2_index/
    |- hg19-blacklist_sorted.bed
    |- hg19_chromosome_sizes_sorted.txt
    |- hg19.refFlat.txt
    |- hg19_all_enhancers_merged_4dmr.bed 



The configuration file below will execute every function, including spike_in_calibration, and in this example the data is thought to not be ATAC-seq data. Here, the peak calling software is set to MACS2 using the IgG control samples as a background reference. 

.. code-block:: bash 

    out_dir="/Users/me/documents/cut_tag_results"

    fastq="/Users/me/documents/in_folder_cut_tag/fastq"

    bowtie2_index="/Users/me/documents/in_folder_cut_tag/hg19_bowtie2_index"

    bowtie2_index_spike_in="/Users/me/documents/in_folder_cut_tag/ecoli_bowtie2_index"

    blacklist="/Users/me/documents/in_folder_cut_tag/hg19-blacklist_sorted.bed"

    atac="False"

    chromosome_sizes="/Users/me/documents/in_folder_cut_tag/hg19_chromosome_sizes_sorted.txt"

    software="macs2"

    genome_size="2.7e9"

    control_index="IgG"

    refflat_ref="/Users/me/documents/in_folder_cut_tag/hg19.refFlat.txt"

    group_A="diabetic-1_rep1 diabetic-1_rep2 diabetic-2_rep1 diabetic-2_rep2"

    group_B="healthy-1_rep1 healthy1-rep2 healthy-2_rep1 healthy-2_rep2 healthy-3_rep1 healthy_rep2"




Running The Ready-made File:
================================

When you have set up your configuration file, you are ready to run the analysis! You do not have to change anything in the example_full_pipeline.sh shell script, unless you wish to exclude some of the function.
Just remember that most of the functions do rely on each other. 

When running the ready-made shell script containing the pipeline, it is important to place the config.txt file and the shell script in the same folder. If you they are not in the same folder, the shell script will not find the configuration file, leading to faliure. 

Then, when you are ready you may open your terminal window, move to the folder you have placed the config.txt file and shell script files. Then run the following command to start the analysis:

.. code-block:: bash

   $ bash example_full_pipeline.sh


