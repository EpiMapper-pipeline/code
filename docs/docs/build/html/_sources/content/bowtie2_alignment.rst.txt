=================
bowtie2_alignment
=================

.. contents::
    :local:

Arguments 
===========

Mapping reads to a reference genome with `Bowtie2 <https://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ alignment of FASTQ sequencing reads files from high-thoughput sequecing, and visulizing results.

**Required arguments:**

- ``-f, --fastq``: Path to directory containing FASTQ files.

*either*

- ``-i, --bowtie2_index_pathway``: Input file folder of Bowtie2 reference genome indexing files. 

*or*

- ``-r, --reference``: Input reference genome FASTA file for the creation of Bowtie2 reference genome. 

**Optional arguments:**

- ``-o, out_dir``: Output directory, default = current working directory.
- ``-s , --spike_in``: If the alignment is spike-in, default = False.
- ``-m , --merge``: Merges technical replicates, default = False.

Example Usage:
===============

Example with existing bowtie2 indexing files:

.. code-block:: bash

    $ epimapper bowtie2_alignment -f /Users/me/my_fastq_files -i /Users/me/in_folder/bowtie2_index_hg38 -o /Users/me/results


Example with reference genome FASTA file:

.. code-block:: bash

    $ epimapper bowtie2_alignment -f /Users/me/my_fastq_files -r /Users/me/in_folder/hg38.fna -o /Users/me/results


Example of spike-in alignment:

.. code-block:: bash

    $ epimapper bowtie2_alignment -f /Users/me/my_fastq_files -s True -i /Users/me/in_folder/bowtie2_index_ecoli -o /Users/me/results



Example usage with technical replicates:

.. code-block:: bash

    $ epimapper bowtie2_alignment -f /Users/me/my_fastq_files -m True -i /Users/me/in_folder/bowtie2_index_hg38 -o /Users/me/results

Output
========

Like all the other functions in EpiMapper Python package, the function will create a main ``Epimapper`` output directiry, if it is not already present in the chosen output directory.
Further, this function will create directories based on if the alignment is considered a "spike-in" alignment (if the ``-s`` parameter is set to True) or not. If not the output directroy structure will as shown below.
Here, the "sam" directory will contain the aligned sam files, and the "bowtie2_summary" folder will contian text files for each sample with alignment statistics. These text files are used to create the "bowtie2_alignment_rep.csv" located in the "summary_tables" folder. Lastly, a PNG image file containing boxplots of sequencing depths, mapped fragments and alignment rate is created and also stored in the "summary_tabels" folder.

.. code-block:: bash

    Epimapper
    |- alignment/
    |   |- sam
    |   |-  bowtie2_summary
    |- summary_tables
    |   |- Sequencing_depth_ref.png
    |   | - bowtie2_alignment_rep.csv



However if the the ``-s`` parameter is set to True, the output folders will have differnt name, but contian the same types of files. 

.. code-block:: bash

    Epimapper
    |- alignment
    |   |- sam_spike_in
    |   |-  bowtie2_summary_spike_in
    |- summary_tables
    |   |- Sequencing_depth_spike_in.png
    |   | - bowtie2_alignment_spike_in.csv


If spike-in alignment is performed after alignment to reference genome the summary table and figure will contain both results.