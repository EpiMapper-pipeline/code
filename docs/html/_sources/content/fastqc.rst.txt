=================
fastqc 
=================

.. contents::
    :local:

Uses the `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ to asses the quality of raw FASTQ reads, and create plots.

Arguments
==========

**Required arguments:**

- ``-f, --fastq``: Path to directory containing FASTQ files with raw reads for quality control.

**Optional arguments:**

- ``-o, out_dir``: Output directory, default = current working directory


Example usage
===============

.. code-block:: bash
    
    $ epimapper fastqc -f /Users/me/my_fastq_files -o /Users/me/results


Output
=======

This function will create multiple folders containing the quality reports. Like all the other functions in EpiMapper Python package, the function will create a main ``Epimapper`` output directiry, if it is not already present in the chosen output directory.
Further, it will create a ``fastqc`` directory, which will contain several sub-folders, one for each sample which were inputed to the function. The HTML files created containts quality control summaries including: Basic Statistics, per base sequence quality, per sequence quality scores, per base sequence content, per sequence GC content, per base N content, sequence length distribution, sequence duplication levels, overrepresented sequences and adapter content.
The zip folder contains files necessary to create the HTML file, so its generally enough to just look at the HTML file conaining the summaries.

.. code-block:: bash

    Epimapper
    |- fastqc
    |   |- "sample_name"
    |   |   | - "sample-name"_R1_fastq.html
    |   |   | - "sample-name"_R2_fastq.html
    |   |   | - "sample-name"_R1_fastq.zip
    |   |   | - "sample-name"_R2_fastq.zip
