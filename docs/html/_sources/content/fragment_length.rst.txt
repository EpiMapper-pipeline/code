=================
fragment_length
=================

.. contents::
    :local:

Evaluation of mapped fragment length distribution of input SAM files exported from high-thoughput sequencing alignment, and visulizing results.

Arguments
==========

**Required arguments:**

- ``-s, -sam``: Input file folder of SAM files exported from alignment for fragment length evaluation.

**Optional arguments:**

- ``-o, out_dir``: Output directory, default = current working directory.

Example Usage
==============

.. code-block:: bash

    $ epimapper fragment_length -s /Users/me/results/Epimapper/alignment/removeDuplicate/sam_duplicates_removed  -o /Users/me/results



Output
=======

Like all the other functions in EpiMapper Python package, the function will create a main ``Epimapper`` output directiry, if it is not already present in the chosen output directory. The function will also create one additional folder, which will contain the output. 
Here, each sample´s fragments lengths will be represented by one text file, which are combined into one CSV file (Fragment_length_all_samples.csv). This CSV file will be used to create two PNG files each containing a plot, one violon plot and one line plot.

.. code-block:: bash

    Epimapper
        |- FragmentLength
        |   |- "sample-name"_fragmentLen.txt
        |   |- Fragment_length_all_samples.csv
        |- summary_tables
        |   |- Fragment_length_violin.png
        |   |- Fragment_length_lineplot.png


