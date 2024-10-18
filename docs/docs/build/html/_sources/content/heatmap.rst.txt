=================
heatmap
=================

.. contents::
    :local:

Visualizes the enrichment of the exprimental target in predefined genomic regions and peaks by creating heatmaps.

Arguments
==========

**Required arguments:**

- ``-b, --bam_dir``: Input file folder of BAM files for target enrichement heatmap.
- ``-p, --peaks``: Input file folder of sorted BED files from peak calling.
- ``-bl, --blacklist``: Input BED file with genomic reagion that should be excluded.
- ``-r, --ref``: Input reference genome RefFlat text (.txt) file.

**Optional arguments:**

- ``-c , --cores``: Defining the number of parallel processes,  default = 8.
- ``-o, out_dir``: Output directory, default = current working directory.


Example Usage
==============


.. code-block:: bash

    $ epimapper  -b /Users/me/results/Epimapper/alignment/bam -p /Users/me/results/Epimapper/peakCalling/seacr/control -bl /Users/me/in_folder/hg38_blacklist.bed -r  /Users/me/in_folder/hg38_refFlat.txt  -o /Users/me/results


Output
=======

Like all the other functions in EpiMapper Python package, the function will create a main ``Epimapper`` output directiry, if it is not already present in the chosen output directory. This function will output heatmap PNG files, one for each sample comparing signal to called peaks, and one combined plot comparing peaks to known genes. 
To create these figures this function will firstly create matrixes which will be stored in the input ``-p, --peaks``` directroy. 



.. code-block:: bash

    Epimapper
    |- summary_tables
    |   |- matrix_heatmap.png
    |   |- "sample-name"_heatmap.png



