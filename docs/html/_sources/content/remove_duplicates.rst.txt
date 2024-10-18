=================
remove_duplicates
=================


.. contents::
    :local:

Uses `Picard <https://broadinstitute.github.io/picard/>`_  to remove duplicated reads mapped to the same place in of a reference genome during alignment, and visulizing results.

Arguments
==========

**Required arguments:**

- ``-s , --sam``: Input file folder of SAM files exported from alignment to remove duplicates reads.

**Optional arguments:**

- ``-o, out_dir``: Output directory, default = current working directory.



Example usage:
===============

.. code-block:: bash

    $ epimapper remove_duplicates -f /Users/me/results/Epimapper/alignment/sam  -o /Users/me/results


Output
=======

Like all the other functions in EpiMapper Python package, the function will create a main ``Epimapper`` output directiry, if it is not already present in the chosen output directory. The function will also create two additional folders, which will contain the output. 
Here, the "sam_duplicated_removed" directory will contain the duplicate removed sam files, and the "picard_summary" folder will contian text files for each sample with dupliaction removal statistics. These text files are used to create the "Duplication_summary.csv" located in the "summary_tables" folder. Lastly, a PNG image file containing boxplots of duplication rate, estimated libary size and number of unique fragments is created and also stored in the "summary_tabels" folder.


.. code-block:: bash

    Epimapper
        |- alignment
        |   |- removeDuplicate
        |   |   |- sam_duplicated_removed
        |   |   |- picard_summary
        |- summary_tables
        |   |- Duplication_summary.csv
        |   |- Duplication_rate.png 



