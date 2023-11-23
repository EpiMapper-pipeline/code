# EpiMapper

EpiMapper: Python Package for Data Analysis of Epigenomic Sequencing Data


usage:  epimapper task args

## Tasks available for using:
<ul>
  <li>
    <code>fastqc</code>
    - Performs quality control on raw reads fastq files from high-thoughput sequencing.
  </li>
  <li>
    <code>bowtie2_alignment</code>
    - Mapping reads to a reference genome with Bowtie2 alignment of fastq sequencing reads files from high-thoughput sequecing, and visulizing results.
  </li>
  <li>
    <code>remove_duplicates</code>
    - Remove duplicated reads mapped to the same place in of a reference genome during alignment, and visulizing results.
  </li>
  <li>
    <code>fragment_length</code>
    - Evaluation of mapped fragment length distribution of input SAM files exported from high-thoughput sequencing alignment, and visulizing results.
  </li>
  <li>
    <code>filtering</code>
    - Performs data filtering for mapped reads based on their alignment quality, and file format conversion before high-level data analysis,before  visulizing reproducibility among biological replicates.
  </li>
  <li>
    <code>spike_in_calibration</code>
    - Removes experimental bias by normalizing fragment counts based on sequencing depth to a spike-in genome and visulizes results.
  </li>
  <li>
    <code>peak_calling</code>
    - Finds enriched regions/calls for peaks from chromatin profiling data with SEACR, then visulizes results.
  </li>
  <li>
    <code>heatmap</code>
    - Visualizes the enrichment of target protein in predefined genomic regions and peaks by creating heatmaps.
  </li>
  <li>
    <code>differntial_analysis</code>
    - Preforms differntial analysis on enriched reagions/peaks before annotating the stastitically significant changes to spesific genomic reagions and visulizing the results. 
  </li>
</ul>



positional arguments:
  task        Pipeline task to run

optional arguments:
  -h, --help  show this help message and exit


