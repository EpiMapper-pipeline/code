# EpiMapper: Python Package for Data Analysis of Epigenomic Sequencing Data
## epimapper Documentation

EpiMapper is a complete pipeline for data analysis of epigenomic sequencing.

[Home](index.md) | [Quality Control of Raw Reads](fastqc.md) | [Bowtie2 Alignment](bowtie2_alignment.md) | [Removal of Duplicated Reads](remove_duplicates.md) | [Evaluation of Fragment Length Distribution](fragment_length.md) | [Quality Filtering and File Conversion](filtering.md) | [Spike-in Normalization](spike_in_calibration.md) | [SEACR Peak Calling](peak_calling.md) | [Enrichment Visualization in Heatmaps](heatmaps.md) | [Differential Analysis and Annotaion](differential_analysis.md)


## filtering

Performs data filtering for mapped reads based on their alignment quality, and file format conversion before high-level data analysis,before  visulizing reproducibility among biological replicates.

<p><strong>Required:</strong></p>
<ul>
  <li><code>-s SAM, --sam SAM </code>: Input file folder of SAM files exported from alignment for filtering </li>

  <li><code>-cs CHROMOSOME SIZES, --chromosome_sizes SAM </code>: Input file of sorted chromosome sizes information </li>


  <li><code> -bl BLACKLIST, --blacklist BLACKLIST </code>: Input file in BED format with genomic reagion that should be excluded </li>
</ul>

<p><strong>Optional, has default values:</strong></p>
<ul>
  <li><code>-o , --out_dir</code>: Output directory, default = current working directory </li>
</ul>
