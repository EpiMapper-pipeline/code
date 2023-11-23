# EpiMapper: Python Package for Data Analysis of Epigenomic Sequencing Data
## epimapper Documentation

EpiMapper is a complete pipeline for data analysis of epigenomic sequencing.

[Home](README.md) | [Quality Control of Raw Reads](docs/fastqc.md) | [Bowtie2 Alignment](docs/bowtie2_alignment.md) | [Removal of Duplicated Reads](docs/remove_duplicates.md) | [Evaluation of Fragment Length Distribution](docs/fragment_length.md) | [Quality Filtering and File Conversion](docs/filtering.md) | [Spike-in Normalization](docs/spike_in_calibration.md) | [Peak Calling](docs/peak_calling.md) | [Enrichment Visualization in Heatmaps](docs/heatmaps.md) | [Differential Analysis and Annotaion](docs/differential_analysis.md)


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
