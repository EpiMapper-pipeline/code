# EpiMapper: Python Package for Data Analysis of Epigenomic Sequencing Data
## epimapper Documentation

EpiMapper is a complete pipeline for data analysis of epigenomic sequencing.



[Home](index.md) | [Quality Control of Raw Reads](docs/fastqc.md) | [Bowtie2 Alignment](docs/bowtie2_alignment.md) | [Removal of Duplicated Reads](docs/remove_duplicates.md) | [Evaluation of Fragment Length Distribution](docs/fragment_length.md) | [Quality Filtering and File Conversion](docs/filtering.md) | [Spike-in Normalization](docs/spike_in_calibration.md) | [Peak Calling](docs/peak_calling.md) | [Enrichment Visualization in Heatmaps](docs/heatmaps.md) | [Differential Analysis and Annotaion](docs/differential_analysis.md)


## fastqc

Performs quality control on raw reads fastq files from high-thoughput sequencing.
<p><strong>Required:</strong></p>
<ul>
  <li><code>-f FASTQ, --fastq FASTQ </code>: Input file folder of fastq files with raw reads for quality control.</li>
  
</ul>

<p><strong>Optional, has default values:</strong></p>
<ul>
  <li><code>-o , --out_dir</code>: Output directory, default = current working directory </li>
</ul>
