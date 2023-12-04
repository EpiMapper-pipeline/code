# EpiMapper: Python Package for Data Analysis of Epigenomic Sequencing Data
## epimapper Documentation

EpiMapper is a complete pipeline for data analysis of epigenomic sequencing.


[Home](index.md) | [Quality Control of Raw Reads](fastqc.md) | [Bowtie2 Alignment](bowtie2_alignment.md) | [Removal of Duplicated Reads](remove_duplicates.md) | [Evaluation of Fragment Length Distribution](fragment_length.md) | [Quality Filtering and File Conversion](filtering.md) | [Spike-in Normalization](spike_in_calibration.md) | [SEACR Peak Calling](peak_calling.md) | [Enrichment Visualization in Heatmaps](heatmaps.md) | [Differential Analysis and Annotaion](differential_analysis.md)



## bowtie2_alignment

Mapping reads to a reference genome with Bowtie2 alignment of fastq sequencing reads files from high-thoughput sequecing, and visulizing results.

<p><strong>Required:</strong></p>
<ul>
  <li><code>-f FASTQ, --fastq FASTQ </code>: Input file folder of fastq files for aligment.</li>

  
</ul>

<p><strong>Optional, has default values:</strong></p>
<ul>
   <li><code>-i --bowtie2_index_pathway </code>: Input file folder of Bowtie2 reference genome indexing files.</li>
   <li><code>-r, --reference: </code>: Input reference genome FASTA file for the creation of Bowtie2 reference genome indexing files if the user does not have indexing files avalible.
  <li><code> -s , --spike_in</code>: If the alignment is spike-in, default = False </li>
  <li><code> -m , --merge</code>: Merges technical replicates, default = False </li>
  <li><code> -o , --out_dir</code>: Output directory, default = current working directory </li>
</ul>
