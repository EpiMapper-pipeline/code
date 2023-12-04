# EpiMapper: Python Package for Data Analysis of Epigenomic Sequencing Data
## epimapper Documentation

EpiMapper is a complete pipeline for data analysis of epigenomic sequencing.

[Home](README.md) | [Quality Control of Raw Reads](docs/fastqc.md) | [Bowtie2 Alignment](docs/bowtie2_alignment.md) | [Removal of Duplicated Reads](docs/remove_duplicates.md) | [Evaluation of Fragment Length Distribution](docs/fragment_length.md) | [Quality Filtering and File Conversion](docs/filtering.md) | [Spike-in Normalization](docs/spike_in_calibration.md) | [Peak Calling](docs/peak_calling.md) | [Enrichment Visualization in Heatmaps](docs/heatmaps.md) | [Differential Analysis and Annotaion](docs/differential_analysis.md)


# heatmap
Visualizes the enrichment of the exprimental target in predefined genomic regions and peaks by creating heatmaps.



<p><strong>Required:</strong></p>
<ul>
  <li><code>-b BAM, --bam_dir BAM </code>: Input file folder of BAM files for target enrichement heatmap. </li>

  <li><code>-p, PEAKS --peaks PEAKS </code>: Input file folder of sorted BED files from peak calling. </li>

  <li><code> -bl BLACKLIST, --blacklist BLACKLIST </code>: Input BED file with genomic reagion that should be excluded.</li>

  <li><code>-r REFERENCE, --ref REFERENCE </code>: Input reference genome RefFlat text (.txt) file.</li>



<p><strong>Optional, has default values:</strong></p>
<ul>
  <li><code>-o , --out_dir</code>: Output directory, default = current working directory. </li>
   
 <li><code>-c , --cores</code>: Defining the number of parallel processes,  default = 8. </li>

</ul>