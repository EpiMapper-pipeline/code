# CUT&Tag-Analyzer: A Complete Pipeline for CUT&Tag Data Analysis
## cut_and_tag_analyzer Documentation

CUT&Tag-Analyzer is a complete pipeline for data analysis of CUT&Tag sequencing data.

[Home](index.md) | [Quality Control of Raw Reads](fastqc.md) | [Bowtie2 Alignment](bowtie2_alignment.md) | [Removal of Duplicated Reads](remove_duplicates.md) | [Evaluation of Fragment Length Distribution](fragment_length.md) | [Quality Filtering and File Conversion](filtering.md) | [Spike-in Normalization](spike_in_calibration.md) | [SEACR Peak Calling](peak_calling.md) | [Enrichment Visualization in Heatmaps](heatmaps.md) | [Differential Analysis and Annotaion](differential_analysis.md)


# heatmap
Visualizes the enrichment of target protein in predefined genomic regions and peaks by creating heatmaps



<p><strong>Required:</strong></p>
<ul>
  <li><code>-b BAM, --bam_dir BAM </code>: Input file folder of BAM files for target enrichement heatmap </li>

  <li><code>-p, PEAKS --peaks PEAKS </code>: Input file folder of  BED files from SEACR peak calling </li>

  <li><code>-c CONTROL BEDGRAPH --bedgraph_control: Input control BedGraph file used to define background siganl in peak calling </li>

  <li><code> -bl BLACKLIST, --blacklist BLACKLIST </code>: Input file with genomic reagion that should be excluded</li>

  <li><code>-r REFERENCE, --ref REFERENCE </code>: Input reference genome BED file </li>



<p><strong>Optional, has default values:</strong></p>
<ul>
  <li><code>-o , --out_dir</code>: Output directory, default = current working directory </li>
   
 <li><code>-c , --cores</code>: Defining the number of parallel processes,  default = 8 </li>

</ul>