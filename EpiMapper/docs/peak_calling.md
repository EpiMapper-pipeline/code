# EpiMapper: Python Package for Data Analysis of Epigenomic Sequencing Data
## epimapper Documentation

EpiMapper is a complete pipeline for data analysis of epigenomic sequencing.


[Home](index.md) | [Quality Control of Raw Reads](fastqc.md) | [Bowtie2 Alignment](bowtie2_alignment.md) | [Removal of Duplicated Reads](remove_duplicates.md) | [Evaluation of Fragment Length Distribution](fragment_length.md) | [Quality Filtering and File Conversion](filtering.md) | [Spike-in Normalization](spike_in_calibration.md) | [SEACR Peak Calling](peak_calling.md) | [Enrichment Visualization in Heatmaps](heatmaps.md) | [Differential Analysis and Annotaion](differential_analysis.md)


# peak_calling
Finds enriched regions/calls for peaks from chromatin profiling data with SEACR, then visulizes results

<p><strong>Required:</strong></p>
<ul>

  <li><code>-f, FRAGMENT BED --fragments  FRAGMENT BED </code>: Input file folder of filterd BED files for peak calling </li>

<p> *If using the SEACR software for peak calling:* </p>

 <li><code>-bg BEDGRAPH, --bedgraph_ex BEDGRAPH </code>: Input file folder of BedGraph files for peak calling with SEACR </li>

 <li><code>-s SEACR PATH, --seacr_path SEACR PATH  </code>: Input shell scurpt file of SEACR software </li>


 <p> *If using the MACS2 software for peak calling:* </p>
<li><code> -b BAM --bam </code>: Input folder containing BAM files for MACS2 peak calling </li>
<li><code> -gs GENOME_size --genome_size</code>: Relative genome size of organism being studied.  About 90% or 70% of the genome size ()





</ul>

<p><strong>Optional, has default values:</strong></p>
<ul>
    <li><code>-c CONTROL_INDEX --control_index: Indexes of control files (i.e, 'control' 'igG' ect.) that should be used as a background-noise reference during peak calling. If not provided the peak calling softwares will select peaks based on a cut-off value (i.e p-value for SEACR) Default = False </li>

   <li><code> -p, --percentage </code> Cut-off percentage for peak calling software, Default= 0.01 for SEACR, and for MACS2 it is 0.05 q-value
  <li><code>-tbl ALIGNMENT SUMMARY TABLE, --fragment_table ALIGNMENT SUMMARY TABLE: Input CSV file with the following columns = ["Sample",	"Replication", "SequencingDepth", "MappedFragments", "AlignmentRate"] with corresponding sample information, default = “bowtie2_alignment_ref.csv” exported by this pipeline function: <code> bowtie2_alignment </code> </li>
  <li><code>-o , --out_dir</code>: Output directory, default = current working directory </li>
</ul>