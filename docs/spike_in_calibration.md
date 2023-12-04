# EpiMapper: Python Package for Data Analysis of Epigenomic Sequencing Data
## epimapper Documentation

EpiMapper is a complete pipeline for data analysis of epigenomic sequencing.


[Home](README.md) | [Quality Control of Raw Reads](docs/fastqc.md) | [Bowtie2 Alignment](docs/bowtie2_alignment.md) | [Removal of Duplicated Reads](docs/remove_duplicates.md) | [Evaluation of Fragment Length Distribution](docs/fragment_length.md) | [Quality Filtering and File Conversion](docs/filtering.md) | [Spike-in Normalization](docs/spike_in_calibration.md) | [Peak Calling](docs/peak_calling.md) | [Enrichment Visualization in Heatmaps](docs/heatmaps.md) | [Differential Analysis and Annotaion](docs/differential_analysis.md)

## spike_in_calibration
Removes experimental bias by normalizing fragment counts based on sequencing depth to a spike-in genome and visulizes results.

<p><strong>Required:</strong></p>
<ul>
  <li><code>-b BED, --bed BED </code>: Input file folder of filterd BED files for normalization. </li>

  <li><code>-ss SAM SPIKE IN --sam_spike_in</code>: Input file folder of SAM files exported from alignment to a spike in genome. </li>

  <li><code>-cs CHROMOSOME SIZES, --chromosome_sizes</code>: Input file of sorted chromosome sizes information. </li>

</ul>

<p><strong>Optional, has default values:</strong></p>
<ul>

  <li><code>-tbl ALIGNMENT SUMMARY TABLE, --fragment_table ALIGNMENT SUMMARY TABLE</code>: Input CSV file containing the following columns = ["Sample",	"Replication", "SequencingDepth", "MappedFragments", "AlignmentRate", "MappedFragments_SpikeIn",	"AlignmentRate_SpikeIn"] with corresponding sample information , default = “bowtie2_alignment_ref_and_spike_in.csv” exported by this pipeline function: <code>bowtie2_alignment</code> </li>

  <li><code>-o , --out_dir</code>: Output directory, default = current working directory. </li>

</ul>