# EpiMapper: A New Tool for Analyzing High Thotughput Sequencing from CUT&Tag
## epimapper Documentation

EpiMapper is a complete pipeline for data analysis of epigenomic sequencing.


[Home](README.md) | [Quality Control of Raw Reads](docs/fastqc.md) | [Bowtie2 Alignment](docs/bowtie2_alignment.md) | [Removal of Duplicated Reads](docs/remove_duplicates.md) | [Evaluation of Fragment Length Distribution](docs/fragment_length.md) | [Quality Filtering and File Conversion](docs/filtering.md) | [Spike-in Normalization](docs/spike_in_calibration.md) | [Peak Calling](docs/peak_calling.md) | [Enrichment Visualization in Heatmaps](docs/heatmaps.md) | [Differential Analysis and Annotaion](docs/differential_analysis.md)


# peak_calling
Finds enriched regions/calls for peaks from chromatin profiling data with [SEACR](https://github.com/FredHutch/SEACR) or [MACS2](https://github.com/macs3-project/MACS), then visulizes results.

Both of these peak calling software have an alternative where control samples are used to create a background noise level before finding enriched reagions, peaks. However, not every experiment does include control samples, therefore this a optional feature in the EpiMapper pipeline. If your experiment includes control samples and you wish it perform peak calling based on these, you may utilize the "-c, --control_index" parameter to insert the names of the control sampels. The peak_calling function will separate all the samples containing the "-c, --control_index", and use these at input for background-noise. An example would be if you have the samples: [H3K4me3_rep1, H3K4me3_rep2, H3K27me3_rep1, H3K27me3_rep2, IgG_rep1, IgG_rep1], you would input “IgG” into the "-c, --control_index". If you have the samples [healthy_rep1, healthy_rep2, cancer_rep1, cancer_rep2, control_rep1, control_rep2] you would input “control” into the "-c, --control_index".



<p><strong>Required:</strong></p>
<ul>

  <li><code>-f, FRAGMENT BED --fragments  FRAGMENT BED </code>: Input file folder of filterd BED files for peak calling. </li>

  <li><code>-soft, --software </code>:  input has to be either “macs2” or “seacr”. Here you decide which peak calling software you would like to use. Generally, MACS2 is more used and prefered for samples with higher background noise ( ie. ATAC-seq and CHiP-seq), while SEACR may be used for samples with lower noise. </li>

<p> *If using the SEACR software for peak calling:* </p>

 <li><code>-bg BEDGRAPH, --bedgraph BEDGRAPH </code>: Input file folder of BedGraph files for peak calling with SEACR. </li>

 <li><code>-s SEACR PATH, --seacr_path SEACR PATH  </code>: Input shell scurpt file of SEACR software. </li>


 <p> *If using the MACS2 software for peak calling:* </p>
<li><code> -b BAM --bam </code>: Input folder containing BAM files for MACS2 peak calling. </li>
<li><code> -gs GENOME_size --genome_size</code>: Relative genome size of organism being studied.  About 90% or 70% of the genome size. (i.e, 2.7e9 for humans, 1.87e9 for mice, 9e7 for <i>Caenorhabditis elegans</i>, or 1.2e8 for fruitfly).





</ul>

<p><strong>Optional, has default values:</strong></p>
<ul>
    <li><code>-c CONTROL_INDEX --control_index:</code> Indexes of control files (i.e, 'control' 'igG' ect.) that should be used as a background-noise reference during peak calling. If not provided the peak calling softwares will select peaks based on a cut-off value (i.e p-value for SEACR) Default = False. </li>

   <li><code> -p, --percentage </code> Cut-off percentage for peak calling software, Default= 0.01 for SEACR, and for MACS2 it is 0.05 q-value.
  <li><code>-tbl ALIGNMENT SUMMARY TABLE, --fragment_table ALIGNMENT SUMMARY TABLE</code>: Input CSV file with the following columns = ["Sample",	"Replication", "SequencingDepth", "MappedFragments", "AlignmentRate"] with corresponding sample information, default = “bowtie2_alignment_ref.csv” exported by this pipeline function: <code> bowtie2_alignment </code>. </li>
  <li><code>-o , --out_dir</code>: Output directory, default = current working directory. </li>
</ul>