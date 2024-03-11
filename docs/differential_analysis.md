# EpiMapper: A New Tool for Analyzing High Thotughput Sequencing from CUT&Tag
## epimapper Documentation

EpiMapper is a complete pipeline for data analysis of epigenomic sequencing.

[Home](README.md) | [Quality Control of Raw Reads](docs/fastqc.md) | [Bowtie2 Alignment](docs/bowtie2_alignment.md) | [Removal of Duplicated Reads](docs/remove_duplicates.md) | [Evaluation of Fragment Length Distribution](docs/fragment_length.md) | [Quality Filtering and File Conversion](docs/filtering.md) | [Spike-in Normalization](docs/spike_in_calibration.md) | [Peak Calling](docs/peak_calling.md) | [Enrichment Visualization in Heatmaps](docs/heatmaps.md) | [Differential Analysis and Annotaion](docs/differential_analysis.md)


# differntial_analysis
Preforms differntial analysis on enriched reagions/peaks before annotating the stastitically significant changes to spesific genomic regions and visulizing the results. 

For ATAC-seq, fold enrichment is utilized since ATAC-seq experiments usually do not include control samples. Thus, the read counts from alignment are affected by sequencing depth, hence, the samples should not be compared directly.

The annotation step is derived from another epigenetic Python Package: [HMST-seq-Analyzer](https://hmst-seq.github.io/hmst/) and allows for annotation to gene regions, transcription start sites (TSS), transcription end sites (TES), and end 5´ distance regions. Further, annotation to enhancer regions may be done, however, a separate enchancer BED must provided in the -e, --enhancer parameter. These files must be processed to have this format:



<p><strong>Required:</strong></p>
<ul>
  <li><code>-p, PEAKS --peaks PEAKS </code>: Input file folder of  sorted  BED files from peak calling. </li>

  <li><code>-cs CHROMOSOME SIZES, --chromosome_sizes </code>: Input file of sorted chromosome sizes information. </li>

  <li><code> -bl BLACKLIST, --blacklist BLACKLIST </code>: Input file with genomic reagion that should be excluded.</li>

  <li><code> -r REFERENCE --refrence_refFlat</code>: RefFlat text (.txt) file to which the peaks will be annotated to. The reference file will be used for the division of genomic regions (i.e., Transcription start site/promotor, intrageneic etc.)". </li>

  <li><code>-la LIST_A, --list_a LIST_A </code>: List of sample names that will be used to do differential analysis. </li>

   <li><code>-lb LIST_B, --list_b LIST_B </code>: List of sample names that will be used to do differential analysis. </li>



<p><strong>Optional, has default values:</strong></p>
<ul>
  <li><code>-o , --out_dir</code>: Output directory, default = current working directory </li>

  <li><code>-fold , --fold_enrichment </code>: Input must be either "True" or "False". The function will use fold enrichemnt as well as normalization before differential analysis. Default=False</li>

<li><code>-an , --annotate</code>: Input must be either "True" or "False" Weither to annotate the differntial peaks or not. This may be beneficial for some analysis protocols, however it may be a slightly time consuming taske. Default=False  </li>

  <li><code>-bg BEDGRAPH, --bedgraph BEDGRAPH </code>: Input file folder of BedGraph fragment files for differntial analysis. No default, this input is required if the data is not ATAC-seq data. </li>
  <li><code> -cut, --p_value_cutoff </code>: Cut-off p-value for differential analysis, default =0.05 </li>

  <li><code>-n , --normalize </code>: Whether or not to quantile normalize the data. Might be beneficial if spike in calibration is not preformed. Default = False. </li>

  <li><code>-X </code>: The number of upstream basepaires from TSS, TES, gene, when creating genomic region files. Default = 1000. </li>
  <li><code>-Y </code>: The number of downstream basepairs from TSS, TES, gene, when creating genomic region files. Default = 1000. </li>
  <li><code>-M </code>: The number of basepairs from gene start site, 5dist, when creating genomic region files. Default=10000.
  <li><code>-N </code>: The number of basepairs from gene start site, 5dits, when creating genomic region files. Default=1000000.</li>
  <li><code>-l , --minIntergenicLen</code>: Minimim intergenic region distance. Default = 2000. </li>
  <li><code>-xL, --maxIntergenicLen</code>: Maximum intergentic region distance. Default = 10000000. </li>
  <li><code>-i, --intergenicBTGenes </code>: Wheither intergenic regions is considrered between gene body regions (True), or betweeen TSS and TES (False). Default=True. </li>
  <li><code>-e, --enhancer </code>: Sorted BED file with defined enhancer regions for annotation. </li>

</ul>


