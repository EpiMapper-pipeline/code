# CUT&Tag-Analyzer: A Complete Pipeline for CUT&Tag Data Analysis
## cut_and_tag_analyzer Documentation

CUT&Tag-Analyzer is a complete pipeline for data analysis of CUT&Tag sequencing data.

[Home](index.md) | [Quality Control of Raw Reads](fastqc.md) | [Bowtie2 Alignment](bowtie2_alignment.md) | [Removal of Duplicated Reads](remove_duplicates.md) | [Evaluation of Fragment Length Distribution](fragment_length.md) | [Quality Filtering and File Conversion](filtering.md) | [Spike-in Normalization](spike_in_calibration.md) | [SEACR Peak Calling](peak_calling.md) | [Enrichment Visualization in Heatmaps](heatmaps.md) | [Differential Analysis and Annotaion](differential_analysis.md)

# differntial_analysis
Preforms differntial analysis on enriched reagions/peaks before annotating the stastitically significant changes to spesific genomic reagions and visulizing the results. 



<p><strong>Required:</strong></p>
<ul>
  <li><code>-bg BEDGRAPH, --bedgraph BEDGRAPH </code>: Input file folder of BedGraph fragment files for differntial analysis </li>

  <li><code>-p, PEAKS --peaks PEAKS </code>: Input file folder of  BED files from SEACR peak calling </li>

  <li><code>-cs CHROMOSOME SIZES, --chromosome_sizes SAM </code>: Input file of sorted chromosome sizes information </li>

  <li><code> -bl BLACKLIST, --blacklist BLACKLIST </code>: Input file with genomic reagion that should be excluded</li>

  <li><code> -r REFERENCE --refrence_refFlat</code>: RefFlat text (.txt) file to which the peaks will be annotated to. The reference file will be used for the division of genomic regions (i.e., Transcription start site/promotor, intrageneic etc." </li>

  <li><code>-la LIST_A, --list_a LIST_A </code>: List of sample names that will be used to do differential analysis </li>

   <li><code>-lb LIST_B, --list_b LIST_B </code>: List of sample names that will be used to do differential analysis </li>



<p><strong>Optional, has default values:</strong></p>
<ul>
  <li><code>-o , --out_dir</code>: Output directory, default = current working directory </li>

  <li><code>-n , --normalize </code>: Whether or not to quantile normalize the data. Might be beneficial if spike in calibration is not preformed. Default = False </li>

  <li><code>-X </code>: The number of upstream basepaires from TSS, TES, gene, when creating genomic region files. Default = 1000 </li>
  <li><code>-Y </code>: The number of downstream basepairs from TSS, TES, gene, when creating genomic region files. Default = 1000 </li>
  <li><code>-M </code>: The number of basepairs from gene start site, 5dist, when creating genomic region files. Default=10000.
  <li><code>-N </code>: The number of basepairs from gene start site, 5dits, when creating genomic region files. Default=1000000</li>
  <li><code>-l , --minIntergenicLen</code>: Minimim intergenic region distance. Default = 2000 </li>
  <li><code>-xL, --maxIntergenicLen</code>: Maximum intergentic region distance. Default = 10000000 </li>
  <li><code>-i, --intergenicBTGenes </code>: Wheither intergenic regions is considrered between gene body regions (True), or betweeen TSS and TES (False). Default=True </li>
  <li><code>-e, --enhancer </code>: BED file with defined enhancer regions for annotation. </li>

</ul>


