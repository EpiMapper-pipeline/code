This contains the following functions to use for Cut&tag analysis:

* For every function activated the folder CutAndTagAnalyzer will be made as an outfielder for all steps.
* Appropriate folders will be made inside of the CutAndTagAnalyzer directory, depending on the function.
* Input reads must have the following name structure: [sample name/histone]_rep[replication number]_[R1/R2)]_[(technical replicates]].fastq
	Ex: 
	H3K4me3_rep1_R1.fastq
	H3K4me3_rep1_R1_001.fastq

*All plots will be found in the summary_tabels folder


1. fastqc.py 
	- Outputs in the «fastqc» folder located inside CutAndTagAnalyzer

	- Required input:

		* --fastq (-f): The full pathway to fastq, filenames must contain either R1 or R2 depending on the read number.
	
	- Optional input:
		* --out_dir (-o): The full pathway to the desired output directory, default = current directory
		
	- Output: folders with fastqc statistics for the input fastq files
	

2. bowtie2_alignment.py


	- Outputs in the «sam», «bowtie2_summary» or «sam_spike_in» and «bowtie2_summary_spike_in» as well as the «summary_tables» folder
	 located inside CutAndTagAnalyzer, depending on if the alignment is spike in alignment or not. 


	 -Required input:


        
		* --fastq_dir (-f): The full pathway to fastq, filenames must contain either R1 or R2 depending on the read number.
        
		* --bowtie2_index_pathway(-i): The full pathway to a bowtie2 index directory, either for spike in or normal alignment
        
        
	-Optional input:

       
        	* --out_dir(-o): The full pathway to the desired output directory, where the CutAndTagAnalyzer directory will be made
        
		* --spike_in(-s): bool, default = False, if the alignment is spike in or not 
        
       		* --merge_technical_replicates(-m), default = False, if there need for merging technical replicates


    
	- Output: SAM files, bowtie2 summary files, CVS summary table of the alignment, and plots of various alignment statistics


3. remove_duplicates.py


	     
    Required input:
        
        * --sam (-s), str, the full pathway to the directory with sam files



    Optional input: 
        
        * --out_dir (-o), the full pathway to a chosen output directory, default is the current directory



      
    - Outputs in the «removeDuplicate», «picard_summary» directory as well as in «summary_tables»


    - Output: SAM files, Picard summary txt files, CSV summary table, and plots of duplication rate


4.fragment_length.py


	  - Required input:

            * --sam(-s), str, the full path to the sam directory containing the sam files being analyzed
            

        - Optional input:

            * --out_dir(-o), str, The full pathway to the desired output directory, where the CutAndTagAnalyzer directory will be made


	-Outputs in "summary_tabels"

	-Outputs: summary table of the fragment lengths, 2 summary plots

5. filtering.py


	 - Required input:


            * --sam (-s), str, the full path to the sam directory containing the sam files being analyzed
            
            * --chromosome_sizes (-cs), str, the full path to the file containing chromosome size information about the genome
            
            * --blacklist (-bl), str, the full path to bed file with blacklisted regions of the genome
            
            
            
        - Optional input:

            * --out_dir (-o), str, The full pathway to the desired output directory, where the 	CutAndTagAnalyzer directory will be made
	
	-Output in the "summary_tabels" directory


	- Output: Two correlation plots, one log2-correlation and one log2-correlation matrix filtered with the removal of bins with lower than 8 counts.

6. spike_in_calibration.py


    - Implements bedtools genomecov to calibrate bed files based on alignment to a spike in genome

    - Creates normalized bedgraph files.

    - Creates some basic plots based on the calibration

    - Required input:
            
            * --sam_spike_in(-ss), str, full path to sam directory containing the sam files bowtie2 alignment to spike in genome
            
            * --chromosome_sizes (-cs), str, full path to file containing chromosome size information about the genome
            
            * --bed (-b), str, full path to directory with bed files being analyzes
            
            
    - Optional input:
            
            * --out_dir(-o), str, The full pathway to desired output directory, where CutAndTagAnalyzer directory will be made
            
            * --fragment_table, (-tbl), str, full path to table with information about number of mapped fragments, column names = ["Sample",	"Replication", "SequencingDepth", "MappedFragments", "AlignmentRate", "MappedFragments_SpikeIn",	"AlignmentRate_SpikeIn"]
                             Default: will be collected from summary_tables, made in python script "bowtie2_alignment.py"



7. peak_calling.py
 
    - Implements SEACR to find peaks in bed files.

     -Required input:
        
        * --bedgraph_ex (-bg), full pathway to bedgraph directory with files being analyzed
        
        * --bedgraph_control (-c), full pathway to bedgraph file that will be used as control 
        
        * --seacr_path (-s), full pathway to seacr shell script
        
        * --fragments (-f), full pathway to directory with fragment bed files being used to find FRiPs (Fragment in peaks)
        
        
	-Optional input:
        
        * --out_dir(-o): The full pathway to desired output directory, where CutAndTagAnalyzer directory will be made

         *--fragment_table, (-tbl), str, full path to table with information about number of mapped fragments, column names = ["Sample",	"Replication", "SequencingDepth", "MappedFragments", "AlignmentRate", "MappedFragments_SpikeIn",	"AlignmentRate_SpikeIn"]
                          Default: will be collected from summary_tables, made in python script "bowtie2_alignment.py"
    

8. heatmap.py 

    - Implements deeptools to create heatmaps of transcription units and cut&tag peaks 

     -Required input:
        
        * --bam_dir (-b), full pathway to directory with bam files being used to make bigwig files used to create matrixes
        
        * --ref (-r), full pathway to sorted reference bed file of genome 
        
        * --peaks (-p), full pathway to directory with bed files containing peaks from SEACR peak calling
        
        * --blacklist (-bl), full pathway to bed file containing blacklisted parts of the genome
        
        
	-Optional input:
        
        * --out_dir(-o): The full pathway to desired output directory, where CutAndTagAnalyzer directory will be made

        * --cores (-c), number of cores being used in matrix calculation, default = 8

  7. differential_analysis.py - Main Python script for performing the DAR analysis

   -Required input:
        
        
        * --peaks (-p): The full pathway to the directory containing files from seacr peak calling
        
        * --bedgraph (-b): The full pathway to the directory containing the bedgraph files
        
        * --chromosome_sizes (-cs): The full pathway to the file containing chromosome size information
        
        * --blacklist (-bl): The full pathway to the bed file containing blacklisted regions of the genome
        
        * --list_a (-la): List of samples that will be compared to list b
        
        * --list_b (-lb): List of samples that will be compared to list a
        
         
     -Optional input:
        
        * --out_dir(-o): The full pathway to desired output directory, where CutAndTagAnalyzer directory will be made


     -Output:

          *Creates new directories, if not present, in the output directory:
        
        	 - CutAndTagAnalyzer
        	 - CutAndTagAnalyzer/summary_tables
       		 - CutAndTagAnalyzer/differential_analysis
       		 - CutAndTagAnalyzer/differential_analysis/out_combined_files
        	 -  CutAndTagAnalyzer/differential_analysis/DAR
        
	* hg38.100b.windows.bed - blacklist window file
	* hg38.100.b.windows.BlackListFiltered.bed - blacklist filtered hg38 reference genome file
	* XXX_repX__hg38_XY.100b.windows.BlackListFiltered.bed.gz - blacklist filtered bed files from samples
	* combined_peaks_seacr_top01.bed - combined seacr peaks from all samples
	* combined_peaks_seacr_top01_merged.bed - combined and merged seacer peaks from all samples
	* combined_signals_100b.bed.gz - combined signals from all samples 
	* combined_signals_100b.head - combined signals from all samples header
	* combined_signals_100b_quantLog.bed.gz - quantile and log transformed combined signals from all samples
	* mapped_peaks_100bp_top01_merged.bed.gz - peaks mapped to signals
	* /genome 
		- containing the mapped peaks overlap with reference genome files found in /data
		- files describing from which part of the genome the peaks are found
		- CutAndTag_DAR_ttest_pval_0.001.csv - the t test results
		- CutAndTag_DAR_ttest_pval_0.001.pdf - plot of t test results
	*combined_peaks_seacr_top01_merged_pvals0.05.csv - file with annotated peaks 

	dar_files - Folder with necessary files for performing DAR analysis 
	- /genomes/hg38 
		* hg38.chrom.sizes.clear.sorted - chromosome sizes
		* hg38.refFlat.txt - hg38 genome reference
		* blacklist.bed - blacklisted parts of hg38


	- /data 
	 	* 5dist_Down1000000_Down5000removedShort.bed - reference bedfile 5dist down
		* 5dist_Up1000000_Up5000removedShort.bed - reference bedfile 5dist up 
		* gene_Up5000_Down1000removedShort.bed - reference bedfile for genes
		* hg38_all_enhancers_merged_hglft_genome_327b3_4dmr.bed - reference bedfile for enhancers 
		*hg38.refFlat_clean_sorted.bed - reference hg38 genome bedfile
		*hg38.refFlat_clean_sorted.txt - reference hg38 gene/RNA name reference file
		*intergenic_uniqueSorted_betweenTSS_TES_genes_minLen500.bed - reference bedfile for intergenics
		*removed_regions_all_TSS_TES_5dist_geneBodyLessThan0.bed - reference bedfile 
		* TES_Up5000_Down1000removedShort.bed - reference bedfile
		*TSS_Up5000_Down1000_removedShort.bed -reference bedfile
		* list_region_files.txt - text file with paths to region files     
