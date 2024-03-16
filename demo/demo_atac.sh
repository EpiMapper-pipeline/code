# EpiMapper demo run on human ATAC-seq data (only chr21)



# 1. fastqc

epimapper fastqc -f fastq -o out

# 2. bowtie2_alignment 

epimapper bowtie2_alignment -f fastq -i in/hg19_chr21_bowtie2_index -o out

# 3. remove_duplicates

epimapper remove_duplicates -s out/Epimapper/alignment/sam -o out

# 4. fragment_length

epimapper fragment_length -s out/Epimapper/alignment/removeDuplicate/sam_duplicates_removed -o out

# 5. filtering

epimapper filtering -s /Users/eier/Documents/demo/ATAC/out/Epimapper/alignment/removeDuplicate/sam_duplicates_removed \
-cs in/hg19_chromosome_sizes_sorted.txt -bl in/hg19-blacklist_sorted.bed -atac True -o /Users/eier/Documents/demo/ATAC/out

# 6. peak_calling 

epimapper peak_calling -soft macs2 -f /Users/eier/Documents/demo/ATAC/out/Epimapper/alignment/bed -b /Users/eier/Documents/demo/ATAC/out/Epimapper/alignment/bam \
-gs 2.7e9  -o /Users/eier/Documents/demo/ATAC/out

# 7. heatmaps

epimapper heatmap -b out/Epimapper/alignment/bam -bl in/hg19-blacklist_sorted.bed \
-p out/Epimapper/peakCalling/macs2/top_peaks -r in/hg19.refFlat_chr21.txt -o /Users/eier/Documents/demo/ATAC/out


#8. differntial_analysis 

epimapper differential_analysis -p out/Epimapper/peakCalling/macs2/top_peaks \
-r in/hg19.refFlat_chr21.txt  -bl in/hg19-blacklist_sorted.bed -cs in/hg19_chromosome_sizes_sorted_filtered.txt \
-fold True -an True -e in/hg19_all_enhancers_merged_4dmr.bed -o out \
-la diabetic-1_rep1 diabetic-2_rep1 diabetic-3_rep1 diabetic-4_rep1 diabetic-5_rep1 diabetic-6_rep1 \
-lb healthy-1_rep1 healthy-2_rep1 healthy-3_rep1 healthy-4_rep1 healthy-5_rep1 healthy-6_rep1 healthy-7_rep1 healthy-8_rep1 healthy-9_rep1
