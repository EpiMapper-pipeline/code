

# EpiMapper Demo run on human modification data


# 1. fastqc

epimapper fastqc -f fastq -o out

# 2. bowtie2_alignment - To reference genome Hg38

epimapper bowtie2_alignment -f fastq -i in/bowtie2_index_hg38  -m True -o out

# bowtie2_alignment - To spike-in genome (E.coli)

epimapper bowtie2_alignment -f fastq -s True -i in/bowtie2_index_ecoli -m True -o out

# 3. remove_duplicates

epimapper remove_duplicates -s out/Epimapper/alignment/sam -o out

# 4. fragment_length

#epimapper fragment_length -s out/Epimapper/alignment/removeDuplicate/sam_duplicates_removed -o out

# 5. filtering

epimapper filtering -s out/Epimapper/alignment/removeDuplicate/sam_duplicates_removed \
-cs in/hg38.chrom.sizes.clear.sorted -bl in/blacklist.bed  -sn True -o out 

# 6. spike_in_calibration

epimapper spike_in_calibration -b out/Epimapper/alignment/bed -cs in/hg38.chrom.sizes.clear.sorted \
-ss out/Epimapper/alignment/sam_spike_in -o out

# 7. peak_calling

epimapper peak_calling  -soft seacr -f out/Epimapper/alignment/bed -bg out/Epimapper/alignment/bedgraph \
-c IgG -o out

# 8. heatmap

epimapper heatmap -b out/Epimapper/alignment/bam  -p out/Epimapper/peakCalling/seacr/control \
-bl in/blacklist.bed -r in/hg38.refFlat.txt  -o out

# 9. differential_analysis

epimapper differential_analysis -p out/Epimapper/peakCalling/seacr/control \
-bg out/Epimapper/alignment/bedgraph \
-bl in/blacklist.bed -r in/hg38.refFlat.txt -cs in/hg38.chrom.sizes.clear.sorted \
-la H3K27me3_rep1 H3K27me3_rep2 -lb H3K4me3_rep1 H3K4me3_rep2 -an True \
-e  in/hg38_all_enhancers_merged_hglft_genome_327b3_4dmr.bed -o out


