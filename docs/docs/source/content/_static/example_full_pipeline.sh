#!/bin/bash

# Uses the config.txt file for configuration
source config.txt

# Below is the script for running the pipeline, you do not have to change anything here.
# However, it is a good example on how the different functions may be called!
# If you just want to use a few of the functions you may remove the the ones you dont want to use, 
# just remember that most of the functions to depent on eachother. I.e, you cannot skip straight to 
# peak_calling based on just FASTQ reads.


# Function to check if a path exists
path_exists() {
    [[ -e "$1" ]]
}


# Function to perform FASTQC quality control
perform_fastqc() {
    epimapper fastqc -f "$fastq" -o "$out_dir"
}

# Function to perform Bowtie2 alignment
perform_bowtie2_alignment() {
    if [ -n "$bowtie2_index" ]; then
        epimapper bowtie2_alignment -f "$fastq" -i "${bowtie2_index}" -o "$out_dir"
    elif [ -n "$fasta_reference" ]; then
        epimapper bowtie2_alignment -f "$fastq" -r "${fasta_reference}" -o "$out_dir"
    fi

   if [ -n "$bowtie2_index_spike_in" ]; then
        epimapper bowtie2_alignment -f "$fastq" -i "${bowtie2_index_spike_in}" -o "$out_dir" -s True
    elif [ -n "$fasta_reference_spike_in" ]; then
        epimapper bowtie2_alignment -f "$fastq" -r "${fasta_reference_spike_in}" -o "$out_dir" -s True
    fi
}



# Function to remove duplicates with Picard
remove_duplicates_with_picard() {
    sam="${out_dir}/Epimapper/alignment/sam"
    epimapper remove_duplicates -s "$sam" -o "$out_dir"
}


# Function to assess fragment length
assess_fragment_length() {
    if [ "$rm_duplicate" = true ]; then
        sam="${out_dir}/Epimapper/alignment/removeDuplicate/sam_duplicates_removed"
    else
        sam="${out_dir}/Epimapper/alignment/sam"
    fi 
    
    epimapper fragment_length -s "$sam" -o "$out_dir"
}


# Function for filtering and file conversion
filtering_and_conversion() {
    if [ "$rm_duplicate" = true ]; then
        sam="${out_dir}/Epimapper/alignment/removeDuplicate/sam_duplicates_removed"
    else
        sam="${out_dir}/Epimapper/alignment/sam"
    fi

    epimapper filtering -s "$sam" -bl "$blacklist" -cs "$chromosome_sizes" -o "$out_dir"
}

# Function for spike-in calibration
spike_in_calibration() {
    path_to_check_spike_in="${out_dir}/Epimapper/alignment/sam_spike_in"
    if path_exists "$path_to_check_spike_in"; then
        echo "Performing spike-in calibration"
        bed="${out_dir}/Epimapper/alignment/bed"
        epimapper spike_in_calibration -ss "$path_to_check_spike_in" -cs "$chromosome_sizes" -b "$bed" -o "$out_dir"
    else
        echo "Do not perform spike-in calibration"
    fi
}

# Function for peak calling
perform_peak_calling() {
    bed="${out_dir}/Epimapper/alignment/bed"
    if [ "$software" == "seacr" ]; then
        bedgraph="${out_dir}/Epimapper/alignment/bedgraph"
        if [ -n "$control_index" ]; then
            epimapper peak_calling -soft seacr -f "$bed" -bg "$bedgraph" ${percentage:+-p "$percentage"} -c "$control_index" -o "$out_dir"
        else
            epimapper peak_calling -soft seacr -f "$bed" -bg "$bedgraph" ${percentage:+-p "$percentage"} -o "$out_dir"
        fi
    elif [ "$software" == "macs2" ]; then
        bam="${out_dir}/Epimapper/alignment/bam"
        if [ -n "$control_index" ]; then
            epimapper peak_calling -soft macs2 -f "$bed" -b "$bam" -gs "$genome_size" ${percentage:+-p "$percentage"} -c "$control_index" -o "$out_dir"
        else
            epimapper peak_calling -soft macs2 -f "$bed" -b "$bam" -gs "$genome_size" ${percentage:+-p "$percentage"} -o "$out_dir"
        fi
    fi
}


# Function for generating heatmaps
generate_heatmaps() {
    bam="${out_dir}/Epimapper/alignment/bam"
    case $software in
        macs2)
            peaks="${out_dir}/Epimapper/peakCalling/macs2/$(if [ -n "$control_index" ]; then echo "control"; else echo "top_peaks"; fi)"
            ;;
        seacr)
            peaks="${out_dir}/Epimapper/peakCalling/seacr/$(if [ -n "$control_index" ]; then echo "control"; else echo "top_${percentage:-0.01}"; fi)"
            ;;
    esac
    epimapper heatmap -b "$bam" -bl "$blacklist" -r "$refflat_ref" -p "$peaks" -c "$cores" -o "$out_dir"
}

# Function for differential analysis
perform_differential_analysis() {
    bedgraph="${out_dir}/Epimapper/alignment/bedgraph"
    epimapper differential_analysis -bg "$bedgraph" -p "$peaks" -r "$refflat_ref" -cs "$chromosome_sizes" -bl "$blacklist" -la "$group_A" -lb "$group_B" -o "$out_dir"
}

# Main execution
perform_fastqc
perform_bowtie2_alignment

if [ "$rm_duplicate" = true ]; then
   # remove_duplicates_with_picard
fi

assess_fragment_length
filtering_and_conversion
spike_in_calibration
perform_peak_calling
generate_heatmaps
perform_differential_analysis
