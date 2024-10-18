#!/bin/bash

# Uses the config.txt file for configuration

source config.txt

# Below is the script for running the pipeline, you do not have to change anything here.
# However, it is a good example on how the different functions may be called!
# If you just want to use a few of the functions you may remove the the ones you dont want to use, 
# just remember that most of the functions to depent on eachother. I.e, you cannot skip straight to 
# peak_calling based on just FASTQ reads.


#Define project directory 
if [ -n "$out_dir" ]; then
    project_path="${out_dir}"

else
    
    project_path=$(pwd)
fi


# 1. Do FASTQC quality control:
epimapper fastqc -f $fastq -o "${project_path}"


# 2. Do Bowtie2 alignment:

# To reference genome:
if [ -n "$bowtie2_index" ]; then
    epimapper bowtie2_alignment -f "${fastq}" -i "${bowtie2_index}" -o "${project_path}"
    
elif [ -r "$bowtie2_index" ]; then
    epimapper bowtie2_alignment -f "${fastq}" -r "${fasta_reference}" -o "${project_path}"

fi

# To spike in genome:
if [ -n "$bowtie2_index_spike_in" ]; then
    epimapper bowtie2_alignment -f "${fastq}" -i "${bowtie2_index}" -s True -o "${project_path}"
    
elif [ -n "$fasta_reference_spike_in" ]; then
    epimapper bowtie2_alignment -f "${fastq}" -r "${fasta_reference}" -s True -o "${project_path}"
fi



# 3. Remove duplicates with picard:
sam="${project_path}/Epimapper/alignment/sam"
epimapper remove_duplicates -s "${sam}" -o "${project_path}"





# 4. Fragment length assesment:


# Check if duplication removal was preformed:
path_to_check="${project_path}/Epimapper/alignment/removeDuplicate/sam_duplicates_removed"
if [ -e "$path_to_check" ]; then
    sam="${project_path}/Epimapper/alignment/removeDuplicate/sam_duplicates_removed"
else
    sam="${project_path}/Epimapper/alignment/sam"
fi



epimapper fragment_length -s "${sam}" -o "${project_path}"




# 5. Filtering and file conversion:

epimapper filtering -s "${sam}" -bl "${blacklist}" -cs "${chromosome_sizes}"



# 6. Spike in calibration: 
path_to_check_spike_in="${project_path}/Epimapper/alignment/sam_spike_in"
if [ -e "$path_to_check_spike_in" ]; then
    echo "Preforms spike in calibration"
    bed="${project_path}/Epimapper/alignment/bed"
    epimapper spike_in_calibration -ss "${path_to_check_spike_in}" -cs "${chromosome_sizes}" -b "${bed}"
else
    echo "Do not preform spike in calibration"
fi

# 7. peak calling

# Check if percentage is defined 
if [ -n "$percentage" ]; then
    # Check if the variable is a number
    if [[ $percentage =~ ^[0-9]+$ ]]; then
        new_percentage=1
    else 
        new_percentage=0
    fi
else
    new_percentage=0
fi

if [ "$software" == "seacr" ]; then
    bedgraph="${project_path}/Epimapper/alignment/bedgraph"
    if [ "$new_percentage" -eq 1 ]; then
        if [ -n "$control_index" ]; then
            epimapper peak_calling -soft seacr -f "${bed}" -bg "${bedgraph}" -p "${percentage}" -c "${control_index}"
        else
            epimapper peak_calling -soft seacr -f "${bed}" -bg "${bedgraph}" -p "${percentage}"
        fi 
    elif [ "$new_percentage" -eq 0 ]; then

        if [ -n "$control_index" ]; then 
            epimapper peak_calling -soft seacr -f "${bed}" -bg "${bedgraph}" -c "${control_index}"
        else 
            epimapper peak_calling -soft seacr -f "${bed}" -bg "${bedgraph}" 
        fi
    fi


elif [ $"software" == "macs2" ]; then
    bam="${project_path}/Epimapper/alignment/bam"
    if [ "$new_percentage" -eq 1 ]; then
        if [ -n "$control_index" ]; then
            epimapper peak_calling -soft macs2 -f "${bed}" -b "${bam}" -gs "${genome_size}" -p "${percentage}" -c "${control_index}"
        else
            epimapper peak_calling -soft macs2 -f "${bed}" -b "${bam}" -gs "${genome_size}" -p "${percentage}"
        fi

    elif [ "$new_percentage" -eq 0 ]; then
        if [ -n "$control_index" ]; then 
            epimapper peak_calling -soft macs2 -f "${bed}" -b "${bam}" -gs "${genome_size}" -c "${control_index}"
        else 
            epimapper peak_calling -soft macs2 -f "${bed}" -b "${bam}" -gs "${genome_size}"
        fi
    fi
fi


# 8. Heatmaps
# defining bam directroy  
bam="${project_path}/Epimapper/alignment/bam"

# findinf the peaks files
if [ $"software" == "macs2" ] && [ -n "$control_index" ]; then
    peaks="${project_path}/Epimapper/peakCalling/macs2/control"

elif [ $"software" == "macs2" ] && [ -n "$control_index" ];then
    peaks="${project_path}/Epimapper/peakCalling/macs2/top_peaks"

elif [ $"software" == "seacr" ] && [ -n "$control_index" ]; then
    peaks="${project_path}/Epimapper/peakCalling/seacr/control"

elif [ $"software" == "seacr" ] && ! [ -n "$control_index" ] && [ "$new_percentage" -eq 1 ]; then
    peaks="${project_path}/Epimapper/peakCalling/seacr/top_${percentage}"
elif [ $"software" == "seacr" ] && ! [ -n "$control_index" ] && [ "$new_percentage" -eq 0] ; then
    peaks="${project_path}/Epimapper/peakCalling/seacr/top_0.01"


epimapper heatmap -b "${bam}" -bl "${blacklist}" -r "${refflat_ref}" -p "${peaks}" -c "${cores}"


# 9. Differntial analysis 

# Defining bedgraph directory
bedgraph="${project_path}/Epimapper/alignment/bedgraph"

epimapper differntial_analysis -bg "${bedgraph}" -p "${peaks}" -r "${refflat_ref}" -cs "${chromosome_sizes}" -bl "${blacklist}"\
 -la "${group_A}" -lb "${group_B}"


