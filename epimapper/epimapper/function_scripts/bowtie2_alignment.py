#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import os
import glob
import pathlib as pl
import re 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mlp
import subprocess


# Add to every plot function
mlp.use("agg")


def set_parser(parser):
    
    """
    
    Creates and returns variables from terminal inputs
    
    Function input from shell or terminal:
    
    -Required input:
        
		* --fastq (-f): The full pathway to fastq, filenames must contain either R1 or R2 depending on the read number.
        
        
        
	-Optional input:
        
        * --out_dir(-o): The full pathway to desired output directory, where Epimapper directory will be made
        
		* --spike_in(-s): bool, default = False, if the alignment is spike in or not 
        
        * --merge_technical_replicates(-m), default = False, if there need for merging technical replicates
    
        * --bowtie2_index_pathway(-i):, no default, The full pathway to a bowtie2 index directory, either for spike in or normal alignment
        
        * --reference_file (-r):, no default, must be provided if not bowtie2 indexing files are not already created
    Function output:
        
        * args, object, containing the input data mentioned above
        
    
    """
    
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False, description="Input files for alignment")
    required_name=parser.add_argument_group("required arguments")
    optional_name = parser.add_argument_group("optinal arguments")
    required_name.add_argument('-f','--fastq', help='Input directory name', required=True, type=str)
    optional_name.add_argument("-o", "--out_dir", help="Output directory", required=False, type=str)
    optional_name.add_argument("-i", "--bowtie2_index_pathway", help="Input path to bowtie2 index direcotry", required= False, type=str)
    optional_name.add_argument("-r", "--reference_file", help="Inputh path to reference genome file if no bowtie indexing files are made",required=False, type = str)
    optional_name.add_argument("-s", "--spike_in", help="True if the data is spike in alignment", required=False, type = bool, default=False)
    optional_name.add_argument("-m","--merge_replicates", required=False, type=bool, default=False)

    
    
    return parser


        

def merge(files, fastq, fastq_merged): 
        """
        
        Merges technical duplicates, by adding them to one single file.
        
        
        Function input:
            
            * files, list, lists of fastq files going to get merged
            
            * fastq, str, full path to where the fastq files are stored
            
        
        Function output:
            
            * This function only outputs the the fastq directory from input
            
            * Outputs the merged fastq files
            
            
        """
        tmp_files = glob.glob(files)
        
        for file in tmp_files:
            
            tmp_file_name = re.split(r"_\d+[.]",os.path.basename(file))[0]
            
            cmd = "cat " + fastq +"/" +tmp_file_name+"*"+".fastq* > "+ fastq_merged + "/" + tmp_file_name +".fastq" 
            
            exit_code = subprocess.run(cmd,shell=True)
            if not exit_code.returncode ==0:
                print("Error in merging technical replicates")
                exit(1)



def get_names(files):
            
    """
    Creates a list and appends the sample names from input fastq files.
    Remove duplicate sample names, and sorts them.
    
    
    Function input:
        
        * files, list, list of fastq files being aligned
        
                File names must follow this pattern:
                    
                "[sample name/histone]_rep[replication number]_[R1/R2)].fastq"
                
                
        
    Function output:
        
        * sample_names, list, list of sample names from the input files. 
            
    """
    tmp_files = glob.glob(files)
    sample_names = []
    for file in tmp_files :
        tmp_file_name = pl.PurePath(file).name
        samp = re.split(".fastq", tmp_file_name)[0]                       # Finds the names of the samples by splitting the filename 
        sample_names.append(samp)                                         # Attaches each sample name to the sample_names list 
    sample_names = list(dict.fromkeys(sample_names))                      # Removes sample name duplicates
    sample_names.sort()
    return sample_names



def bowtie2_build(ref_file, align):
    """

    """
    indexing =  os.path.basename(ref_file).split(".")[0]
    
    ref = os.path.join(align, indexing)
    
    os.mkdir(ref)
    
    cmd_build = "bowtie2-build " + ref_file + " "+  ref + "/"+indexing
    
    exit_code = subprocess.run(cmd_build,shell=True)
    
    if exit_code.returncode == 0:
        print("Succeeded in creating Bowtie2 indexing files \n Indexing files avalible at: " + ref )
    else: 
        print("Error: Could not create Bowtie2 indexing files. \n Please check reference file.")
        exit(1)
    
    return ref


def bowtie2(fastq, spike_in, ref, projPath, files):
    
    """
    
    Calles for bowtie2 alignment with:
        * End to end
        
        * Very sensitive
        
        * No mixed
    
    
    Outputs sam files to Epimapper/alignment/sam
    
    Outputs bowtie2 txt summary files to either Epimapper/alignment/bowtie2_summary or Epimapper/alignment/bowtie2_summary_spike_in 
    depending on if the alignment is spike in or not
    
    Function input:
        
        * fastq, str, full pathway to directory where the fastq files being aligned are
        
        * ref, str, full pathway to directory with bowtie2 indexing files
        
        * projPath, str,  full pathway to main output directory
        
        * files, list, list with full pathway to fastq files being aligned
        
    
    Function output:
        
        * sam files from each alignment in Epimapper directory
        
    
    """
    indexing_name = os.path.basename(glob.glob(os.path.join(ref,"*.bt2"))[0])
    indexing = re.split(r"[.]",indexing_name)[0]
    sample_names = get_names(files)
    k=0
    for i in range(int(len(sample_names)/2)):
        tmp_name = re.split("_R1", sample_names[k])[0]
        R1 = str(sample_names[k])
        R2 = str(sample_names[k+1])
        cmd="bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 8 -x " + ref + "/"+ indexing+" -1 " + fastq +"/" +R1 +".fastq" +\
            " -2 " +fastq+"/"+R2+".fastq -S " + projPath +"/Epimapper/alignment/sam/"+tmp_name+".sam &> "+ projPath+"/Epimapper/alignment/bowtie2_summary/"+tmp_name+".txt"
        
        cmd_spike_in="bowtie2 --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 8 -x " + ref +"/"+indexing+ " -1 " +  fastq +"/" +R1 +".fastq" +\
            " -2 " +fastq+"/"+R2+".fastq -S " + projPath +"/Epimapper/alignment/sam_spike_in/"+tmp_name+"_spike_in.sam &> "+ projPath+"/Epimapper/alignment/bowtie2_summary_spike_in/"+tmp_name+"_spike_in.txt"
        
        if spike_in:
            
            exit_code = subprocess.run(cmd_spike_in,shell=True)
            if not exit_code.returncode==0:
                print("Error Bowtie2 in spike in alignment.")
                exit(1)
            
        else:
            exit_code =subprocess.run(cmd,shell=True)
            if not exit_code.returncode==0:
                print("Error in Bowtie2 alignment")
                exit(1)
            
        k=k+2
        
        
def plot_summary(summary_tables, table, spike_in):
    
    """
    
    
    Plots statisitcal summary boxplots from bowtie2 alignment, depending if spike in alignment information is avalible.
    
    Saves figures to  Epimapper/summary_tables
    
    Uses subplots, creats figure with size (15,15)
    
    Plots:
        
        1. Sequencing depth per million
        
        2. Alignable fragments, mapped fragments per million
        
        3. Alignment rate (%)
        
        4. Alignent rate to spike in (if avalible information)
    
    Function input:
        
        * summary_tables, str, full path to output directory, where plots are saved
        
        * table, pd.DataFrame, dataframe containing the data to be plotted
        
        * spike_in, bool, if the data is from spike in alignment or not
    
    
    """
    if len(table) != 0:
        sns.set_style("whitegrid")
        sns.set_palette("husl")
        if spike_in:
            fig, axes= plt.subplots(nrows=1, ncols=2, figsize=(15, 15))
            
            plt1 = sns.boxplot(x = table.Sample, y = table.MappedFragments_SpikeIn, hue = table.Sample, dodge=False, ax = axes[0])
            plt1.set_title("Aligneble Fragments", fontsize = 25.0)
            plt1.set_ylabel("Mapped Fragments per Million", fontsize = 25.0)
            plt1.set_xlabel("")
            plt1.legend([],[], frameon=False)
            plt.setp(plt1.xaxis.get_majorticklabels(), rotation=90, fontsize =25)
            
            plt2 = sns.boxplot(x = table.Sample, y = table.AlignmentRate_SpikeIn.str.strip("%").astype(float), hue = table.Sample, dodge = False, ax = axes[1])
            plt2.set_title("Alignment Rate to Spike In", fontsize = 25.0)
            plt2.set_ylabel("% of Mapped Fragments", fontsize = 25.0)
            plt2.set_xlabel("")
            plt2.legend([],[], frameon=False)
            plt.setp(plt2.xaxis.get_majorticklabels(), rotation=90, fontsize =25)
            
            fig.tight_layout(pad=2.0)
            
            fig.savefig(os.path.join(summary_tables, "Sequencing_depth_spikeIn.png"))
    
        elif len(table.columns) == 7: 
        
            fig, axes= plt.subplots(nrows=2, ncols=2, figsize=(15, 15))
            
            plt1 = sns.boxplot( x = table.Sample, y = table.SequencingDepth/1000000, hue = table.Sample, dodge=False, ax=axes[0,0])
            plt1.set_title("Sequencing Depth", fontsize = 25.0)
            plt1.set_ylabel("Sequencing Depth per Million",fontsize = 25.0)
            plt1.set_xlabel("")
            plt.setp(plt1.xaxis.get_majorticklabels(), rotation=90, fontsize =25)
            plt1.legend([],[], frameon=False)
            
            
            plt2 = sns.boxplot(x = table.Sample, y = table.MappedFragments/1000000, hue = table.Sample, dodge=False, ax = axes[0,1])
            plt2.set_title("Aligneble Fragments", fontsize = 25.0)
            plt2.set_xlabel("")
            plt2.set_ylabel("Mapped Fragments per Million", fontsize = 25.0)
            plt.setp(plt2.xaxis.get_majorticklabels(), rotation=90, fontsize =25)
            plt2.legend([],[], frameon=False)
            
            
            plt3 = sns.boxplot(x = table.Sample, y = table.AlignmentRate.str.strip("%").astype(float), hue = table.Sample, dodge=False, ax = axes[1,0])
            plt3.set_title("Alignment Rate",fontsize = 25.0)
            plt3.set_xlabel("")
            plt3.set_ylabel("% of Mapped Fragments",fontsize = 25.0)
            plt.setp(plt3.xaxis.get_majorticklabels(), rotation=90, fontsize =25)
            plt3.legend([],[], frameon=False)
            
            
            plt4 = sns.boxplot(x = table.Sample, y = table.AlignmentRate_SpikeIn.str.strip("%").astype(float), hue = table.Sample, dodge = False, ax = axes[1,1])
            plt4.set_title("Alignment Rate to Spike In",fontsize = 25.0)
            plt4.set_xlabel("")
            plt4.set_ylabel("% of Mapped Fragments",fontsize = 25.0)
            plt4.legend([],[], frameon=False)
            plt.setp(plt4.xaxis.get_majorticklabels(), rotation=90, fontsize =25)
            
            fig.tight_layout(pad=2.0)
            
            
            fig.savefig(os.path.join(summary_tables, "Sequencing_depth.png"))
        
        
        else:
        
            fig, axes= plt.subplots(nrows=1, ncols=3, figsize=(15, 15))
            
            plt1 = sns.boxplot( x = table.Sample, y = table.SequencingDepth/1000000, hue = table.Sample, dodge=False, ax=axes[0])
            plt1.set_title("Sequencing Depth",fontsize = 25.0)
            plt1.set_xlabel("")
            plt1.set_ylabel("Sequencing Depth per Million", fontsize = 25.0)
            plt.setp(plt1.xaxis.get_majorticklabels(), rotation=90, fontsize =25)
            plt1.legend([],[], frameon=False)
            
            
            plt2 = sns.boxplot(x = table.Sample, y = table.MappedFragments/1000000, hue = table.Sample, dodge=False, ax = axes[1])
            plt2.set_title("Aligneble Fragments", fontsize = 25.0)
            plt2.set_xlabel("")
            plt2.set_ylabel("Mapped Fragments per Million", fontsize = 25.0)
            plt.setp(plt2.xaxis.get_majorticklabels(), rotation=90, fontsize =25)
            plt2.legend([],[], frameon=False)
            
            
            plt3 = sns.boxplot(x = table.Sample, y = table.AlignmentRate.str.strip("%").astype(float), hue = table.Sample, dodge=False, ax = axes[2])
            plt3.set_title("Alignment Rate",fontsize = 25.0)
            plt3.set_xlabel("")
            plt3.set_ylabel("% of Mapped Fragments",fontsize = 25.0)
            plt.setp(plt3.xaxis.get_majorticklabels(), rotation=90, fontsize =25)
            plt3.legend([],[], frameon=False)
            
            fig.tight_layout(pad=2.0)
            
            fig.savefig(os.path.join(summary_tables, "Sequencing_depth_ref.png"))



def summary_table(bowtie2_summary, spike_in):
    
    """
    
    
    Prosseses the information from bowtie2 alignment summary txt files and creats a summary table
    
    Saves the summary table to summary_tables directory inside Epimapper directory
    
    
    Function input:
        
        * bowtie2_summay, str, full path to the directory containing bowtie2 summary txt-files
        
        * spike_in, bool, if the summary txt-files are from spike in calibration or not. 
    
    
    Function output:
        
        * summay, table containing informatinf from the bowtie2 summary txt files:
            
            - Columns ref: Sample, Replication, SequencingDepth, MappedFragmanets, AlignmnetRate
            
            - Columns spike in: Sample, Replication, SequencingDepth, MappedFragments_SpikeIn, AlignmentRate_SpikeIn
    
    """

    file_name = os.path.join(bowtie2_summary, "*.txt")
    files = glob.glob(file_name)
    files.sort()
    if spike_in:
        summary= pd.DataFrame( columns = ["Sample","Replication", "SequencingDepth", "MappedFragments_SpikeIn","AlignmentRate_SpikeIn"])
    else:
        summary = pd.DataFrame(columns = ["Sample", "Replication", "SequencingDepth", "MappedFragments","AlignmentRate"])
    for file in files:
        
        sample_name = pl.PurePath(file).name
        Sample = sample_name.split("_")[0]
        rep = "".join(re.findall("\Brep\d+", sample_name))
        
        with open(file, "r") as t:
    
            tbl = [line.split() for line in t.readlines() if not line.startswith("Warning")]
            clean_tbl = [i for i in tbl if i[0][0].isdigit()]
            if not clean_tbl:
                break
            MappedFragNum = int(clean_tbl[3][0])+ int(clean_tbl[4][0])
            alignResults = [Sample, rep, int(clean_tbl[0][0]),MappedFragNum,clean_tbl[5][0]]
        
        t.close()
           
        summary.loc[len(summary)] = alignResults
    return summary
        

        
def make_summary(bowtie2_summary, bowtie2_summary_spike_in, summary_tables, spike_in): 
    """
    
    Runs functions depending on which data is avalible in  directories : bowtie2_summary and bowtie2_summary_spike_in:
        
        * summary_table()
        
        * plot_summary()
        
    
   Function input:
       
       * bowtie2_summary, str, full path to where bowtie2 alignment txt files are present
       
       * bowtie2_summary_spike_in, str, full path to where bowtie2 alignment spike in  txt files are present
       
       * summary_tables, str, full path to outpur directory, where plots and summary table are created.
       
   
    Function output:
        
        * Functions outputs plots and summary table to summary_tables directory 
    
    """
    if os.path.exists(bowtie2_summary) and os.path.exists(bowtie2_summary_spike_in) and len(os.listdir(bowtie2_summary_spike_in))!= 0 and len(os.listdir(bowtie2_summary))!= 0 :
        sum_1 = summary_table(bowtie2_summary, spike_in = False)
        sum_2 = summary_table(bowtie2_summary_spike_in, spike_in=True)
        if len(sum_1)  and len(sum_2) != 0:
            final_summary = pd.merge(sum_1, sum_2, on = ["Sample","SequencingDepth","Replication"])
            
            if final_summary.shape[0] > 1:
                print("Succeeded in creating Bowtie2 alignment summary report \n Report avalible at: " + os.path.join(summary_tables, "bowtie2_alignment_ref_and_spike_in.csv"))
                
                final_summary.to_csv(os.path.join(summary_tables, "bowtie2_alignment_ref_and_spike_in.csv"),index =False)
            
                plot_summary(summary_tables,final_summary, spike_in = False)
                print("Done with plotting. Plots avalible at: "+summary_tables)
            else:
                print("Error in creating summary report")
                exit(1)
        
    elif spike_in:
        
        if len(os.listdir(bowtie2_summary_spike_in))!= 0:
            
            summary = summary_table(bowtie2_summary_spike_in, True)
            if summary.shape[0]>1:
                print("Succeeded in creating Bowtie2 alignment summary report \n Report avalible at: " +os.path.join(summary_tables,"bowtie2_alignment_spike_in.csv"))
                summary.to_csv(os.path.join(summary_tables,"bowtie2_alignment_spike_in.csv"), index =False )
                
                plot_summary(summary_tables,summary,True)
                print("Done with plotting. Plots avalible at: "+summary_tables)
            else:
               print("Error in creating summary report")
               exit(1) 
        
        
    else:
        if len(os.listdir(bowtie2_summary))!= 0:
            summary = summary_table(bowtie2_summary, spike_in = False)
            if summary.shape[0]>1:
                print("Succeeded in creating Bowtie2 alignment summary report \n Report avalible at: "+os.path.join(summary_tables,"bowtie2_alignment_ref.csv"))
                summary.to_csv(os.path.join(summary_tables,"bowtie2_alignment_ref.csv") , index = False)
                
                plot_summary(summary_tables, table = summary, spike_in = False )
                print("Done with plotting. Plots avalible at: "+summary_tables)
            else:
                print("Error in creating summary report")
                exit(1)
    

def check_input(args):
    fastq = args.fastq
    if os.path.exists(fastq):
        
        path = os.path.join(fastq,"*.fastq")
        files = glob.glob(path)
        if files == []:
            print("Chosen fastq directory: " + fastq+" is empty or does not contain any fastq files. \n  Please check your directory or select another one.")
            exit(1)
    else: 
        print("Chosen fastq directory: " + fastq+" does not exist. \n Please select another directory.")
        exit(1)
        
    if args.bowtie2_index_pathway: 
        if  os.path.exists(args.bowtie2_index_pathway):
            files = glob.glob(args.bowtie2_index_pathway+"/*.bt2")
            if  files== []:
                print("Chosen bowtie2 indexing directory: " + args.bowtie2_index_pathway + " is empty or does not contain bowtie2 indexing files. \n Please check your directory or select another one.")
                exit(1)
        else:
            print("Chosen bowtie2 indexing directory: " + args.bowtie2_index_pathway+" does not exist. Please select another directory")
            exit(1)
            
    if args.reference_file:
        if not os.path.exists(args.reference_file):
            print("Chosen reference genome file: " + args.reference_file + " is non-existent. \n Please check your file or select another one.")
            exit(1)
            
        if not os.path.isfile(args.reference_file):
            print("Chosen reference genome file: " + args.reference_file + " is not a file. \n Please check your file or select another one.")
            exit(1)
            
        elif not os.path.getsize(args.reference_file) > 0:
            print("Chosen reference genome file: " + args.reference_file + " is empty. \n Please check your file or select another one.")
            exit(1)
            
    return True


def run(args):
    """
    
    Creates new directories, if not present, in the output directory:
        
        - Epimapper
        
        - Epimapper/alignment
        
        - Epimapper/alignment/bowtie2_summary
        
        - Epimapper/alignment/bowtie2_summary_spike_in (if spike in alignment)
        
        - Epimapper/alignment/sam
        
        - Epimapper/alignment/sam_spike_in (if spike in alignment)
        
        - Epimapper/summary_tables
        
        
        
        
        Function input:
            
            * args, class, containing the input from shell script or command line
            
            
            
    
    Runs functions created above with input data from args:
        
        * bowtie2()
        
        * make_summary()
        
    
    
    """
    
    if args.out_dir:
        
        projPath = args.out_dir
        
        if os.path.exists(args.out_dir) and os.path.isdir(args.out_dir):
            print("Chosen output directory: " +args.out_dir)
            projPath = args.out_dir
        else:
            print("Chosen output directory: " + args.out_dir + " does not exist or is not a directory. \n Please chose another directory or create current one.")
            exit(1)
    else:
        projPath = os.getcwd()
        print("Current output directory is: " + projPath)
            

    
    path=os.path.join(projPath,"Epimapper")

    if not os.path.exists(path):
        os.mkdir(path)
        
        
    summary_tables=os.path.join(path,"summary_tables")

    if not os.path.exists(summary_tables):
        os.mkdir(summary_tables)    
    
    align = os.path.join(path, "alignment")

    if not os.path.exists(align):
        os.mkdir(align)
    bowtie2_summary_spike_in = os.path.join(align,"bowtie2_summary_spike_in")
    bowtie2_summary = os.path.join(align,"bowtie2_summary")
    sam = os.path.join(align,"sam")
    if not os.path.exists(sam):
        os.mkdir(sam)
        
        
   
    if args.spike_in:
        bowtie2_summary = os.path.join(align,"bowtie2_summary")
        
        if  not os.path.exists(bowtie2_summary_spike_in):
            os.mkdir(bowtie2_summary_spike_in)
            
            
        sam_spike_in = os.path.join(align,"sam_spike_in")
        if not os.path.exists(sam_spike_in):
            os.mkdir(sam_spike_in)
        
    else:
        if not os.path.exists(bowtie2_summary):
               os.mkdir(bowtie2_summary)
        
        
        sam = os.path.join(align,"sam")
        if not os.path.exists(sam):
            os.mkdir(sam)


    if check_input(args):
        spike_in = args.spike_in
    
        fastq = args.fastq
        
    
        files = os.path.join(fastq,"*.fastq")
        
        if not args.bowtie2_index_pathway:
            
            ref_file = args.reference_file
            
            ref = bowtie2_build(ref_file, align)
            
        else: ref = args.bowtie2_index_pathway
        
        
        if args.merge_replicates:
            
            ex_folder = fastq.split(os.path.basename(fastq))[0]
            
            name = os.path.basename(fastq) +"_merged"
            
            fastq_merged = os.path.join(ex_folder,name)
            
            if not fastq_merged.exists():
                
                os.mkdir(fastq_merged)
            
            merge(files, fastq, fastq_merged)
            
            fastq = fastq_merged
            
            files = os.path.join(fastq_merged,"*.fastq*")
        
        
        bowtie2(fastq, spike_in, ref, projPath,files)
        
        
        make_summary(bowtie2_summary, bowtie2_summary_spike_in, summary_tables, spike_in)
    
    

if __name__=='__main__':
    
    """
    Exceutes the script if filename is called in terminal/ shell script. 
    
    
    Runs functions:
        
        * run()
    
    """
    
    
    
    args = set_parser(argparse.ArgumentParser('python bowtie2_alignment.py')).parse_args()
    
    run(args)
    
    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
