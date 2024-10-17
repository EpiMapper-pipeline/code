#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import os
import glob
import pathlib as pl
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mlp
import shutil
import subprocess
mlp.use("agg")



def set_parser(parser):
    
    """

    Creates and returns variables from terminal inputs
  
    Function input from shell or terminal:
        
        
        
    Required input:
        
        * --sam (-s), str, full pathway to directory with sam files



    Optional input: 
        
        * --out_dir (-o), full pathway to chosen output directory, default is current directory
      
    
    Function output:
        
        * args, object, containing the input data mentioned above
        
         
    """    
    


    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False, description="remove duplicates")
  
    required_name=parser.add_argument_group("required arguments")
    
    optional_name = parser.add_argument_group("optinal arguments")
    
    optional_name.add_argument("-o", "--out_dir", help="Input project pathway", required= False, type=str)
  
    required_name.add_argument('-s','--sam', help='Input directory pathway', required=True, type=str )

    
  
    return parser

 

    
    
    

def picard(files, sam, sorted_sam, removeDuplicate, picard_summary):
    
    
    
    """
    
    
     Uses every file in input sam directory that ends with .sam
     
     1. Sorts the sam files
     
     2. Markes the duplicates
     
     3. Removes the duplicates and makes picard summary txt files
     
     
     Function input:
         
         * files, list, list of files being analyzed
         
         * sam, str, full pathway to directory with sam files
         
         * sorted_sam, str, full pathway to directory where sorted sam files will be outputed
         
         * removeDuplicate, str, full pathway to directory where output from picard will be put
         
         * picard_summary, str, full pathway to directory where picard summary txt files will be put
         
     
    """
    
    files = os.path.join(sam,"*.sam")
    tmp_files = glob.glob(files)
    #test jbw
    TMP_folder=os.path.join(sam,'TMP')
    TMP_folder=os.path.abspath(TMP_folder)
    if not os.path.exists(TMP_folder):
        os.makedirs(TMP_folder)
        print('Create temporary folder, ', TMP_folder)
    #end test

    for file in tmp_files:
        
        new_name = pl.PurePath(file).name.split(".")[0]
        tmp_file_name = pl.PurePath(file).name
        
        #test jbw
        #cmd1 = "picard  -Dpicard.useLegacyParser=false SortSam -I " +sam +"/" +tmp_file_name + \
        #    " -O "+sorted_sam+"/"+new_name+".sorted.sam --SORT_ORDER coordinate"
        ##jbw 2024
        cmd1 = "picard  -Dpicard.useLegacyParser=false SortSam -I " +sam +"/" +tmp_file_name +  " --VALIDATION_STRINGENCY SILENT "\
            " -O "+sorted_sam+"/"+new_name+".sorted.sam --SORT_ORDER coordinate " + " -TMP_DIR " + TMP_folder   
     
        #end test

        #cmd2 = "picard -Dpicard.useLegacyParser=false  MarkDuplicates -I "+sorted_sam+"/"+new_name+".sorted.sam  \
            #-O "+removeDuplicate+"/"+new_name+".sorted.dupMarked.sam \
               # -METRICS_FILE "+picard_summary+"/"+new_name+"_picard.dupMark.txt"
                
        #test jbw 14.06        
        #cmd3 = "picard -Dpicard.useLegacyParser=false  MarkDuplicates -I "+sorted_sam+"/"+new_name+".sorted.sam \
        #    -O "+removeDuplicate+"/"+new_name+".sorted.rmDup.sam \
        #        -REMOVE_DUPLICATES true -METRICS_FILE "+picard_summary+"/"+new_name+"_picard.rmDup.txt"

        #jbw 2024
        cmd3 = "picard -Dpicard.useLegacyParser=false  MarkDuplicates -I "+sorted_sam+"/"+new_name+".sorted.sam \
            -O "+removeDuplicate+"/"+new_name+".rmDup.sam \
            -REMOVE_DUPLICATES true -METRICS_FILE "+picard_summary+"/"+new_name+"_picard.rmDup.txt"+ " --VALIDATION_STRINGENCY SILENT "

        #jbw 2024
        cmd4= "picard  -Dpicard.useLegacyParser=false SortSam -I " + removeDuplicate +"/" + new_name + ".rmDup.sam" + " --VALIDATION_STRINGENCY SILENT " \
            " -O "+ removeDuplicate + "/"+new_name+".sorted.rmDup.sam --SORT_ORDER queryname " + " -TMP_DIR " + TMP_folder
       
        exit_code1 = subprocess.run(cmd1,shell=True)
        #exit_code2 =subprocess.run(cmd2,shell=True)
        exit_code3 =subprocess.run(cmd3,shell=True)
        exit_code4= subprocess.run(cmd4,shell=True)
        #end test


    
def plot_summary(duplication_summary, summary_tables):
    """
    
    Plots boxplots in one figure based on data from input duplication_summary:
        
        
        1. Duplication rate
        
        2. Estimated library size
        
        3. Number of uniqe fragments
        
        
    
    Function input:
        
        * duplication_summary, pd.DataFrame, dataframe containing summary information from picard duplication removal and make_summary()
        
        * summary_tables, str, full pathway to output directory where plots will be put
    
    
    
    Function output:
        
        * one png file in summary_tables contaning the three mentioned plots 
        
        
        
    """
    #jbq 2024 fill nan as 0
    duplication_summary.fillna(0,inplace=True) 
    
    sns.set_style("whitegrid")
    sns.set_palette("husl")
    fig, axes= plt.subplots(nrows=1, ncols=3)
     
    plt1 = sns.boxplot( x = duplication_summary.Sample, y = duplication_summary.DuplicationRate, hue = duplication_summary.Sample, dodge=False, ax=axes[0])
    plt1.set_ylabel("Duplication Rate (*100%)",fontsize = 11, fontweight='bold')
    plt1.set_xlabel("Sample",fontsize = 10, fontweight='bold')
    plt1.tick_params(axis='x', which='major', labelsize=9)  
    plt1.tick_params(axis='y', which='major', labelsize=7) 
    plt.setp(plt1.xaxis.get_majorticklabels(), rotation=90)
    plt1.legend([],[], frameon=False)
     
    plt2 = sns.boxplot(x = duplication_summary.Sample, y = duplication_summary.EstimatedLibrarySize, hue = duplication_summary.Sample, dodge=False, ax = axes[1])
    plt2.set_ylabel("Estimated Library Size",fontsize = 11, fontweight='bold')
    plt2.set_xlabel("Sample",fontsize = 10, fontweight='bold')
    plt2.tick_params(axis='x', which='major', labelsize=9)  
    plt2.tick_params(axis='y', which='major', labelsize=7) 
    plt.setp(plt2.xaxis.get_majorticklabels(), rotation=90)
    plt2.legend([],[], frameon=False)
     
     
    plt3 = sns.boxplot(x = duplication_summary.Sample, y = duplication_summary.UniqueFragNumber, hue = duplication_summary.Sample, dodge=False, ax = axes[2])
    plt3.set_ylabel("Number of Uniqe Fragments",fontsize = 11, fontweight='bold')
    plt3.set_xlabel("Sample",fontsize = 11, fontweight='bold')
    plt3.tick_params(axis='x', which='major', labelsize=9)  
    plt3.tick_params(axis='y', which='major', labelsize=7) 
    plt.setp(plt3.xaxis.get_majorticklabels(), rotation=90)
    plt3.legend([],[], frameon=False)

    fig.tight_layout(pad=3.0)
     
    fig.savefig(os.path.join(summary_tables, "Duplication_rate.png"), dpi=300)
    
    print("Done with plotting. Plots avalible at: "+summary_tables)





def make_summary(picard_summary,summary_tables):
    """
    
    Goes thorugh the picard summaries and collects for each sample:
        
        * Duplication rate
        
        * Estimated library size
        
        * Unique fragments mapped
    
    Returns boxplots of these data in summary_tables directory
    
    Function input:
        
        * picard_summary, str, full path to where picard summary txt files are saved
        
        * summary_tables, str, full path to output directory where Duplication_summary summary table and plots will be saved 
        
    
    
    Runs functions:
        
        * plot_summary() with pd.DataFrame created as input
    
    
    Function output:
        
        * summary_tables/Duplication_summary.csv, summary table containting information from picard duplication removal
        
        * summary_tabels/Duplicatin_rate.png, plot output from plot_summary(), containing boxplots
    
    
    """
    p_files = os.path.join(picard_summary,"*.rmDup.txt")
    tmp_files = glob.glob(p_files)
    tmp_files.sort()
     
    duplication_summary = pd.DataFrame(columns = ["Sample", "Replicate", "MappedFragmets_ref", "DuplicationRate" ,"EstimatedLibrarySize", "UniqueFragNumber"])
    for file in tmp_files:
    

        sample_name = pl.PurePath(file).name
        
        sample = sample_name.split("_")[0]
        
        rep = "".join(re.findall("\Brep\d+", sample_name))
        print(file) 
        tbl = pd.read_table(file, header = 5, nrows = 1, decimal= ",")
        uniquemapped = round(float(tbl.at[0,"READ_PAIRS_EXAMINED"]) * (1-float(tbl.at[0,"PERCENT_DUPLICATION"])))
    
        
        res = [sample, rep, (tbl.at[0,"READ_PAIRS_EXAMINED"]),round(float(tbl.at[0,"PERCENT_DUPLICATION"])*100,3),tbl.at[0,"ESTIMATED_LIBRARY_SIZE"], uniquemapped]
        
        duplication_summary.loc[len(duplication_summary)] = res
        
    duplication_summary.sort_values("Sample")
    
    if duplication_summary.shape[0] > 0:
        print("Succeeded in creating duplication summary. \n Summary avalible at: "+os.path.join(summary_tables, "Duplication_summary.csv"))
        duplication_summary.to_csv(os.path.join(summary_tables, "Duplication_summary.csv"), index=False)
        plot_summary(duplication_summary, summary_tables)
    else:
        print("Error in creating duplication summary.")
        exit(1)






def run(args):
    
    """
    
    Creates new directories, if not present, in the output directory:
       
           - Epimapper
           
           - Epimapper/alignment
           
           - Epimapper/summary_tables
           
           - sam/sorted_sam 
           
           - Epimapper/alignment/removeDuplicate
           
           - Epimapper/alignment/removeDuplicate/picard_summary
           



    Function input:
          
        * args, class, containing the input from shell script or command line
     
    Runs functions created above with input data from args:
        
        * picard()
        
        * make_summary()


    """
    if args.out_dir:                        
        
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
        
    align = os.path.join(path, "alignment")


    summary_tables=os.path.join(path,"summary_tables")

    if not os.path.exists(summary_tables):
        os.mkdir(summary_tables) 

    if not os.path.exists(align):
        os.mkdir(align)
        
        
    removeDuplicate = os.path.join(align, "removeDuplicate")
    if not os.path.exists(removeDuplicate):
        os.mkdir(removeDuplicate)
        
        
    sam = args.sam
    if  os.path.isdir(sam):
        
        sorted_sam = os.path.join(sam, "sorted_sam")
        if not os.path.exists(sorted_sam):
            os.mkdir(sorted_sam)
        
        sam_duplicates_removed = os.path.join(removeDuplicate,"sam_duplicates_removed")
        if not os.path.exists(sam_duplicates_removed ):
            os.mkdir(sam_duplicates_removed)
        
        picard_summary = os.path.join(removeDuplicate, "picard_summary")
        if not os.path.exists(picard_summary):
            os.mkdir(picard_summary)
    
    else:
        print("Chosen sam folder: "+ sam+" is not a directory or does not exist.\n Please select another directory.")
        exit(1)
    
    
    files = os.path.join(sam, "*.sam")
    
    if not glob.glob(files) ==[]:
    
        picard(files, sam, sorted_sam, removeDuplicate, picard_summary)
        
        make_summary(picard_summary, summary_tables)
        
        cmd_mv = "mv " + os.path.join(removeDuplicate,"*sorted.rmDup.sam") + " " + sam_duplicates_removed
        
        subprocess.run(cmd_mv,shell=True)
       
        #test jbw 13.06
        cmd_rm = "rm " + os.path.join(removeDuplicate, "*.sam")
        subprocess.run(cmd_rm,shell=True)
        #cmd_rm2 = "rm " + os.path.join(sam, "*.sam")
        #subprocess.run(cmd_rm2,shell=True)
        #end test
        
    else:
        print("Chosen sam folder: " + sam+" is empty or does not contain any sam files. \n Please select another directory")
        exit(1)


if __name__=='__main__':
    
    """
    Exceutes the script if filename is called in terminal/ shell script. 
    
    
    Runs functions:
        
        * run()
    
    """
    
    
    
    args = set_parser(argparse.ArgumentParser('python remove_duplicates.py')).parse_args()
    
    run(args)
    














    
