#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import pathlib as pl
import re
import subprocess
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mlp
import plotnine as pn
import numpy as np



mlp.use("agg")


def set_parser(parser):
    
    
    """
    
      
    Creates and returns variables from terminal inputs
    
    Function input from shell or terminal:
    
        - Required input:
            * --sam(-s), str, full path to sam directory containing the sam files being analyzed
            
        - Optinal input:
            * --out_dir(-o), str, The full pathway to desired output directory, where Epimapper directory will be made
        
    
        Function output:
            
            * args, object, containing the input data mentioned above
            
    """
    
    
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False, description="Input sam files for fragment length analysis")
    
    required_name=parser.add_argument_group("required arguments")
    
    optional_name = parser.add_argument_group("optinal arguments")
    
    required_name.add_argument("-o", "--out_dir", help = "Input project pathway", required=False, type=str)
    
    optional_name.add_argument("-s", "--sam", help = "Input directory with sam files", required=True, type=str)
    
    
    
    return parser



    
def frag_len(files,sam, fraglen):
    
    """
    
    Extracts the 9th column from the alignment sam file which is the fragment length, and saves them to fragmentLen.txt files
    
    Function input:
        
        * files, list, list with full path of the files being analyzes
        
        * sam, str, full path to directory containting the files being analyzed
        
        * fraglen, str, full path to where fragmentLen.txt files will be stored
    
    """
    
    
    
    
    tmp_files = glob.glob(files)
    
    for file in tmp_files:
        new_name = pl.PurePath(file).name.split(".")[0]
        
        tmp_file_name = pl.PurePath(file).name
        
        cmd = "samtools view -F \
            0x04 "+ sam +"/"+tmp_file_name+" | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk \
                -v OFS='\t' '{print $2, $1/2}' >" +fraglen+"/"+new_name+"_fragmentLen.txt"
                
        subprocess.run(cmd,shell=True)
    




def plot_len(frag_length_summary, summary_tables):
    
    
    """
    
    Plots violin plot and line plot of fragment lenght and count.
    
    
    
    Function input: 
        * frag_length_summary, pd.DataFrame, containing summary informatin from fragmentLen txt files, made in summary() function
        
        * summary_tables, str, full path to output directory, where the summary table and plots are saved
        
        
        
    
    Function out_put:
        
        * Fragment_violin.png, png file, containing violon plots of fragment length and weighted fragment count for each sample
        
        * Fragment_lineplot.png, png file, coitaining line plots of fragment length and count for each sample
        
    
    
    
    """
    
    
    sample_n = len(set(frag_length_summary.Sample))
    fig, plot = (pn.ggplot() + pn.aes(x = frag_length_summary.Sample.to_list(), y = frag_length_summary.Fragment_length.to_list(), weight = frag_length_summary.Weight.to_list(), fill = frag_length_summary.Sample.to_list()) + \
    pn.geom_violin(bw = 5) +\
    pn.scale_y_continuous(breaks = np.arange(0,850,50)) + \
    pn.theme_bw(base_size = 55) + \
    pn.theme(axis_text_x=pn.element_text(rotation=20, hjust=1,size=40),axis_text_y=pn.element_text(size=40),axis_title=pn.element_text(size=50))+ \
    pn.ylab("Fragment Length") + pn.theme(legend_position=(.5, 0.95), legend_direction='horizontal', legend_title=pn.element_blank())  +\
    pn.theme(figure_size=(5*sample_n, 30)) + \
    pn.xlab("")).draw(show=False, return_ggplot=True)
    
        
    
    fig.savefig(os.path.join(summary_tables, "Fragment_length_violin.png") )
    
    
    plt.figure(figsize=(15,8))
    frag_length_summary.reset_index(drop=True,inplace=True)
    plt2 = sns.lineplot(x =  frag_length_summary.Fragment_length.to_list(), y = frag_length_summary.Fragment_count.to_list(), hue = frag_length_summary.Sample,)
    plt2.set_ylabel("Count",fontsize=20)
    plt2.set_xlabel("Fragment Length (bp)",fontsize=20)  
    plt2.tick_params(labelsize=18) #axis font
    plt.setp(plt2.xaxis.get_majorticklabels(), rotation=90)
    plt2.legend(bbox_to_anchor=(1.005, 0.5), loc="center left", borderaxespad=0, fontsize=20)
    plt.tight_layout()
    
    plt.savefig(os.path.join(summary_tables, "Fragment_length_lineplot.png"))


    
    
  






def summary(fraglen,summary_tables):
    
    """
    
      
    Prosseses the information from fragmentLen txt files and creats a summary table
    
    Saves the summary table to summary_tables directory inside Epimapper directory
    
    
    Function input:
        
        * fraglen, str, full path to the directory containing fragmentLen txt-files
        
        * summary_tables, str, full path to output directory, where the summary table and plots are saved
        
        
    
    Runs function:
        
        *plot_len()
    
    
    Function output:
        
        * Fragment_length.cvs, cvs file,  containing summary information about all the samples from the fragmentLen txt-files
          found in summary_tables
        
        
    
    
    """
    
    
    
    len_files = glob.glob(os.path.join(fraglen, "*_fragmentLen.txt"))
    frag_length_summary = pd.DataFrame(columns = ["Fragment_length", "Fragment_count", "Weight", "Rep", "Sample"])
    
    len_files.sort()
    for file in len_files:
        
        sample = pl.PurePath(file).name.split("_fragmentLen.txt")[0]
        
        
        
        rep = "".join(re.findall("\Brep\d+", sample))
        
        tbl = pd.read_table(file, names = ["Fragment_length", "Fragment_count"])
        tbl["Sample"] = sample
        tbl["Weight"] = tbl["Fragment_count"] / sum(tbl["Fragment_count"])
        tbl["Rep"] = rep
        frag_length_summary = pd.concat(objs=[frag_length_summary,tbl])
        
    if frag_length_summary.shape[0] > 0 :  
        print("Succeeded in creating fragment length summary")
        frag_length_summary.to_csv(os.path.join(fraglen, "Fragment_length_all_samples.csv"), index=False)
    else:
        print("Error in creating fragment length summay")
        exit(1)

    
    
    plot_len(frag_length_summary, summary_tables)
    






def run(args):
    
    
    """
    
    Creates new directories, if not present, in the output directory:
        
        - Epimapper
        
        - Epimapper/fragmentLength
        
        - Epimapper/summary_tables
        
    
        
        
        
        
    Function input:
            
        * args, class, containing the input from shell script or command line
            
            
            
    
    Runs functions created above with input data from args:
        
        * fraglen()
        
        * summary()
        
    
        
    
    
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
     
    summary_tables=os.path.join(path,"summary_tables")
    
    if not os.path.exists(summary_tables):
        os.mkdir(summary_tables)    
     
     
    fraglen = os.path.join(path,"FragmentLength")
    if not os.path.exists(fraglen):
        os.mkdir(fraglen)
     
       
    sam = args.sam
    if  os.path.isdir(sam):
        
       files = glob.glob(sam+"/"+"*.sam")
       if files ==[]:
           print("Chosen sam folder: " + sam +" is empty or does not contain any sam files \n Please check your directory or select another one")
           exit(1)
           
       files = os.path.join(sam,"*.sam")
     
       frag_len(files, sam, fraglen)
     
       summary(fraglen, summary_tables)
    
    
    else:
        print("Chosen sam folder: " + sam +" is not a directory or does not exist.\n Please select another directory.")
        exit(1)    
    

    
    
    
if __name__=='__main__':
    
    """
    Exceutes the script if filename is called in terminal/ shell script. 
    
    
    Runs functions:
        
        * run()
    
    """
    
    
    
    args = set_parser(argparse.ArgumentParser('python fragment_length.py')).parse_args()
    
    run(args)    
    
        
        
        
        
        
        
        
    
