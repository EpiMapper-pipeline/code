#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import glob
import pathlib as pl
import re
import argparse
import pandas as pd

import subprocess

def set_parser(parser):
    
    """
    

            
    """
    
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False, description="Input sam files for fragment length analysis")
    
    required_name=parser.add_argument_group("required arguments")
    
    optional_name = parser.add_argument_group("optional")
    
    optional_name.add_argument("-o", "--out_dir", help = "Output directory", required = False, type = str)
    
    required_name.add_argument("-f", "--fastq", help = "Directory with paired-end fastq files for adaptertrimming", required = True, type = str)
    
    required_name.add_argument("-fa", "--forward_adapter", help="Forward adaptor for R1")
    
    required_name.add_argument("-ba", "--backward_adapter", help="Backward adaptor for R2")
    
    optional_name.add_argument("-c","--cores", help="The number of cores being used in the process.", default="8")

    
    return parser



def get_names(tmp_files):
                              
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
    
    
    sample_names = [] 
    for file in tmp_files :
        tmp_file_name = pl.PurePath(file).name
        samp = re.split(".fastq", tmp_file_name)[0]                       # Finds the names of the samples by splitting the filename 
        sample_names.append(samp)                  # Attaches each sample name to the sample_names list 
    sample_names = list(dict.fromkeys(sample_names))            # Removes sample name duplicates
    sample_names.sort()
    
    return sample_names

    

def run_cutadapt(sample_names,fastq,clean_fastq,forward_adapt,back_adapt,cores): 
    k=0  
    
    
    for i in range(int(len(sample_names)/2)):
    
        tmp_name = re.split("_R1", sample_names[k])[0]
        R1 = glob.glob(os.path.join(fastq,str(sample_names[k])+".*"))[0]

        R2 = glob.glob(os.path.join(fastq,str(sample_names[k+1]+".*")))[0]
    
        out_R1= R1.replace(fastq,clean_fastq)
  
     
        out_R2= R2.replace(fastq,clean_fastq)
     
   
        cmd = "cutadapt --pair-adapters -j "+cores+" -a "+ forward_adapt + " -A " +back_adapt + " -o " +out_R1 + " -p "+ out_R2 + " " + R1 +" " + R2
        subprocess.run(cmd, shell=True)
        k=k+2
        






def run(args):
    fastq=args.fastq
    forward_adapt = args.forward_adapter
    back_adapt = args.backward_adapter
    cores= args.cores
    if not args.out_dir:
        clean_fastq = os.path.join(fastq,"clean")
        if not os.path.exists(clean_fastq):
            os.mkdir(clean_fastq)
    else: 
        clean_fastq=args.out_dir
        
    
    files = os.path.join(fastq,"*.fastq*")  
    tmp_files = glob.glob(files)
        
    sample_names =  get_names(tmp_files) 
    run_cutadapt(sample_names,fastq,clean_fastq,forward_adapt,back_adapt,cores)
    




if __name__=='__main__':
    
    """
    
    Exceutes the script if filename is called in terminal/ shell script. 
    
    
    Runs functions:
        
        * run()
    
    """
    
    
    
        
    args = set_parser(argparse.ArgumentParser('python cutadapt.py')).parse_args()
    
    run(args)
 
   
            
