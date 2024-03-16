#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import os
import glob
import pathlib as pl
import re 
import subprocess
import shutil

def set_parser(parser):
    
    
    """
    
    Creates and returns variables from terminal inputs
    
    Function input from shell or terminal:
    
    
	- Requied input:
		* --fastq (-f): The full pathway to fastq, filenames must contain either R1 or R2 depending on the read number.

    - Optinal input:
        * --out_dir (-o): The full pathway to disired output directory, default = current directory
        
        
        
     Function output:
         
         * args, object, containing the input data mentioned above
         
 
    """
    
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False, description="Input files for quality control")
    required_name=parser.add_argument_group("required arguments")
    
    optional_name=parser.add_argument_group("optional argument")
    
    required_name.add_argument('-f','--fastq', help='Input directory name', required=True, type=str)
    
    optional_name.add_argument("-o", "--out_dir", help="Output directory", required=False, type=str)

    return parser


      
def run(args):
    """
    Creates new directories, if not present, in the output directory:
        
        - Epimapper
        
        - Epimapper/fastqc
        
    
    Function input:
        
        * args, object, containing the input from shell script or command line
        
        
    Runs functions created above with input data from args:
        
        * get_names()
        
        * run_fastq()
        
    
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
    
    
    fastq = args.fastq         
        
    path=os.path.join(projPath,"Epimapper")       

    if not os.path.exists(path):
        os.mkdir(path)
        
        
    fastqc = os.path.join(path, "fastqc")                  
    
    if not os.path.exists(fastqc):
        os.mkdir(fastqc)
    
    
    if os.path.exists(fastq): 
        files = os.path.join(fastq,"*.fastq*")     
        
        tmp_files = glob.glob(files)
        
        sample_names =  get_names(tmp_files) 
        
        if not tmp_files ==[]:
        
            run_fastqc(sample_names,fastq, fastqc)
        else:
            print("Chosen fastq directory: " + fastq+" is empty or does not contain any fastq files. \n  Please check your directory or select another one.")
            exit(1)
    else:
        print("Chosen fastq directory: " + fastq +" does not exist. \n Please select another directory.")
        exit(1)
        
    
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
    
    
    
    
    
    
def run_fastqc(sample_names,fastq, fastqc):  
                                                   
    """
    
    
    
     Runs through every name in sample_names, finds R1 and R2 with of the same sample
     
     Makes directroy with same name as the current sample name in the fastqc directory
     
     If there is already a directory with current sample name, it will be deleted
     
     Runs fastqc on both R1 and R2 of current sample and puts output in the newly created sample directory 
     
     Expects only two reads for each sample
     
     
     
     
     Function input:
         
         * sample_names, list, list of names of sample present in the fastq directory
         
         * fastq, str, full path to directory where fastq files are present
         
         *fastqc, str, full path to directory where fastqc reports will be outputed
     
        
    """
    k=0  
    for i in range(int(len(sample_names)/2)):
        tmp_name = re.split("_R1", sample_names[k])[0]
        R1 = str(sample_names[k])
        R2 = str(sample_names[k+1])
        sample_dir = os.path.join(fastqc, tmp_name)
        if os.path.exists(sample_dir):
            shutil.rmtree(sample_dir)
        os.mkdir(sample_dir)
        cmd1 = "fastqc -o "+ sample_dir + "  -f fastq " + fastq +"/" +R1 +".fastq*" 
        cmd2 = "fastqc -o "+ sample_dir + "  -f fastq " + fastq +"/" +R2 +".fastq*" 
        exit_code1 = subprocess.run(cmd1,shell=True)
        exit_code2 = subprocess.run(cmd2, shell=True)
        k=k+2
        if not exit_code1.returncode==0 and exit_code2.returncode==0:
            print("Error in running FastQC")
            exit(0)
            
    
    print("Fastqc quality control successful. \n Quality summary reports avalible at: "+ fastqc)
        
if __name__=='__main__':
    
    
    
    """
    Exceutes the script if filename is called in terminal/ shell script. 
    
    
    Runs functions:
        
        * run()
    
    """
    
    
    
    args = set_parser(argparse.ArgumentParser('python fastqc.py')).parse_args()  
    
    run(args)
    
    
        
        