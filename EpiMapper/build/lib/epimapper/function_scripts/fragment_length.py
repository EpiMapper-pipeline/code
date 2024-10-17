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

from scipy import interpolate
from scipy.ndimage import gaussian_filter1d



mlp.use("agg")


def set_parser(parser):
    
    
    """
    
      
    Creates and returns variables from terminal inputs
    
    Function input from shell or terminal:
    
        - Required input:
            * --sam(-s), str, full path to sam directory containing the sam files being analyzed
            
        - Optinal input:
            * --out_dir(-o), str, The full pathway to desired output directory, where Epimapper directory will be made

            * --se_method(-sm), int, The method for calculating the standard error when computing confidence intervals

            * --local_pct(-lp), float, Specify the percentage of local data points to consider when calculating the local standard error
        
    
        Function output:
            
            * args, object, containing the input data mentioned above
            
    """
    
    
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False, description="Input sam files for fragment length analysis")
    
    required_name=parser.add_argument_group("required arguments")
    
    optional_name = parser.add_argument_group("optinal arguments")
    
    required_name.add_argument("-s", "--sam", help = "Input directory with sam files", required=True, type=str)

    optional_name.add_argument("-o", "--out_dir", help = "Input project pathway", required=False, type=str)

    optional_name.add_argument("-sm", "--se_method", help="Select the method for calculating the standard error when computing confidence intervals: 0: Use raw data without smoothing, 1: Use local standard error based on neighboring data points, 2: Use global standard error based on the entire dataset. Default value is 2.", required=False, type=int, default=2)

    optional_name.add_argument("-lp", "--local_pct", help = "Specify the percentage of local data points to consider when calculating the local standard error. This value determines the window size around each data point used to compute the local mean and standard deviation. Default is 0.05, meaning 5% of the total data points.", required=False, type=float, default = 0.05)
    
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
    

def calculate_global_se(y):
   ''' compute mean and standard deviation of all dataset
   '''
   smooth_mean=np.mean(np.array(y))
   smooth_median=np.median(np.array(y))
   #here compute 95% confidence interval 1.96, 98% confidence interval is 2.33, 99% confidence interval is 2.58
   se = np.std(np.array(y))/np.sqrt(len(y))

   return smooth_mean, smooth_median, se

def calculate_local_se(y, local_pct):
    ''' compute mean and standard deviation of each data point
    '''
    se = []
    window_size = int(local_pct * len(y))
    for i in range(len(y)):
        start = max(0, i - window_size)
        end = min(len(y), i + window_size + 1)
        local_window = y[start:end]

        local_mean = np.mean(local_window)
        local_std = np.std(local_window, ddof=1)
        local_se = local_std / np.sqrt(len(local_window))
        se.append(local_se)
    
    return np.array(se)


def smooth_fragment_length_distribute(frag_length_sample, sample_name, local_pct, se_method):
        x = frag_length_sample["Fragment_length"]
        y = frag_length_sample["Fragment_count"]

        ## calculate the data size before interpolation
        delt_length=max(x)-min(x)

        if delt_length> 10000:
                step_size=250
        elif delt_length>5000 and delt_length<=10000:
                step_size=50
        elif delt_length>100 and delt_length <= 5000:
                step_size=10
        elif delt_length>20  and delt_length<=100:
                step_size=2
        else:
                step_size=1

        num_of_data_size=int(delt_length/step_size)
        xnew = np.linspace(min(x), max(x), num_of_data_size)

        ## interpolation
        f_interp = interpolate.interp1d(x, y, kind='linear', fill_value="extrapolate")
        ynew = f_interp(xnew)

        ## gaussian_filter
        ynew = np.array(ynew, dtype=np.float64)
        ynew_smooth = gaussian_filter1d(ynew, sigma=3)

        frag_length_smooth = pd.DataFrame(columns=['Sample', 'Fragment_length_new', 'Fragment_count_interpolate', 'Fragment_count_smooth','Fragment_count_ci_lower' ,'Fragment_count_ci_higher'])
        ###### calculate the confidence_interval based on the whole dataset
        if se_method == 1: # local se  
                se = calculate_local_se(ynew, local_pct)
        elif se_method == 2: # global se
                smooth_mean, smooth_median, se = calculate_global_se(ynew)
        else:
                raise ValueError("Invalid method for standard error calculation. Please choose 1 for whole SE or 2 for local SE.")

        confidence_interval = 1.96 * se

        frag_length_smooth['Fragment_length_new'] = xnew
        frag_length_smooth['Fragment_count_interpolate'] = ynew
        frag_length_smooth['Fragment_count_smooth'] = ynew_smooth
        frag_length_smooth['Fragment_count_ci_lower'] = ynew_smooth - confidence_interval
        frag_length_smooth['Fragment_count_ci_higher'] = ynew_smooth + confidence_interval
        frag_length_smooth['Sample'] = sample_name
        
        return frag_length_smooth


def plot_len_violin(frag_length_summary, summary_tables):
    """
    Plots violin plot of fragment lenght and count.
    
    Function input: 
        * frag_length_summary, pd.DataFrame, containing summary informatin from fragmentLen txt files, made in summary() function
        * summary_tables, str, full path to output directory, where the summary table and plots are saved

    Function out_put:
        * Fragment_violin.png, png file, containing violon plots of fragment length and weighted fragment count for each sample
    """
    sample_n = len(set(frag_length_summary.Sample))
    #jbw 2024
    (fig, plot) = (pn.ggplot() + pn.aes(x = frag_length_summary.Sample.to_list(), y = frag_length_summary.Fragment_length.to_list(), weight = frag_length_summary.Weight.to_list(), fill = frag_length_summary.Sample.to_list()) + \
    pn.geom_violin(bw = 5) +\
    pn.scale_y_continuous(breaks = np.arange(0,850,50)) + \
    pn.theme_bw(base_size = 60) + \
    pn.theme(axis_text_x=pn.element_text(rotation=25, hjust=1,size=45),axis_text_y=pn.element_text(size=45),axis_title=pn.element_text(size=50),axis_title_y=pn.element_text(size=60,weight='bold'))+ \
    pn.ylab("Fragment Length") + pn.theme(legend_position=(.5, 0.95), legend_direction='horizontal', legend_title=pn.element_blank(),legend_text=pn.element_text(size=50))  +\
    pn.theme(figure_size=(30, 33)) + \
    pn.xlab("")).draw(show=False, return_ggplot=True)
        
    
    fig.savefig(os.path.join(summary_tables, "Fragment_length_violin.png"), dpi=300)
    
def plot_len_lineplot(frag_length_summary, summary_tables):
    """
    Plots line plot of fragment lenght and count based on the raw data.

    Function input: 
        * frag_length_summary, pd.DataFrame, containing summary informatin from fragmentLen txt files, made in summary() function
        * summary_tables, str, full path to output directory, where the summary table and plots are saved

    Function out_put:
        * Fragment_lineplot.png, png file, coitaining line plots of fragment length and count for each sample
    """
    plt.figure(figsize=(15,8))
    frag_length_summary.reset_index(drop=True,inplace=True)
    plt2 = sns.lineplot(x =  frag_length_summary.Fragment_length.to_list(), y = frag_length_summary.Fragment_count.to_list(), hue = frag_length_summary.Sample,)
    plt2.set_ylabel("Count",fontsize=22, fontweight='bold')
    plt2.set_xlabel("Fragment Length (bp)",fontsize=22, fontweight='bold')  
    plt2.tick_params(labelsize=18) #axis font
    plt.setp(plt2.xaxis.get_majorticklabels(), rotation=90)
    plt2.legend(bbox_to_anchor=(1.005, 0.5), loc="center left", borderaxespad=0, fontsize=20)
    plt.tight_layout()
    
    plt.savefig(os.path.join(summary_tables, "Fragment_length_lineplot.png"), dpi=300)


def plot_len_lineplot_smooth(frag_length_smooth, summary_tables):
    """
    Plots line plot of fragment lenght and count based on the smoothing data.

    Function input: 
        * frag_length_smooth, pd.DataFrame, containing infromation after smoothing
        * summary_tables, str, full path to output directory, where the summary table and plots are saved

    Function out_put:
        * Fragment_lineplot.png, png file, coitaining line plots of fragment length and count for each sample
    """
    plt.figure(figsize=(15,8))
    plt2 = sns.lineplot(x='Fragment_length_new', y='Fragment_count_smooth', hue='Sample', data=frag_length_smooth)

    sample_n = frag_length_smooth['Sample'].unique()

    for sample in sample_n:
        sample_data = frag_length_smooth[frag_length_smooth['Sample'] == sample]
        color = sns.color_palette("tab10")[list(sample_n).index(sample)]
        plt.fill_between(sample_data['Fragment_length_new'], 
                        sample_data['Fragment_count_ci_lower'], 
                        sample_data['Fragment_count_ci_higher'], 
                        color=color, alpha=0.2)
        
    plt2.set_ylabel("Count",fontsize=22, fontweight='bold')
    plt2.set_xlabel("Fragment Length (bp)",fontsize=22, fontweight='bold')  
    plt2.tick_params(labelsize=18) #axis font
    plt.title("Smoothing Fragment Length Distribution with 95% CI",fontsize=22, fontweight='bold') #axis font
    plt2.legend(bbox_to_anchor=(1.005, 0.5), loc="center left", borderaxespad=0, fontsize=20)
    plt.tight_layout()

    plt.savefig(os.path.join(summary_tables, "Fragment_length_lineplot.png"), dpi=300)



def summary(fraglen, summary_tables, local_pct, se_method):
    
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

    plot_len_violin(frag_length_summary, summary_tables)

    if se_method == 0:
        print("Plot the original fragment length distribution with confidence intervals")
        plot_len_lineplot(frag_length_summary, summary_tables)

    else:
        print("Plot the smoothed fragment length distribution with confidence intervals")
        frag_length_smooth = pd.DataFrame(columns=['Sample', 'Fragment_length_new', 'Fragment_count_interpolate', 'Fragment_count_smooth','Fragment_count_ci_lower' ,'Fragment_count_ci_higher'])
        for sample in list(frag_length_summary.Sample.unique()):
            frag_length_sample = frag_length_summary[frag_length_summary['Sample'] ==sample][["Fragment_length","Fragment_count"]]
            frag_length_sample.sort_values(by = "Fragment_length", inplace=True)
            frag_length_sample_smooth = smooth_fragment_length_distribute(frag_length_sample, sample, local_pct, se_method)
            if frag_length_smooth.empty:
                frag_length_smooth = frag_length_sample_smooth
            else:
                frag_length_smooth = pd.concat([frag_length_smooth, frag_length_sample_smooth])

        plot_len_lineplot_smooth(frag_length_smooth,summary_tables)
    
    
    






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
    se_method = args.se_method
    local_pct = args.local_pct
    if  os.path.isdir(sam):
        
       files = glob.glob(sam+"/"+"*.sam")
       if files ==[]:
           print("Chosen sam folder: " + sam +" is empty or does not contain any sam files \n Please check your directory or select another one")
           exit(1)
           
       files = os.path.join(sam,"*.sam")
     
       frag_len(files, sam, fraglen)
     
       summary(fraglen, summary_tables, local_pct, se_method)

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
    
        
        
        
        
        
        
        
    
