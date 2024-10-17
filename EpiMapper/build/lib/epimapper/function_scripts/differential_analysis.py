#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
#jbw 2024
import sys
import pandas as pd

import glob

import subprocess

from scipy.stats import zscore
import re

import pathlib as pl

import argparse

import numpy as np

from scipy.stats import ttest_ind, norm

import multiprocessing as mp

from  .dar_scripts.find_dar import find_col_index4sample, parallel_do_DAR

from .dar_scripts.dar_analysis_may import make_pvalue_files, plot_pca4samples, make_genome_annotation_files_by_hmst, pie_plot_of_dars

from .dar_scripts.dar_peaks2genome import annotation

from .dar_scripts.enrichment import enrichment_main

from .dar_scripts.make_region_files import main as mk_region

import matplotlib as mlt

mlt.use('Agg')






def set_parser(parser):
    
    """
    
    Creates and returns variables from terminal inputs
    
    Function input from shell or terminal:
    
    -Required input:
        
        
        * --peaks (-p): The full pathway to the directory containing files from seacr peak calling
        
        * --bedgraph (-bg): The full pathway to the directory containing the bedgraph files
        
        * --chromosome_sizes (-cs): The full pathway to the file containing chromosome size information
        
        * --blacklist (-bl): The full pahtway to the bed file containing blacklisted reagions of the genome
        
        * --list_a (-la): List of samples that will be compared to list b
        
        * --list_b (-lb): List of samples that will be compred to list a
        
        * --normalize (-n): Bool, if the combined signals should be normalized or not. 
        
	-Optional input:
        
        * --out_dir(-o): The full pathway to desired output directory, where Epimapper directory will be made

    Function output:
        
        * args, object, containing the input data mentioned above
    
    
    """
    
    
    
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False, description="Input sam files for fragment length analysis")
    
    required_name=parser.add_argument_group("required arguments")
    
    optional_name = parser.add_argument_group("optional")
    
    required_name.add_argument("-p", "--peaks", required = True, type = str)
   
    #test jbw 23.06
    optional_name.add_argument("-bg", "--bedgraph", required = False, type = str)
    
    required_name.add_argument("-cs", "--chromosome_sizes", required = True, type = str)
    
    required_name.add_argument("-bl", "--blacklist", required = True, type = str)
    
    required_name.add_argument("-la", "--list_a",  nargs='+', required = True)
    
    #test jbw
    optional_name.add_argument("-an", "--annotate", required=False, type =str, default=None, help="Map peaks to various genomic regions , default=None, use True to enable this function")
    
    required_name.add_argument ("-lb", "--list_b", nargs='+', required= True)
    
    required_name.add_argument("-r", "--reference_refFlat", required=True, help = "Path to a reference refFlat genome file.")
    
    optional_name.add_argument("-o", "--out_dir", required=False, type = str)
    
    optional_name.add_argument("-fold", "--fold_enrichment", required=False, type = bool, default=False, help="Use either fold enrichment or bedgraph for differential analysis, default = bedgraph")
    
    optional_name.add_argument("-cut", "--p_value_cutoff", required=False, type=float)

    #test lgg 16.08
    optional_name.add_argument("-tm", "--test_methods", required = False, type = str, default = "ttest", help="there are four methods, ttest, kstest, mannwhitneyu and ranksumtest")
    #test jbw
    #optional_name.add_argument("-n", "--normalize", required=False, type = bool, default = False, help="Whether to normlaize input reads, default= False")
    optional_name.add_argument("-n", "--normalize", required=False, type =str, default = None , help="Whether to normlaize input reads, default= None, use True if need to normalize the data")


    optional_name.add_argument("-X",
                        metavar='',
                        type=int,
                        default=1000,
                        help="the number of upstream bp, TSS, TES, gene. default=1000")
    optional_name.add_argument("-Y",
                        metavar='',
                        type=int,
                        default=1000,
                        help="the number of downstream bp, TSS, TES, gene. default=1000.")
    optional_name.add_argument("-M",
                        metavar='',
                        type=int,
                        default=10000,
                        help="the number of bp from gene start site, 5dist. default=10000.")
    optional_name.add_argument("-N",
                        metavar='',
                        type=int,
                        default=1000000,
                        help="the number of bp from gene start site, 5dist. default=1000000.")
    optional_name.add_argument("-l", "--minIntergenicLen",
                    metavar='',
                    default=2000,
                    type=int,
                    help="minimum intergenic region distance. default=2000")
    
    optional_name.add_argument("-xL","--maxIntergenicLen",
                          metavar='',
                          default=10000000,
                          type=int,
                          help="maximum intergenic region distance, default=10000000")
    optional_name.add_argument("-i", "--intergenicBTGenes",
                    metavar='',
                    default=True,
                    help="intergenic regions is between gene body regions 'True', or between TSS and TES 'False'. default=True")
    
    
    optional_name.add_argument("-e", "--enhancer", help = "Enchaner file bedfile for annotation")
    
    return parser






def make_100_windows(chromosome_sizes, genome_blacklist, LEN, diff_dir):
    
    """
    
    Uses bedtools makewindows and intersect to create a blacklist window file, based on chromosome sizes
    
    
    Function input:
        
        * chromosome_sizes, str, full path to file containing chromosome size information about the genome
        
        * genome_blacklist, str, full path to bed file with blacklist areas of the genome
        
        * LEN, str, string of a number defineing how many basepairs the  blacklist windows should be
        
        * dir_diff, str, full path to where the blacklist window file will be stored¨
        
   
    
   
    """
    
    cmd1  = "bedtools makewindows -g "+ chromosome_sizes + " -w "+LEN+"  > "+ diff_dir+ "/" + LEN+"b.windows.bed"
    
    out_code1 =subprocess.run(cmd1,shell=True)
    
    cmd2 = "bedtools intersect -v -a "+diff_dir+ "/" + LEN+"b.windows.bed" + " -b "+ genome_blacklist+ " > "+diff_dir+"/"+ LEN+".b.windows.BlackListFiltered.bed"
    
    out_code2 = subprocess.run(cmd2,shell=True)
    
    if out_code1.returncode ==0 and out_code2.returncode==0:
        print("Done with - making bp blacklist removed window files")
        
    else: 
        print("Error in making 100bp blacklist removed window files")
        exit(1)
    
    return diff_dir+"/"+ LEN+".b.windows.BlackListFiltered.bed"






def make_master_peak(peak_files, diff_dir, out_combined_files, list_a, list_b):
    
    """
    Takes many bed files containing peak information, sorts them and
    
    creates one master bed file conatining all the singals before merging the 
    
    overlapping signals with bedtools merge. 
    
    
    Function input: 
        
        * peak_files,list, list of strings to the bed files containing peak signals.
        
        * diff_dir, str, full path to output directory containing out_combined_files
        
        * out_combined_files, str, full path to output directory
    
    
    Runs functions:
        
        *sort_bed_file_df()
    
    """
    
    
    master_peak = pd.DataFrame()
    for file in peak_files:
        new_name = os.path.basename(file).split("_")[0]+"_"+os.path.basename(file).split("_")[1]
        if new_name in list_a or new_name in list_b:
            print(file)
            #jbw 2024
            if (pl.Path(file).stat().st_size<1):
              print("Empty in : ", file )
              continue
            else:
              peak_df = pd.read_table(file, header=None, sep = "\t")
              peak_df2 = peak_df[[0, 1, 2, 3]].copy()
            
            if len(master_peak) == 0:
                master_peak = peak_df2.copy()
            else: 
                master_peak = pd.concat([master_peak, peak_df2])
    #jbw 2024 
    if master_peak.shape[0]==0:
           print(new_name)
           print(list_a)
           print(list_b)
           print("Mismatch between input file name and list_a/list_b cause no input peak data was found, I stop!!")
           sys.exit(new_name + " not in "+ " list_a: "+ " ".join(list_a) + " or list_b: " + " ".join(list_b) )
    #end 
          
    master_peak.columns =["chrs", "pos_start", "pos_end", "type"]
    if master_peak.shape[0]  < 1:
        
        print("Peak filenames not in list a or list b. \n Please check your spelling or peak directory.")
        exit(1)
        
    
    columns_name = master_peak.columns
    
    master_peak.to_csv(os.path.join(out_combined_files,"master_peak.bed"), header = False, index = False)
    
    cmd_sort = "sort -k1,1V -k2,2n -k3,3n " + os.path.join(out_combined_files,"master_peak.bed") + " > " + os.path.join(out_combined_files, "master_peak_sorted.bed")
    
    out_sort=subprocess.run(cmd_sort,shell=True)
    
    sorted_masterpeak = pd.read_csv(os.path.join(out_combined_files, "master_peak_sorted.bed"), names=columns_name)
    
    chrm = sorted_masterpeak[sorted_masterpeak["chrs"].str.contains("M")]

    sorted_masterpeak = sorted_masterpeak.drop(list(chrm.index))

    #jbw 2024
    #sorted_masterpeak = sorted_masterpeak.append(chrm)
    sorted_masterpeak = pd.concat([sorted_masterpeak, chrm], ignore_index=True).copy()

    
    sorted_masterpeak["type"] = sorted_masterpeak[["chrs", "pos_start", "pos_end"]].astype(str).T.agg(':'.join).copy()
    
    
    out_file=os.path.join(out_combined_files, "combined_peaks.bed")
    
    sorted_masterpeak.to_csv(out_file,sep='\t',index=False,header=None)
    
    out_file2= out_file.replace('.bed','_merged.bed') 
    
    cmd3 = "bedtools merge  -c 4  -o collapse -i " + out_file+ " > "+ out_file2
    
    out3 = subprocess.run(cmd3,shell=True)
    if out_sort.returncode==0 and out3.returncode==0:
        print("Done with - Making master peak")
    else:
        print("Error in making master peak")
        exit(1)
    
    cmd4 = "rm " + os.path.join(out_combined_files,"*master*")
    subprocess.run(cmd4,shell=True)
    
    return 
 








def map_bg_window(bedgraph, bdg_files, diff_dir, chromosome_sizes, LEN, list_a, list_b):
    
    """
    Used bedtools map to create window bed files based on bedgraph, blacklisted removed genome file, 
    
    and chromosome sizes.
    
    Creates zipped files.
    
    
    Function input:
        
        * bedgraph, str, full pathway to the directory containing bedgraph files
        
        * bdg_files, list, list of strings to the bedgraph files
        
        * diff_dir, str, full pathway to out_put directory
        
        * chromosome_sizes, str, full pathway to file containing chromosome size information
        
        * genome, str, name of genome
        
        * LEN, str, string of number of basepairs the windows should be
        
    
    """
    for f in bdg_files:
    
        tmp_in_file = f
    
        new_name = re.split(".fragments",pl.PurePath(tmp_in_file).name)[0]
        
        if "sorted" in pl.PurePath(tmp_in_file).name:
            
            continue
        
        if new_name in list_a or new_name in list_b:
        
            cmd_sort = "sort -k1,1V -k2,2n -k3,3n " + f +" > " + os.path.join(bedgraph, os.path.basename(f.replace(".bedgraph", "_sorted.bedgraph")))
            
            subprocess.run(cmd_sort,shell=True)
            
            in_pd = pd.read_table(os.path.join(bedgraph, os.path.basename(f.replace(".bedgraph", "_sorted.bedgraph"))), names = ['chrs','pos_start','pos_end','score'])
            
            chrm = in_pd[in_pd["chrs"].str.contains("M")]

            in_pd = in_pd.drop(list(chrm.index))
           
            #jbw 2024
            #new_pd = in_pd.append(chrm)
            new_pd = pd.concat([in_pd, chrm], ignore_index=True).copy()

            new_pd.to_csv(os.path.join(bedgraph,os.path.basename(f.replace(".bedgraph", "_sorted.bedgraph"))), index = False, header = False, sep = "\t")
            
            
    sorted_bdg_files = glob.glob(os.path.join(bedgraph, "*_sorted.bedgraph"))
    
    
    for file in sorted_bdg_files:
        
    
        tmp_in_file=file
        
        new_name = re.split(".fragments",pl.PurePath(tmp_in_file).name)[0]
        
        cmd= 'bedtools map -g ' + chromosome_sizes+ \
                ' -a ' + diff_dir+"/"+ LEN+".b.windows.BlackListFiltered.bed" \
                ' -b '+ tmp_in_file + ' -c 4 -o mean > ' + diff_dir+ "/"+new_name + "_" + LEN+ 'b.windows.BlackListFiltered.bed'
                
        cmd2=  'gzip  -f '+ diff_dir +"/" +new_name + '_'+ LEN + 'b.windows.BlackListFiltered.bed'
        
        
        out1=subprocess.run(cmd,shell=True)
        
        out2=subprocess.run(cmd2,shell=True)
     
        if not out1.returncode==0 and out2.returncode==0:
            print("Error in mapping bedgraph files.")
            exit(1)
            
        
    print("Done with - Mapping bedgraph windows")
 
        








def combine_windows(diff_dir, out_combined_files):
    """
    Combines all window bin signals from multiple samples to one file.
    
    Function input:
        
        * diff_dir, str, full path to where the input files are located
        
        * out_combined_files, full path to where the combined signal file will be outputed
        
        
    
    """
    bed_files = glob.glob(os.path.join(diff_dir, "*100b.windows.BlackListFiltered.bed.gz"))
    
    all_tmp_df=[]
    all_file=[]
 
    for in_f in bed_files:
  
        tmp_df=pd.read_table(in_f,sep = "\t",header=None,compression='gzip').replace('.',0)
        tmp_df.columns=['chrs','pos_start','pos_end','signal']
        tmp_df['id']=tmp_df['chrs'].astype(str)+':'+tmp_df['pos_start'].astype(str)+ ':'+tmp_df['pos_end'].astype(str)
        all_tmp_df.append(tmp_df[['id','signal']].copy())    
        all_file.append(re.split("_100", pl.PurePath(in_f).name)[0])
        

    merged_df=[]
    for in_df in all_tmp_df:
        if len(merged_df)==0:
            merged_df=in_df.copy()
        else:
            merged_df=pd.merge(merged_df,in_df,how='left',on='id')


    merged_df.columns=['id']+all_file
    merged_df.drop(merged_df[merged_df.iloc[:,1:].apply(lambda x: x.astype(float)==0).sum(axis=1)==merged_df.shape[1]-1].index, inplace = True)
    
    merged_df[['chrs','pos_start','pos_end']]=merged_df.id.str.split(':',expand=True)
    
    cols = ["chrs", "pos_start", "pos_end", "id"] + all_file
    
    out_merged_df=merged_df[cols]
    
    out_file='combined_signals_100b.bed.gz'
    
    out_merged_df.to_csv(os.path.join(out_combined_files,out_file),sep='\t',index=False,header=None ,compression='gzip')

    out_file_head=out_file.replace('.bed.gz','.head')
    
    out_merged_df.head(1).to_csv(os.path.join(out_combined_files, out_file_head),sep='\t',index=False)
    
    print("Done with - Combining windows")
    
    return 



def quantile_normalization(df_input):
    """
    Does quantile normalization on a input directory before sorting it.
    
    Function input:
        
        * df_input, dataFrame, input dataframe that is being quantile normalized
        
        
    Function output:
        
        * df, dataFrame, normalized and sortet df_input
        
    
    """
    df = df_input.copy()
    
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    
    
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
        
    return df.copy()


def combine_signal_enrichment(peaks_folder,blacklist_removed_bin, chromosome_sizes, out_combined_files, searchStr1,searchStr2):
    #jbw 2024 
    enrichment_main(peaks_folder,blacklist_removed_bin, chromosome_sizes, out_combined_files, searchStr1,searchStr2)
    
    merged_peaks = os.path.join(out_combined_files, "combined_peaks_merged.bed")
    
    rm_y = "grep -v ^chrY " + merged_peaks+ " > " +merged_peaks.replace(".bed", "_cleaned.bed")
    subprocess.run(rm_y,shell=True)
    
    rm_m = "grep -v ^chrM " + merged_peaks.replace(".bed", "_cleaned.bed") + " > "+ merged_peaks
    subprocess.run(rm_m,shell=True)
    
    rm = "rm " + merged_peaks.replace(".bed", "_cleaned.bed")
    subprocess.run(rm,shell=True)
    
    
    return



def do_normalization(out_combined_files):
    """
    Quantile normalizes and log transforms the combined signals file.
    
    
    Function input:
        
        * out_combined_files, str, full path to diectory containing combined signal file
    
    
    Runs function:
        
        * quantile_normalization()
    
    """
    in_file = os.path.join(out_combined_files,"combined_signals_100b.bed.gz")
    
    in_df=pd.read_csv(in_file,sep='\t',compression='gzip',header=None)
    
    data_df=in_df.iloc[:,4:]

    #Quantile normalize
    
    quant_a_df = quantile_normalization(data_df+1)
    log_quant_a_df= quant_a_df.apply(lambda x: np.log(x))

    log_quant_a_df.corr(method='pearson')

    out_df=pd.concat([in_df.iloc[:,[0,1,2,3]],log_quant_a_df],axis=1).copy()
    out_file=in_file.replace('.bed.gz', '_quantLog.bed.gz')
    out_df.to_csv(out_file,sep='\t',compression='gzip',index=False,header=None)
    
    print("Done with - Normalization")
    return 





def map_peaks_in_wind(out_combined_files, normalize):
    """
    Uses bedtools intersect to find the window bins overlapping with the
    
    quntile normalized combined signals file and stores them to a new file :
        mapped_peaks_100bp_top01_merged.bed
    
    Function input:
        
        * Out_combined_files, str, full path to where the combined signal and combined peak 
        file is present
    
    
    
    
    """
    if normalize:
        combined_signals  = os.path.join(out_combined_files, "combined_signals_100b_quantLog.bed.gz")
    else:
        combined_signals = os.path.join(out_combined_files, "combined_signals_100b.bed.gz")
        
    combined_peaks = os.path.join(out_combined_files, "combined_peaks_merged.bed")
    
    combined_signals_unzip = combined_signals.split(".gz")[0]
    
    cmd1 = "gunzip " + combined_signals
    subprocess.run(cmd1,shell=True)
    
    cmd2 = "bedtools  intersect -wa -wb -loj  -a " + combined_peaks + " -b " + combined_signals_unzip + " > " + out_combined_files + "/" + "mapped_peaks_100bp_merged.bed"
    subprocess.run(cmd2,shell=True)
    
    cmd3 = "gzip -f " + out_combined_files+ "/" + "mapped_peaks_100bp_merged.bed"
    subprocess.run(cmd3,shell=True)
    
    cmd4 = "gzip -f " + combined_signals_unzip
    
    subprocess.run(cmd4,shell=True)
    print("Done with - Map peaks in windows")
    return







def do_dar_analysis(diff_dir, searchStr1, searchStr2, out_combined_files,cutoff,test_methods):
    """
    Uses the imported fdar and dar modules to preform DAR analysis on two groups of samples.
    
    Outputs the p-values files to DAR directory located inside diff_dir
    
    Fuction input:
        
        * diff_dir, str, full path to input directroy
        
        * searchStr1, list, list of strings of the names of the a group of samples 
        
        * searchStr2, list, list of strings of the names of the a group of samples 
        
        * out_combined_files, str, full path to input directory 
    """
    #number of processes and export DAR p.value
    num_of_process=15
    
    export_cutoff=cutoff
    
    #output data folder and names and pval cutoff for plot
    in_peak_mapped_signal_file = os.path.join(out_combined_files, "mapped_peaks_100bp_merged.bed.gz")
    
    in_ar_file= os.path.join(out_combined_files,"combined_peaks_merged.bed")
    
    in_ar_head_file=os.path.join(out_combined_files, "combined_signals_100b.head")
    
    tmp_peakSignal_df=pd.read_csv(in_peak_mapped_signal_file,sep='\t',header=None,compression='gzip')
    #
    #remove rows with empy value such as . 
    total_sample_size=len(searchStr1) + len(searchStr2)
    
    #check missing values in rows
    is_NAN=tmp_peakSignal_df.iloc[:,-total_sample_size:]=='.'
    
    row_has_NAN=is_NAN.any(axis=1)
    
    in_peakSignal_df=tmp_peakSignal_df[~row_has_NAN].copy()
    
    in_ar_df=pd.read_csv(in_ar_file,sep='\t',header=None)
    
    in_ar_head_df=pd.read_csv(in_ar_head_file,sep='\t')
    #
    #make new id name
    in_ar_df['id']=in_ar_df.iloc[:,0].astype(str)+':'+ in_ar_df.iloc[:,1].astype(str)+':' + in_ar_df.iloc[:,2].astype(str)
    in_peakSignal_df['id']=in_peakSignal_df.iloc[:,0].astype(str)+':'+ in_peakSignal_df.iloc[:,1].astype(str)+ ':'+ in_peakSignal_df.iloc[:,2].astype(str)
    in_peakSignal_df=in_peakSignal_df.replace('.',0).copy()
    #
    #DAR finding
    tmp_xStr1_idx, tmp_yStr2_idx =find_col_index4sample(in_ar_head_df, searchStr1, searchStr2)
    
    out_all_pval, out_passed_pval = parallel_do_DAR((num_of_process,in_ar_df,in_peakSignal_df,tmp_xStr1_idx,tmp_yStr2_idx,test_methods))
    #
    out_all_df=pd.DataFrame.from_dict(out_all_pval,orient='index')
    
    out_file2=in_ar_file.replace('.bed','_pvals'+str(export_cutoff)+'.csv')
    
    out_file=in_ar_file.replace('.bed','_pval.csv')
    
    out_file = out_file.replace("out_combined_files", "DAR")
    
    out_file2 = out_file2.replace("out_combined_files", "DAR")
    
    out_all_df.columns=['tval','pval','max_x_group_A','max_y_group_B','num_x_group_A','num_y_group_B' ]
    
    out_all_df.to_csv(out_file,sep='\t',index=True)
    
    out_all_df2=out_all_df[out_all_df.pval<export_cutoff].copy()
    
    out_all_df2.to_csv(out_file2,sep='\t',index=True)
    
    print('Done with - DAR')









def make_pvalue_format(DAR, out_combined_files,cut_off):
    """
    
    Uses imported dar module to create p value formated file.
    
    Function input:
        
        * DAR, str, full path to directory containing p-value file
        
        * out_combined_files, str, full path to directroy containing combined peaks bed file
        
    Function output:
        * out_file, str, full path to the DMR file
        

    """
    in_dar_file = os.path.join(out_combined_files, "combined_peaks_merged.bed")
    in_dar_pval_file = os.path.join(DAR, "combined_peaks_merged_pval.csv")
    pval_cutoff = cut_off
    out_file, merged_df, out_df=make_pvalue_files(in_dar_file,in_dar_pval_file,pval_cutoff)
    print('Done with - Making dar file')
    
    return out_file
  



def make_region_files(diff_dir,reference, chromosome_sizes,X,Y,M,N,intergenic_between_genes, min_intergenic_len, max_intergenic_len,enhancer):
    out_folder=os.path.join(diff_dir,"region_files")
    remove_short = True
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    
    out_region_list = mk_region(reference,  chromosome_sizes, X, Y, M, N, intergenic_between_genes, min_intergenic_len, max_intergenic_len,
             remove_short, out_folder,enhancer)
    
    if not os.path.getsize(out_region_list)>0:
        print("Error in making region files.\nPlease check that your input referece file is in refFlat format.")
        exit(1)
    else:
        print("Done with - making region files")
    return  out_region_list, out_folder

    
    
    

  
def make_genome_files(out_combined_files, out_file,data, sample_names,cutoff):
    """
    Uses the imported dar module to create statistics about where the peaks are located in 
    the genome.
    
    Files are outputed in new "genome" directory in out_combined_files
    
    Function input:
        
        * out_combined_files, str, full path where genome directory will be made
            and the files inputed are located
        

        
        * data, str, full path to directory containing files defining genome areas.
    
    """
    out_gfolder = out_combined_files
 
    in_file= os.path.basename(out_file).split('.bed')[0]
    out_file_name='Epimapper'
    #test jbw
    make_genome_annotation_files_by_hmst(out_gfolder, in_file, 
               out_file_name, data, sample_names,cutoff)
    print('Done with - Genome file')

    return
 



     
def make_pie_plot(out_combined_files, summary_tables,enhancer,cut_off):
    """
    Uses the module dar to create pie plots based on genome files and saves them to
    summary_tables directory
    
    
    
    Function input:
        
        * out_combined_files, str, full path to input directroy
        
        * summary_tables, str, full path to output directory 
        
       
    """
    pval_cutoff = cut_off
    in_folder=os.path.join(out_combined_files, "genome")
    file_name_string='combined*.bed'
    if isinstance(enhancer, str):
        regions=['tss','tes','gene','5dist','enhancers','intergenic']
    else:
        regions=['tss','tes','gene','5dist','intergenic']
    column_idx2pval=8
    out_fig_name=os.path.join(summary_tables,'total_AR_pie.pdf')
    pie_plot_of_dars(in_folder,file_name_string, regions,
    column_idx2pval,pval_cutoff,out_fig_name)
    print('Done with - Pie plot')
   
    return






def make_pca_plot(out_combined_files, summary_tables, searchStr1, searchStr2, normalize):
    """
    Uses the imported dar module to preform PCA analyzis and plots a 3D plot which is
    saved in summary_tables directory
    
    
    Fuction input:
        
        * out_combined_files, str, full path to the input directory containing 
        the files being used in the PCA analysis
            
        * summary_tables, str, full path to the output directory
        
        * searchStr1, list, list of strings of the names of the a group of samples 
        
        * searchStr2, list, list of strings of the names of the a group of samples 
    
    """
    out_data_folder = out_combined_files
    head_file=out_data_folder+'/combined_signals_100b.head'
    if normalize: 
        data_file=out_data_folder+'/combined_signals_100b_quantLog.bed.gz'
    else: 
        data_file = out_data_folder + "/combined_signals_100b.bed.gz"
    #test jbw 06.23    
    out_fig_name=os.path.join(summary_tables,  searchStr1[0].split('_rep')[0]+'_vs_'+searchStr2[0].split('_rep')[0] +'_pca.pdf')
    plot_pca4samples(head_file,data_file, searchStr1,searchStr2,out_fig_name)
    print('Done with - PCA plot')
      
    return






def annotation2genome(diff_dir, DAR, out_combined_files):
    """
    Uses the imported pg module to annotate the signals to the genome regions and create 
    output cvs file with information.
    
    Function input: 
        
        * diff_for, str, full path i input directory
        
        * DAR, str, full path to output directory
        
    
    """
    in_genome_folder = os.path.join(out_combined_files, "genome")
    in_genome_files = glob.glob(os.path.join(in_genome_folder,"*overlap*"))
    annotation(diff_dir, DAR, in_genome_folder, in_genome_files)
    print("Done with - Genome annotaion")
    return





def check_input(args):
    
    if args.peaks is not None:
        
        peak_dir = args.peaks
        if os.path.exists(peak_dir):
             p_files =glob.glob(os.path.join(peak_dir, "*peaks*.bed"))
             if len(p_files) <1:
                 print("Chosen  peak directory: "+peak_dir+" is empty or does not contain any peak bed files. \nPlease check your directory or select another one.")
                 exit(1)
        else:
             print("Chosen peak directory: "+peak_dir+" does not exist. \nPlease select check your path or select another one")
             exit(1)
    else: 
        print("No peak directory selected. \nPlease provide the path to a directory containing bedgraph files in the -p parameter")
        exit(1)
    
    if args.bedgraph is not None:
        
        bedgraph = args.bedgraph
        if os.path.exists(bedgraph):
             b_files =glob.glob(os.path.join(bedgraph, "*.fragments*.bedgraph"))
             if len(b_files) <1:
                 print("Chosen bedgraph directory: "+bedgraph+" is empty or does not contain any bedgraph files. \nPlease check your directory or select another one.")
                 exit(1)
        else:
             print("Chosen bedgraph directory: "+bedgraph+" does not exist. \nPlease select check your path or select another one")
             exit(1)
    else:
        #test jbw 23.06
        bedgraph=None
    #end test

    chromosome_sizes = args.chromosome_sizes
    if  os.path.exists(chromosome_sizes) and os.path.isfile(chromosome_sizes):
        if not os.path.getsize(chromosome_sizes)>0:
            print("Chosen chromosome sizes file: "+chromosome_sizes +" is empty.\n Please check your file or chose another one.")
            exit(1)
    else:
        print("Chosen chromosome sizes file: "+chromosome_sizes +" is not a file or does not exist. \n Please check your file or chose another one.")
        exit(1)
    
    genome_blacklist = args.blacklist
    if os.path.exists(genome_blacklist) and os.path.isfile(genome_blacklist):
      if not os.path.getsize(genome_blacklist)>0:
          print("Chosen genome blacklist file: "+ genome_blacklist+" is empty.\n Please check your file or chose another one.")
          exit(1)
    else:
      print("Chosen genome blacklist file: "+ genome_blacklist+" is not a file or does not exist. \n Please check your file or chose another one.")
      exit(1)

    #test jbw
    if args.normalize is not None:
        #normalize = args.normalize
        normalize = args.normalize.lower() == 'true'
        #if  not isinstance(normalize,str):
        #    print("-n parameter needs either 'True' or 'False' input, please check your spelling.")
        #    exit(1)
    else:
        normalize=False
    

    reference = args.reference_refFlat 
    if  os.path.exists(reference) and os.path.isfile(reference):
        if not os.path.getsize(reference)>0:
            print("Chosen reference refFlat genome file: "+reference +" is empty.\n Please check your file or chose another one.")
            exit(1)
    else:
        print("Chosen reference refFlat genome file: "+reference +" is not a file or does not exist. \n Please check your file or chose another one.")
        exit(1)
    
    if args.p_value_cutoff is not None:
        cut_off = args.p_value_cutoff 
        if not isinstance(cut_off,float):
            print("-cut, --p_value_cutoff needs to be a numeric value.")
            exit(1)
    else:
        cut_off=0.05
    if args.X is not None:   
        X = args.X 
        if not isinstance(X,int):
            print("-X paramameter, the number of upstream bp, TSS, TES, gene. default=1000, needs to be a numeric value.")
            exit(1)
            
    if args.Y is not None:   
        Y = args.Y 
        if not isinstance(Y,int):
            print("-Y paramameter, the number of upstream bp, TSS, TES, gene. default=1000, needs to be a numeric value.")
            exit(1)
    
    if args.N is not None:   
        N = args.N 
        if not isinstance(N,int):
            print("-N paramameter, the number of bp from gene start site, 5dist. default=10000, needs to be a numeric value.")
            exit(1)
            
            
    if args.M is not None:
        M = args.M 
        if not isinstance(M,int):
            print("-M paramameter, the number of downstream bp, TSS, TES, gene. default=1000, needs to be a numeric value.")
            exit(1)
    
    
    if args.M is not None:
        M = args.M 
        if not isinstance(M,int):
            print("-M paramameter, the number of downstream bp, TSS, TES, gene. default=1000, needs to be a numeric value.")
            exit(1)
    
    
    if args.minIntergenicLen is not None:
        minIntergenicLen = args.minIntergenicLen
        if not isinstance(minIntergenicLen,int):
            print("--minIntergenicLen paramameter, minimum intergenic region distance. default=2000, needs to be a numeric value.")
            exit(1)
        
    if args.maxIntergenicLen is not None:
        maxIntergenicLen = args.maxIntergenicLen
        if not isinstance(maxIntergenicLen,int):
            print("--maxIntergenicLen parameter maximum, intergenic region distance, default=10000000, needs to be a numeric value.")
            exit(1)
        
    if args.intergenicBTGenes is not None:
        intergenicBTGenes=args.intergenicBTGenes
        if not isinstance(intergenicBTGenes, bool):
            print("--intergenicBTGenes parameter, intergenic regions is between gene body regions 'True', or between TSS and TES 'False', default=True, needs to be either 'True' or 'False'")
            exit(1)
    if args.enhancer is not None:
        enhancer = args.enhancer
        if  os.path.exists(enhancer) and os.path.isfile(enhancer):
            if not os.path.getsize(chromosome_sizes)>0:
                print("Chosen enhancer file: "+enhancer +" is empty.\n Please check your file or chose another one.")
                exit(1)
        else:
            print("Chosen enhancer file: "+enhancer +" is not a file or does not exist. \n Please check your file or chose another one.")
        
        
    return peak_dir, bedgraph,chromosome_sizes, genome_blacklist, normalize, reference, X,Y,M,N,maxIntergenicLen, minIntergenicLen,intergenicBTGenes,cut_off

  
def run(args):
    """
    
    Creates new directories, if not present, in the output directory:
        
        - Epimapper
        
        - Epimapper/differential_analysis
        
        - Epimapper/differential_analysis/out_combined_files
        
        - Epimapper/differential_analysis/DAR
        
        
    
    Function input:
        
        * args, object, containing the input from shell script or command line
        
        
    Runs functions created above with input data from args:
        
        * make_100_windows()
        
        * make_master_peak()
        
        * map_bg_window()
        
        * combine_windows()
        
        * do_normalization()
        
        * map_peaks_in_wind()
        
        * do_dar_analysis()
        
        * make_dmr_format()
        
        * make_pie_plot()
        
        * make_pca_plot()
        
        * annotation2genome() 
    
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
        os.makedirs(path)
    
    summary_tables=os.path.join(path,"summary_tables")

    if not os.path.exists(summary_tables):
        os.mkdir(summary_tables)  
        
    diff_dir = os.path.join(path, "differential_analysis")
    if not os.path.exists(diff_dir):
        os.mkdir(diff_dir)
    
    out_combined_files = os.path.join(diff_dir, "out_combined_files")
    if not os.path.exists(out_combined_files):
        os.mkdir(out_combined_files)
    
    DAR = os.path.join(diff_dir, "DAR")
    if not os.path.exists(DAR):
        os.mkdir(DAR)
    
    peak_dir, bedgraph,chromosome_sizes, genome_blacklist, normalize, reference, X,Y,M,N,maxIntergenicLen, minIntergenicLen,intergenicBTGenes, cutoff = check_input(args)
   
    #jbw 2024 
    tmp_file_list= glob.glob(os.path.join(peak_dir,"*peaks.*.bed"))
    if len(tmp_file_list)==0:
       tmp_file_list=glob.glob(os.path.join(peak_dir,"*peaks_sorted.bed"))
    #jbw 2024
    peak_files_all = tmp_file_list.copy()
    #peak_files_all = glob.glob(os.path.join(peak_dir,"*peaks*.bed"))
    
    peak_files = [x for x in peak_files_all if "summitRegion" not in x]
    print("all peak files:")
    print(peak_files)

    LEN = "100"
    
    searchStr1 = args.list_a
    
    searchStr2 = args.list_b

    if args.enhancer is not None:
        enhancer = args.enhancer
    else: 
        enhancer = False
        
    blacklist_bin_file = make_100_windows(chromosome_sizes, genome_blacklist, LEN, diff_dir)

    make_master_peak(peak_files, diff_dir, out_combined_files, searchStr1, searchStr2)
    if args.fold_enrichment: 
        print('Use fold enrichment') 
        combine_signal_enrichment(peak_dir, blacklist_bin_file, chromosome_sizes, out_combined_files, searchStr1,searchStr2)
       #test jbw 23.06 
    elif bedgraph is not None:
        #jbw 2024
        print('Use bedgraph file')
        bdg_files = glob.glob(os.path.join(bedgraph, "*.fragments*.bedgraph"))
        map_bg_window(bedgraph, bdg_files, diff_dir, chromosome_sizes, LEN, searchStr1, searchStr2)      
        combine_windows(diff_dir, out_combined_files)
    else:
        print("Neither folder_enrichment nor bedgraph directory selected. \nPlease provide the path to a directory containing bedgraph files in the -bg parameter or use True for folder_enrichment")
        exit(1)
    #end test

    #jbw 2024
    if normalize:
        print('Normalize input data')
        do_normalization(out_combined_files)
    else:
        print('Do not normalize input data')
    
    map_peaks_in_wind(out_combined_files, normalize)
    # 22.08 lgg 
    test_methods = args.test_methods
    
    do_dar_analysis(diff_dir, searchStr1, searchStr2, out_combined_files,cutoff,test_methods)
    
    out_file = make_pvalue_format(DAR, out_combined_files,cutoff)
    
    #test jbw 06.31
    #if args.annotate:
    region_file, data = make_region_files(diff_dir,reference, chromosome_sizes,X,Y,M,N,intergenicBTGenes, minIntergenicLen, maxIntergenicLen,enhancer)
        
    #test jbw 06.23
    #here sample_names only contain the 2 conditions of the two groups samples, respectively
    #sample_names = searchStr1 + searchStr2
    sample_names= [searchStr1[0].split('_rep')[0] , searchStr2[0].split('_rep')[0]]
    #test jbw
    make_genome_files(out_combined_files, out_file, data, sample_names,cutoff)
    #end test

    make_pie_plot(out_combined_files, summary_tables,enhancer,cutoff)
        
    make_pca_plot(out_combined_files, summary_tables, searchStr1, searchStr2, normalize)
    
    #test  lgg
    if args.annotate is not None:
        is_annotation= args.annotate.lower() == 'true'
        if is_annotation:
            print('Do peak annotation ....')
            annotation2genome(diff_dir, DAR, out_combined_files)
        else:
            print('Skip peak annotation to genome ')
    #end test
    return





if __name__=='__main__':
    
    """
    
    Exceutes the script if filename is called in terminal/ shell script. 
    
    
    Runs functions:
        
        * run()
    
    """
    args = set_parser(argparse.ArgumentParser('python differential_analysis.py')).parse_args()
    
    run(args)
    


