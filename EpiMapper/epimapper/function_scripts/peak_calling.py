#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import glob
import pathlib as pl
import re
import argparse
import pandas as pd

import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
import pkg_resources

sns.set_style("whitegrid")

import numpy as np
import matplotlib as mlp
mlp.use("agg")





def set_parser(parser):
    
    """
    
    
    Creates and returns variables from terminal inputs
    
    Function input from shell or terminal:
    
    -Required input:
        
        * --bedgraph_ex (-bg), full pathway to bedgraph directory with files being analyzed
        
        * --bedgraph_control (-c), full pathway to bedgraph file that will be used as control 
        
        * --seacr_path (-s), full pathway to seacr shell script
        
        * --fragments (-f), full pathway to directory with fragment bed files being used to find FRiPs (Fragment in peaks)
        
        
	-Optional input:
        
        * --out_dir(-o): The full pathway to desired output directory, where Epimapper directory will be made

         *--fragment_table, (-tbl), str, full path to table with information about number of mapped fragments, column names = ["Sample",	"Replication", "SequencingDepth", "MappedFragments", "AlignmentRate", "MappedFragments_SpikeIn",	"AlignmentRate_SpikeIn"]
                          Default: will be collected from summary_tables, made in python script "bowtie2_alignment.py"
    
    Function output:
        
        * args, object, containing the input data mentioned above
    
    
    """
    
    
    
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False, description="Input sam files for fragment length analysis")
    
    required_name=parser.add_argument_group("required arguments")
    
    optional_name = parser.add_argument_group("optional")
    
    optional_name.add_argument("-bg", "--bedgraph", required = False, type = str)
    
    optional_name.add_argument("-c", "--control_index", required=False, help = "Indexes of control files (i.e, 'control' 'igG' ect.", default = False)
    
    optional_name.add_argument("-s", "--seacr_path", required=False, type=str)
    
    required_name.add_argument("-f", "--fragments", required=True, type = str)
    
    optional_name.add_argument("-b", "--bam", required=False, type = str, help ="Bam files for MACS2 peakcalling")
    
    optional_name.add_argument("-p", "--percentage", required = False, type = str, help ="decimal of peaks to be picked out, default 0.01")
    
    optional_name.add_argument("-soft", "--software", required=False, type = str, default = "seacr", help = "Select peak calling software, either seacr or macs2")
    
    optional_name.add_argument("-tbl", "--fragment_table", required = False, type= str)
    
    optional_name.add_argument("-o", "--out_dir", required=False, type = str)
  
    #test jbw
    optional_name.add_argument("-gs", "--genome_size", required=False, type =str, help=" The effective genome size of the organism (hs, mm, ce, dm), default = hs for human genome.", default ="hs")
    
    optional_name.add_argument("-norm","--seacr_norm", required=False, type=str, help="Seacr peakcalling normalization option (norm, non), default= non for Seacr peakcalling", default="non")
    #end test
    
    
    
    return parser

def get_line_num(byte_string):
   """
   
   
   Decodes, strips and splits and extracts first word of string.
   
   
   
   Function input:
       
       * byte_string, str, string that will be decoded, stripped and split
    
    Function output:
        
        * out_string, str, first, decoded word of byte_string
        
        
   """
    
   out_string=byte_string.decode().strip().split(' ')[0]
   
   return out_string




 
def seacr_run(tmp_files, seacr_path, control, seacr, percentage, bedgraph, norm):
    
    """
    
    
    Checks if the file is a control file by checking if the first word in sample name is [IgG, igg, control, ctrl].
    Will not run peak analysis on control files.
    
    Runs SEACR on bedgraph files to find peaks.
    
    Runs first with control file and then top 0.01 peaks.
    
    
    
    
    Function input:
        
        * tmp_files, list, list of paths to bedgraph files
        
        * seacr_path, str, full pathway to seacr shell script
        
        * control_file, str, full pathway to bedgraph file that will be used as control 
        
        * seacr, str, full pathway to where the peak files will be saved
    
    
    """
    if isinstance(control,str):
        controls=[]
        conditional =[]
        
        for file in tmp_files:
            
            
            if control in  pl.PurePath(file).name.split(".")[0]:
                controls.append( pl.PurePath(file).name.split(".")[0])
            
            else:
                conditional.append(pl.PurePath(file).name.split(".")[0])
         
        #test jbw 15.06   
        #print(control, controls, conditional ) 
        for sample in conditional:
            name= sample.split("_rep")[0]
            rep = sample.split("_rep")[1]
            #print(sample, name, rep) 
            file = glob.glob(os.path.join(bedgraph, sample+"*"))[0]
            #print(file)
            #print(name + "_" + control + "_rep" + rep )
            #print(name + "_" + control + "_rep" + str(1))
            tmp_control_file=''
            for ci in controls:
                if 'rep'+rep in ci :
                    tmp_control_file=ci
            print('tmp_ctrl, ', tmp_control_file)

            #test jbw 15.06 
            #if name+"_"+control+"_rep"+rep in controls:
            #    #control_file = glob.glob(os.path.join(bedgraph,name+"_"+control+"_rep"+rep+"*"))[0]
            #    control_files = glob.glob(os.path.join(bedgraph,name+"_"+control+"_rep"+rep+"*"))
            #    control_file=control_files[0]
            #    print(1, control_files)
            #    control_cmd = "bash " +seacr_path+" "+ file+ " " + control_file +" " + norm+ " stringent "+ os.path.join(seacr,"control") + "/"+ sample +"_seacr_control_peaks"
            if len(tmp_control_file)>1:
                control_files=  glob.glob(os.path.join(bedgraph,tmp_control_file+"*"))
                control_file=control_files[0]
                tmp_ctrl= os.path.basename(control_file)
                #print(tmp_ctrl)
                #if ('normalize' in tmp_ctrl):
                #    str2norm='norm'
                #else:
                #    str2norm='non'

                print(1, control_files)
                control_cmd = "bash " +seacr_path+" "+ file+ " " + control_file +" " + norm + " stringent "+ os.path.join(seacr,"control") + "/"+ sample +"_seacr_control_peaks"
            elif any(name + "_" + control + "_rep" + str(i) in controls for i in range(1, 10)):
                #control_file = glob.glob(os.path.join(bedgraph,name+"_"+control+"_rep*"))[0]
                control_files=  glob.glob(os.path.join(bedgraph,name+"_"+control+"_rep*"))
                control_file=control_files[0]
                print(2, control_files)
                #if 'normalize' in os.path.basename(control_file):
                #    str2norm='norm'
                #else:
                #    str2norm='non'
                control_cmd = "bash " +seacr_path+" "+ file+ " " + control_file + " " + norm +  " stringent "+ os.path.join(seacr,"control") + "/"+ sample +"_seacr_control_peaks"
            else:
                #control_file=glob.glob(os.path.join(bedgraph,controls[0]+"*"))[0] 
                control_files=glob.glob(os.path.join(bedgraph,controls[0]+"*"))
                control_file=control_files[0]
                print(3, control_files)
                #if 'normalize' in os.path.basename(control_file):
                #    str2norm='norm'
                #else:
                #    str2norm='non'

                control_cmd = "bash " +seacr_path+" "+ file+ " " + control_file + " " + norm + " stringent "+ os.path.join(seacr,"control") + "/"+ sample +"_seacr_control_peaks"  
            
            #top comd
            top_cmd = "bash "  + seacr_path+ " "+file+ " " + percentage+ " " + "non" + " stringent " + os.path.join(seacr,"top_"+percentage) + "/"+sample+"_seacr_top."+percentage+"_peaks"
            print(control_cmd)
            subprocess.run(control_cmd, shell = True)
            print(top_cmd)
            subprocess.run(top_cmd,shell = True)


    else:
        for file in  tmp_files:
            sample = pl.PurePath(file).name.split(".")[0]
            top_cmd = "bash "  + seacr_path+ " "+file+ " " + percentage+ " " + "non" + " stringent " + os.path.join(seacr,"top_"+percentage) + "/"+sample+"_seacr_top."+percentage+"_peaks"
            subprocess.run(top_cmd,shell = True)







def macs2_run(macs2,peakCalling, bam_dir,control,percentage, g_size,macs2_control,macs2_top,is_percent):
   
    #teset jbw 13.06
    tmp_files = glob.glob(os.path.join(bam_dir,"*.BlackListFiltered.bam"))
    if len(tmp_files)==0:
       tmp_files = glob.glob(os.path.join(bam_dir,"*mapped.bam"))
    #end test
    #print(tmp_files)

    if control:
        controls=[]
        conditional =[]
        for file in tmp_files:
            print(file, "\n") 
            #print(control) 
            if control in  pl.PurePath(file).name.split(".")[0]:
                controls.append(file)
            
            else:
                conditional.append(pl.PurePath(file).name.split(".")[0])

        for sample in conditional:
            name= sample.split("_rep")[0]
            rep = sample.split("_rep")[1]
            
            #test jbw 13.06
            #file = glob.glob(os.path.join(bam_dir, sample+"*"))[0]
            file = glob.glob(os.path.join(bam_dir, sample+"*.BlackListFiltered.bam"))[0]
            if len(file)<1:
                file= glob.glob(os.path.join(bam_dir, sample+"*.mapped.bam"))[0]
            #end test
            #print(os.path.basename(file).split("."))

            name = os.path.basename(file).split(".")[0]
            control_str = " ".join(controls)
            
            #test jbw
            #cmd_macs_con = "macs2 callpeak -t " +file +" -f  BAMPE -g "+g_size +" -c " +control_str + "  -n " + name+"_macs2_control --outdir " + macs2_control
            cmd_macs_con = "macs2 callpeak -B --SPMR -p " +  percentage  + " -t " +file +" -f  BAMPE -g "+g_size +" -c " +control_str + "  -n " + name+"_macs2_control --outdir " + macs2_control

            if is_percent:
                cmd_macs_top = "macs2 callpeak -B --SPMR -p "+percentage+" -t " +file +" -f  BAMPE -g "+ g_size+" -n " + name+"_macs2_top_"+percentage+" --outdir " + macs2_top
            else:
                cmd_macs_top = "macs2 callpeak -B --SPMR -t " +file +" -f  BAMPE -g "+ g_size+" -n " + name+"_macs2_top_peaks --outdir " + macs2_top
            #end test

            subprocess.run(cmd_macs_con,shell = True)
            subprocess.run(cmd_macs_top,shell = True)
        types=["macs2_control","macs2_top"]
    
    else:
        for file in tmp_files:
            sample=pl.PurePath(file).name.split(".")[0]
            name= sample.split("_rep")[0]
            #print(file) 
            #test jbw
            if is_percent:
                cmd_macs_top = "macs2 callpeak -B --SPMR -p "+ percentage+" -t " +file +" -f  BAMPE -g "+ g_size+"  -n " + sample+"_macs2_top. --outdir " + macs2_top
            else:
                cmd_macs_top = "macs2 callpeak -B --SPMR -t " +file +" -f  BAMPE -g "+ g_size+"  -n " + sample+"_macs2_top. --outdir " + macs2_top
            subprocess.run(cmd_macs_top,shell = True)
            #end test
        types=["macs2_top"]
    return types


def macs2_summary(macs2,summary_tables, macs2_control, macs2_top):
    
    files = glob.glob(os.path.join(macs2_top,"*.narrowPeak"))
    if os.path.exists(macs2_control):
        files= files + glob.glob(os.path.join(macs2_control,"*.narrowPeak"))
        
    peak_width = pd.DataFrame(columns = ["Sample", "Replication" ,"peakType", "PeakWidth"])
    
    for file in files:

        cmd = "sort -k1,1V -k2,2n -k3,3 "+ file +" > " + file.replace(".narrowPeak","_sorted.bed")
        subprocess.run(cmd,shell=True)
        
        #test jbw
        tmp_in_file=file.replace(".narrowPeak","_sorted.bed")
        print(tmp_in_file)
        check_file=os.path.getsize(tmp_in_file)
        if check_file>0:
        #
          in_pd = pd.read_csv(file.replace(".narrowPeak","_sorted.bed"),sep = "\t", header=None)
        
          chrm = in_pd[in_pd[0].str.contains("M", na=False)]
    
          in_pd = in_pd.drop(list(chrm.index))
    
          new_pd = in_pd.append(chrm)
    
          new_pd[0] =  new_pd[0].apply(lambda x: 'chr' + str(x) if not str(x).startswith('chr') else str(x))
        
          new_pd.to_csv(file.replace(".narrowPeak", "_sorted.bed"),header=False, index=False, sep="\t")
 
          tmp_width = new_pd.iloc[:,:3]
        
          tmp_width["Sample"] = file.split("/")[-1].split("_")[0]
          tmp_width["Replication"] = file.split("/")[-1].split("_")[1]
          tmp_width["peakType"] = file.split("/")[-1].split("_")[2] +"_"+file.split("/")[-1].split("_")[3]
          tmp_width["PeakWidth"] = tmp_width[2] - tmp_width[1] 
          peak_width = pd.concat(objs=[peak_width,tmp_width])
        else:
          print('Empty file, skip - !', file)
        #end test

    peak_summary = pd.DataFrame(columns = ["Sample","Replication", "peakType", "peakN"])
    
    sorted_files =glob.glob(os.path.join(macs2_top,"*_sorted.bed"))
    
    if os.path.exists(macs2_control):
        
        sorted_files= sorted_files + glob.glob(os.path.join(macs2_control,"*_sorted.bed"))
    
    sample_names =[]
    
    reps = []
    
    types = []
    
    for file in sorted_files:
        
        cmd =" wc -l " + file
        out =subprocess.check_output(cmd,shell=True)
        
        peak_n = re.search(r'\d+', out.decode('utf-8')).group()
        
        sample_name = pl.PurePath(file).name
        
        peak_type = sample_name.split("_")[2]+"_"+sample_name.split("_")[3]
     
        rep = "".join(re.findall("\Brep\d+", sample_name))
         
        sample = sample_name.split("_")[0]
         
        rep = sample_name.split("_")[1]
         
        sample_names.append(sample)
         
        types.append(peak_type)
        
        reps.append(rep)
        
        peak_summary.loc[len(peak_summary)] = [sample,rep,peak_type,peak_n]
        
    peak_summary.to_csv(os.path.join(summary_tables, "peak_summary.csv"), index=False)

    #test jbw
    reps=list(set(reps))
    #end test
    return sorted_files, peak_summary, peak_width, reps



def seacr_summary(seacr, summary_tables,seacr_top,seacr_control):
    """
    Creates summary table from SEACR peak calling.
    
    Creates list of sample names present in seacr directory.

    Function input:
        
        * seacr, str, full pathway to directory with bed files containing peaks
        
        * summary_tables, str, full pathway to directory where summary table will be outputed
        
    Function output: 
        
        * sample_names, list, list of sample names from the input files. 
        
        * reps, list, list of sample replications from the input files. 
    

    """
    
    
    
    
    peak_summary = pd.DataFrame(columns = ["Sample","Replication", "peakType", "peakN" ])
    p_files = glob.glob(os.path.join(seacr_top, "*.stringent.bed"))
    if os.path.exists(seacr_control):
       p_files = p_files + glob.glob(os.path.join(seacr_control,"*.stringent.bed"))

    
    sample_names =[]
    
    reps = []
    
    types = []
    
    for file in p_files:
        
        sample_name = pl.PurePath(file).name
        
        peak_type = sample_name.split("_")[2]+"_"+sample_name.split("_")[3]
    
        rep = "".join(re.findall("\Brep\d+", sample_name))
        
        
        
        sample = sample_name.split("_")[0]
        
        rep = sample_name.split("_")[1]
        
        sample_names.append(sample)
        
        types.append(peak_type)
        reps.append(rep)
        
        tbl = pd.read_table(file, names = ["Chromosome", "Start","End", "Total_signal","Max_signal","Max_signal_region"])
        
        peak_summary.loc[len(peak_summary)] = [sample,rep,peak_type, len(tbl)]
    
    if peak_summary.shape[0] ==0:
        print("Error in creating peak summary report. Exiting.")
        exit(1)
    peak_summary.to_csv(os.path.join(summary_tables, "peak_summary.csv"), index=False)
        
    sample_names = list(dict.fromkeys(sample_names)) 
    
    reps = list(dict.fromkeys(reps))
    types = list(dict.fromkeys(types))
    types.sort()
    sample_names.sort()
    reps.sort()    
    return sample_names, reps, types, peak_summary
 


       
    
def peak_width_seacr(sample_names, seacr, reps,types,seacr_top,seacr_control):
    
    """
    Calculates peak width from bed files containing peaks.
    
    Function input:
        
        * sample_names, list, list of names of the present samples
        
        * seacr, str, full pathway to directory with bed files containing peaks
        
        * reps, list, list of replications of the samples being analyzed
        
    
    Function output:
        
        * peak_width, pandas dataframe, table containing information about: Sample, Replication, peakType, PeakWidth
    
    
    
    """
    peak_width = pd.DataFrame(columns = ["Sample", "Replication" ,"peakType", "PeakWidth"])
    for rep in reps:
        for peak in types:
            for Sample in sample_names:
                if "control" in peak:
                    
                    file = seacr_control+"/"+Sample+"_"+rep+"_"+peak+"_peaks.stringent.bed"
                else:
                    file = seacr_top+"/"+Sample+"_"+rep+"_"+peak+"_peaks.stringent.bed"
                
            
                tbl = pd.read_table(file, names = ["Chromosome", "Start" ,"End", "v4", "v5", "v6"])
      
                
                tbl["Sample"] = Sample
                tbl["Replication"] = rep
                tbl["peakType"] = peak
                tbl["PeakWidth"] = tbl["End"] - tbl["Start"] 
                peak_width = pd.concat(objs=[peak_width,tbl])
                
    if peak_width.shape[0] ==0:
        print("Error in finding peak widths. Exiting.")
        exit(1)
    return peak_width
        
   
    
   
    
   
    
   
    

   
    
    
 
#def bedtools_seacr(sample_names, fragments, summary_tables, seacr,summary_table, types,seacr_top,seacr_control):  
    
    """
    Runs bedtools intersect between 2 replications to calculate peak reproducibility. 
    
        * (Overlapping peaks between the two samples)/(total peaks)*100 = % of peaks reproducibility
    
    Runs bedtools intersect between replication and fragments to caluclate FRagment proportion in Peaks regions (FRiPs).
    
        * (Overlapps between fragment bed file and peak file)/(total mapped fragments) = FRagment proportion in Peaks regions (FRiPs).
    
    
    Creates a summary table containing information about: Sample, Replication, fragmentInPeaks, peakN, peakType, Peak_reprod_numb, Frips, PeakReprodRate
    
    
    Function input:
        
        * sample_names, list, list of names of the present samples
        
        * fragments, str, full pathway to directory with fragment bed files being used to find FRiPs (Fragment in peaks)
        
        * summary_tables, str, full pathway to directory where summary table will be outputed
        
        * seacr, str, full pathway to directory with bed files containing peaks
        
        * sum_tbl, str, full path to table with information about number of mapped fragments, column names = ["Sample",	"Replication", "SequencingDepth", "MappedFragments", "AlignmentRate", "MappedFragments_SpikeIn",	"AlignmentRate_SpikeIn"]
        
    
    Function output:
    
        * fragInPeak_df, pandas dataframe, summary table with information about peaks reproducibility and FRiPs
        
        * Saves fragInPeak_df as "peakReproducibility_report" to summary_tables directory
    
    """
    
    """
    
    record_frag_in_peaks=pd.DataFrame(columns= ['Sample','Replication','fragmentInPeaks','peakN', "peakType"])
    peak_overlap = pd.DataFrame(columns= ["Sample", "Replication", "peakType","peakN" ,"Peak_reprod_numb"])

    

    for peak in types:
        for Sample in sample_names:
            if "control" in  peak:
                rep1 = seacr_control+"/"+Sample+"_rep1_"+peak+"_peaks.stringent.bed"
                rep2 = seacr_control+"/"+Sample+"_rep2_"+peak+"_peaks.stringent.bed"
            else:
                rep1 = seacr_top+"/"+Sample+"_rep1_"+peak+"_peaks.stringent.bed"
                rep2 = seacr_top+"/"+Sample+"_rep2_"+peak+"_peaks.stringent.bed"
                
            
            cmd1='wc ' + rep1
            cmd2='wc ' + rep2
            tmp_peakN1=subprocess.check_output(cmd1,shell=True)
            tmp_peakN2=subprocess.check_output(cmd2,shell=True)     
            tmp_peakN1=get_line_num(tmp_peakN1)
            tmp_peakN2=get_line_num(tmp_peakN2)
            cmd='bedtools intersect -a ' + rep1 + ' -b '+ rep2 +'|wc'
            out=subprocess.check_output(cmd,shell=True)
            tmp_peakReprodNum=get_line_num(out)
            
            peak_overlap.loc[len(peak_overlap)] = [Sample,"rep1",peak,tmp_peakN1, tmp_peakReprodNum]
            
            peak_overlap.loc[len(peak_overlap)] = [Sample,"rep2",peak,tmp_peakN2, tmp_peakReprodNum]
            
            tmp_peakReprodNum=get_line_num(out)
                
            frag_file1 = fragments + "/" + Sample + "_rep1.fragments_sorted.bed" 
            frag_file2 = fragments + "/" + Sample + "_rep2.fragments_sorted.bed" 
            
            cmd='bedtools intersect -a ' + frag_file1 + ' -b ' + rep1 + ' -u |wc'
            out=subprocess.check_output(cmd,shell=True)
            tmp_fragment_in_peak=get_line_num(out)
            cmd2='wc ' + rep1
            out2=subprocess.check_output(cmd2,shell=True)
            tmp_peaks=get_line_num(out2)
            record_frag_in_peaks.loc[len(record_frag_in_peaks)]=[Sample,"rep1",tmp_fragment_in_peak,tmp_peaks, peak]
            
            
            cmd='bedtools intersect -a ' + frag_file2 + ' -b ' + rep2 + ' -u |wc'
            out=subprocess.check_output(cmd,shell=True)
            tmp_fragment_in_peak=get_line_num(out)
            
            cmd2='wc ' + rep2
            out2=subprocess.check_output(cmd2,shell=True)
            tmp_peaks=get_line_num(out2)
            record_frag_in_peaks.loc[len(record_frag_in_peaks)]=[Sample,"rep2",tmp_fragment_in_peak,tmp_peaks, peak]
            
             

    out_df = record_frag_in_peaks.copy()
 
    
    algReport_df = summary_table[["Sample", "Replication","SequencingDepth", "MappedFragments", "AlignmentRate" ] ].copy()
    
    
    fragInPeak_df=algReport_df.merge(out_df,on=['Sample',"Replication"],how='inner').copy()
    fragInPeak_df=fragInPeak_df.merge(peak_overlap,on=['Sample',"Replication", "peakType", "peakN"],how='inner').copy()
    
    fragInPeak_df['Frips']=fragInPeak_df.fragmentInPeaks.astype(int)/fragInPeak_df.MappedFragments.astype(int)*100
    
    fragInPeak_df["PeakReprodRate"] =  fragInPeak_df.Peak_reprod_numb.astype(int)/fragInPeak_df.peakN.astype(int)*100
    
    fragInPeak_df.to_csv(os.path.join(summary_tables, "peakReproducibility_report"),sep='\t',index=False)
    
    
    
    return fragInPeak_df
"""




def bedtools_seacr(sorted_files, seacr, peak_summary, fragments,sum_tbl, seacr_control, seacr_top, reps):
    
    sample_data = {}
    
    for file in sorted_files:
        
        sample_name = file.split("/")[-1].split("_")[0]
        file_name = file.split("/")[-1]
        
        if sample_name not in sample_data:
            
            sample_data[sample_name] = []

        sample_data[sample_name].append(file_name)

    sample_df = pd.DataFrame(sample_data)
    
    overlaps = {}
    if len(reps)> 1:
        for col in sample_df:
    
            for index, value in sample_df[col].iteritems():
    
                if not value in overlaps:
                    overlaps[value] = []
                names = sample_df[sample_df.index != index][col].to_list()
                names_filtered = [x  for x in names if value.split("_")[3] in x]
                if "control" in value.split("_")[3]:
                    files = [os.path.join(seacr_control,x) for x in names_filtered]
                    cmd = "bedtools intersect -a " +os.path.join(seacr_control,value)+" -b  "+ " ".join(files) +" -u |wc -l"
                else:
                    files = [os.path.join(seacr_top,x) for x in names_filtered]
                    cmd = "bedtools intersect -a " +os.path.join(seacr_top,value)+" -b  "+ " ".join(files) +" -u |wc -l"
               
                print(cmd)
                out=subprocess.check_output(cmd,shell=True)
    
                overlaps[value]=out
        int_overlaps = {key: int(value.strip()) for key, value in overlaps.items()}
        overlaps_df = pd.DataFrame(int_overlaps, index=[0])
        peak_summary["PeakReprodRate"] = ""
        
        for col in overlaps_df:
            sample_name = col.split("_seacr")[0].split("_")[0]
            rep = col.split("_")[1]
            peak_type=col.split("_")[2]+"_"+col.split("_")[3]
            if "control" in peak_type:
                wc = "wc -l " +os.path.join(seacr_control,col)
            else: 
                wc = "wc -l " +os.path.join(seacr_top,col)
            out = subprocess.check_output(wc,shell=True)
            match = re.search(r'\d+', out.decode('utf-8'))
            #test jbw 06.28
            #replication_rate = round((overlaps_df[col].item()/ int(match.group())) *100,2)
            num_of_match_group=int(match.group())
            if num_of_match_group>0:
                replication_rate = round((overlaps_df[col].item()/num_of_match_group ) *100,2)
            else:
                replication_rate = 0

            #end test
            peak_summary.loc[(peak_summary["Sample"] == sample_name) & (peak_summary["Replication"]==rep) & (peak_summary["peakType"] ==peak_type),"PeakReprodRate"] = replication_rate
        
    peak_summary["Frips"]= ""
    
    for file in sorted_files:
        
        sample_name = file.split("/")[-1].split("_peaks")[0]
        
        samp = sample_name.split("_")[0]
        rep = sample_name.split("_")[1]
        peak_type =  sample_name.split("_")[2]+ "_"+sample_name.split("_")[3]
        
        
        
        frag_file = glob.glob(fragments+"/"+samp+"_"+rep+".fragments_sorted.bed")
        
        if len(frag_file)==0:
            print("Error - Could not find fragment bed file corresponding to "+ samp+"_"+rep+". Exiting.")
            exit(1)
        else:
            frag_file = frag_file[0]
            
        cmd_frips = "bedtools intersect -a " +frag_file +" -b " + file +" -u |wc -l" 
        
        frip_numb = subprocess.check_output(cmd_frips, shell=True)
        
        match_frip = int(re.search(r'\d+', frip_numb.decode('utf-8')).group())
        
        mapped_frags = sum_tbl.loc[(sum_tbl["Sample"] == samp) & (sum_tbl["Replication"]==rep), "MappedFragments"].item()
        
        frips = round(((match_frip/mapped_frags)*100),2)
        
        peak_summary.loc[(peak_summary["Sample"] == samp) & (peak_summary["Replication"]==rep) &(peak_summary["peakType"]==peak_type),"Frips"] = frips        
        
    return peak_summary






def bedtools_macs2(sorted_files, macs2, peak_summary, fragments,sum_tbl, macs2_control, macs2_top,reps):
    
    sample_data = {}
    
    for file in sorted_files:
        
        sample_name = file.split("/")[-1].split("_")[0]
        file_name = file.split("/")[-1]
        
        if sample_name not in sample_data:
            
            sample_data[sample_name] = []

        sample_data[sample_name].append(file_name)

    sample_df = pd.DataFrame(sample_data)
    
    overlaps = {}
    if len(reps)>1:
        for col in sample_df:
    
            for index, value in sample_df[col].iteritems():
                print(value)
                if not value in overlaps:
                    overlaps[value] = []
                names = sample_df[sample_df.index != index][col].to_list()
                names_filtered = [x  for x in names if value.split("_")[3] in x]
                print(names, names_filtered)
                #test jbw bug here??
                if "control" in value.split("_")[3]:
                    files = [os.path.join(macs2_control,x) for x in names_filtered]
                    cmd = "bedtools intersect -a " +os.path.join(macs2_control,value)+" -b  "+ " ".join(files) +" -u |wc -l"
                else:
                    files = [os.path.join(macs2_top,x) for x in names_filtered]
                    cmd = "bedtools intersect -a " +os.path.join(macs2_top,value)+" -b  "+ " ".join(files) +" -u |wc -l"
               
                print(cmd)
                out=subprocess.check_output(cmd,shell=True)
    
                overlaps[value]=out
        int_overlaps = {key: int(value.strip()) for key, value in overlaps.items()}
        overlaps_df = pd.DataFrame(int_overlaps, index=[0])
        peak_summary["PeakReprodRate"] = ""
        
        for col in overlaps_df:
            sample_name = col.split("_macs2")[0].split("_")[0]
            rep = col.split("_")[1]
            peak_type=col.split("_")[2]+"_"+col.split("_")[3]
            if "control" in peak_type:
                wc = "wc -l " +os.path.join(macs2_control,col)
            else: 
                wc = "wc -l " +os.path.join(macs2_top,col)
            out = subprocess.check_output(wc,shell=True)
            match = re.search(r'\d+', out.decode('utf-8'))
            replication_rate = round((overlaps_df[col].item()/ int(match.group())) *100,2)
            peak_summary.loc[(peak_summary["Sample"] == sample_name) & (peak_summary["Replication"]==rep) & (peak_summary["peakType"] ==peak_type),"PeakReprodRate"] = replication_rate
            
    peak_summary["Frips"]= ""
    
    for file in sorted_files:
        
        sample_name = file.split("/")[-1].split("_peaks")[0]
        
        samp = sample_name.split("_")[0]
        rep = sample_name.split("_")[1]
        peak_type =  sample_name.split("_")[2]+ "_"+sample_name.split("_")[3]
        
        
        
        frag_file = glob.glob(fragments+"/"+samp+"_"+rep+".fragments_sorted.bed")[0]

        cmd_frips = "bedtools intersect -a " +frag_file +" -b " + file +" -u |wc -l" 
        
        frip_numb = subprocess.check_output(cmd_frips, shell=True)
        
        match_frip = int(re.search(r'\d+', frip_numb.decode('utf-8')).group())
        
        mapped_frags = sum_tbl.loc[(sum_tbl["Sample"] == samp) & (sum_tbl["Replication"]==rep), "MappedFragments"].item()
        
        frips = round(((match_frip/mapped_frags)*100),2)
        
        peak_summary.loc[(peak_summary["Sample"] == samp) & (peak_summary["Replication"]==rep) &(peak_summary["peakType"]==peak_type),"Frips"] = frips        
    
    return peak_summary
        
def plot(fragInPeak_df,summary_tables, width,types,reps):
    """
    Creates 4 plots 
    
        1. Boxplot of number of peaks per sample
        
        2. Log2 of peak width per sample
        
        3. % of Fragments in Peaks per sample
        
        4. Peak Reproduction (%)
        
    
    Function input: 
        
        *  fragInPeak_df, pandas dataframe, summary table with information about peaks reproducibility and FRiPs
        
        * sample_names, list, list of names of the present samples
        
        * summary_tables, str, full pathway to directory where plots will be outputed
        
        * peak_width, pandas dataframe, table containing information about: Sample, Replication, peakType, PeakWidth
    
    
    Function output:
        
        * Four plots in summary_tables directory: "Peak_numbers.png", "Peak_width.png", "Frips.png", "Peak_reprod_numb.png"
    
    
    
    """
    fragInPeak_df.Sample=fragInPeak_df.Sample.astype('category')
    fragInPeak_df.peakN=fragInPeak_df.peakN.astype(int)

    if len(types)>1:
        pal = sns.color_palette("husl", len(fragInPeak_df["Sample"].unique()))
    
        ctrl_df=fragInPeak_df[fragInPeak_df['peakType'].str.contains('control')]
        
        
        ctrl_df.reset_index(inplace=True,drop=True)
    
        
    
    
        fragInPeak_df["peakN"] = fragInPeak_df['peakN'].apply(lambda x: float(x))
      
        if fragInPeak_df.shape[0]>0: 
            g = sns.FacetGrid(fragInPeak_df, col = "peakType")
            
            g.map_dataframe(sns.boxplot,x= "Sample",y = "peakN", dodge=False)
            
            g.add_legend()
            g.fig.subplots_adjust(hspace=0.9, wspace=.15)
            g.savefig(os.path.join(summary_tables,"Peak_numbers.png"))
            
            
            
            width["PeakWidth"] = width["PeakWidth"].apply(lambda x: float(x))
            
            width["PeakWidth"] = width["PeakWidth"].apply(lambda x: np.log(x))
            
            
            
            
            s = sns.FacetGrid(width, col = "peakType", row= "Replication")
            
            s.map_dataframe(sns.violinplot,x = "Sample", y = "PeakWidth", hue = "Sample", dodge=False)
            
            s.set_ylabels("Peak width (log10)")
            
            s.add_legend()
            s.fig.subplots_adjust(hspace=0.40, wspace=.80)
            s.savefig(os.path.join(summary_tables,"Peak_width.png"))
            
            
            
            
            
            v = sns.FacetGrid(fragInPeak_df, col = "peakType")
            
            v.map_dataframe(sns.boxplot,x = "Sample", y = "Frips",dodge=False, hue = "Sample")
            
            v.set_xlabels("% of Fragments in Peaks")
            
            v.add_legend()
            
          
            v.savefig(os.path.join(summary_tables,"Frips.png"))
            
            
            
            
            if len(reps)>1:
                l = sns.FacetGrid(fragInPeak_df, col = "peakType", row= "Replication")
                
                l.map_dataframe(sns.barplot,x = "Sample", y = "PeakReprodRate", hue = "Sample", dodge=False)
                
                l.add_legend()
                
                l.set_axis_labels("", "Peak Reproduction (%)")
                
                plt.subplots_adjust(bottom=0.25)
                l.fig.subplots_adjust(hspace=0.4, wspace=1)
                l.savefig(os.path.join(summary_tables,"Peak_reprod_numb.png"))
        

        else: 
            print(fragInPeak_df.shape)
            print(fragInPeak_df)
            
    elif len(types) ==1:
        
        if fragInPeak_df.shape[0]>0: 
            
            

            # Convert PeakWidth to float and apply log transformation
            width["PeakWidth"] = width["PeakWidth"].apply(lambda x: float(x))
            #width["PeakWidth"] = width["PeakWidth"].apply(lambda x: np.log(x))
            
            # Initialize a color palette for different samples
            pal = sns.color_palette("husl", len(fragInPeak_df["Sample"].unique()))
            
            # Create and save each plot in separate files
            
            # Plot 1: Box Plot for Peak Numbers
            plt.figure(figsize=(15, 8))
            plt1 = sns.boxplot(x=fragInPeak_df["Sample"], y=fragInPeak_df["peakN"], palette=pal)
            plt1.set_ylabel("Number of Peaks")
            plt1.set_xlabel("Sample")
            plt1.set_xticklabels(plt1.get_xticklabels(), rotation=90)
            plt1.grid(True)
            plt.tight_layout()
            plt.grid()
            plt1.get_figure().savefig(os.path.join(summary_tables, "Peak_numbers.png"))
            
            # Plot 2: Violin Plot for Peak Width
            plt.figure(figsize=(15, 8))
            plt2 = sns.violinplot(x=width["Sample"], y=width["PeakWidth"], palette=pal, dodge=False)
            plt2.set_ylabel("Peak Width")
            plt2.set_xlabel("Sample")
            plt2.set_xticklabels(plt2.get_xticklabels(), rotation=90)
            plt2.grid(True)
            plt.tight_layout()
            plt.grid()
            plt2.get_figure().savefig(os.path.join(summary_tables, "Peak_width.png"))
            
            # Plot 3: Box Plot for Fragment Proportions (frips)
            plt.figure(figsize=(15, 8))
            plt3 = sns.boxplot(x=fragInPeak_df["Sample"], y=fragInPeak_df["Frips"], palette=pal, dodge=False)
            plt3.set_ylabel("Fragment Proportions in Peaks (%)")
            plt3.set_xlabel("Sample")
            plt3.set_xticklabels(plt3.get_xticklabels(), rotation=90)
            plt3.grid(True)
            plt.tight_layout()
            plt3.get_figure().savefig(os.path.join(summary_tables, "frips.png"))
            
            # Plot 4: Bar Plot for Peak Reproducibility Rate
            if len(reps)>1:
                plt.figure(figsize=(15, 8))
                plt4 = sns.barplot(x=fragInPeak_df["Sample"], y=fragInPeak_df["PeakReprodRate"], palette=pal, dodge=False)
                plt4.set_ylabel("Peak Reproducibility Rate (%)")
                plt4.set_xlabel("Sample")
                plt4.set_xticklabels(plt4.get_xticklabels(), rotation=90)
                plt4.grid(True)
                plt.tight_layout()
                plt4.get_figure().savefig(os.path.join(summary_tables, "peaks_reproducibility_rate.png"))
                
                
    

    
    #test jbw 16.06
def peakcall_seacr(seacr_path,peakCalling, summary_tables,sum_tbl,bedgraph,control, fragments, percentage,skip_plot,norm): 

    seacr = os.path.join(peakCalling, "seacr")
    if not os.path.exists(seacr):
        os.mkdir(seacr)
    seacr_control = os.path.join(seacr,"control")
    
    seacr_top = os.path.join(seacr,"top_"+percentage)
    
    if isinstance(control,str):
        if not os.path.exists(seacr_control):
            os.mkdir(seacr_control)
        if not os.path.exists(seacr_top):
            os.mkdir(seacr_top)
    else:
        seacr_top = os.path.join(seacr,"top_"+percentage)
        if not os.path.exists(seacr_top):
            os.mkdir(seacr_top)
        
    
    tmp_files = glob.glob(os.path.join(bedgraph,"*.bedgraph"))
    
    #test jbw 16.06
    #norm = "non"
    
    if skip_plot:
        print("Performing peak calling without calculating peak reproducibility or generating plots.")
        seacr_run(tmp_files, seacr_path, control, seacr, percentage, bedgraph, norm)
        print("Done with calling peaks using SEACR software.\n Called peaks avalible at: "+ seacr)
        exit(0)
    else:
        seacr_run(tmp_files, seacr_path, control, seacr, percentage, bedgraph, norm)
        sample_names, reps, types, peak_summary = seacr_summary(seacr, summary_tables,seacr_top,seacr_control)
    
        width = peak_width_seacr(sample_names, seacr, reps,types,seacr_top,seacr_control)
        
        p_files = glob.glob(os.path.join(seacr_top,"*.stringent.bed"))
        if isinstance(control, str):
            p_files+= glob.glob(os.path.join(seacr_control,"*.stringent.bed"))
        
        fragInPeak_df = bedtools_seacr(p_files, seacr, peak_summary, fragments,sum_tbl, seacr_control, seacr_top,reps)
    
        plot(fragInPeak_df,summary_tables, width,types,reps)
    
     
    
def peakcall_macs2(peakCalling, bam_dir,control,percentage,summary_tables, fragments,sum_tbl, skip_plot,genome_size,is_percent):
    
    
    
    macs2 = os.path.join(peakCalling, "macs2")
    
    if not os.path.exists(macs2):
        
        os.mkdir(macs2)
    macs2_control=os.path.join(macs2,"control")
    macs2_top = os.path.join(macs2, "top_peaks")
                               
    if isinstance(control, str):
        if not os.path.exists(macs2_control):
            os.mkdir(macs2_control)
        if not os.path.exists(macs2_top):
             os.mkdir(macs2_top)
    else:
         if not os.path.exists(macs2_top):
             os.mkdir(macs2_top)
    
    if skip_plot:
        print("Performing peak calling without calculating peak reproducibility or generating plots.")
        types =macs2_run(macs2,peakCalling, bam_dir,control,percentage, genome_size,macs2_control,macs2_top,is_percent)
        #test jbw 2024 sort peak files for heatmap plot??
        #sorted_files, peak_summary, peak_width,reps = macs2_summary(macs2,summary_tables, macs2_control,macs2_top)
        #end test
        print("Done with calling peaks using MACS2 software.\nCalled peaks avalible at: "+ macs2)
        exit(0)
    else:
        types =macs2_run(macs2,peakCalling, bam_dir,control,percentage,genome_size, macs2_control, macs2_top,is_percent)
        
        sorted_files, peak_summary, peak_width,reps = macs2_summary(macs2,summary_tables, macs2_control,macs2_top)
       
        print('reps -> ',reps)
        full_peak_summary=bedtools_macs2(sorted_files, macs2, peak_summary, fragments,sum_tbl,macs2_control, macs2_top,reps)
       
        plot(full_peak_summary, summary_tables,peak_width, types,reps)



def check_input(args,summary_tables):
    skip_plot=False
    
    if args.percentage is not None:
        percentage = args.percentage
        is_percent = True
        try:
            percentage_value = float(percentage)
            if not 0 <= percentage_value <= 1:
                print("Chosen percentage: " + percentage + " is not within the range [0, 1]. Please choose a value in that range.")
                exit(1)
        except ValueError:
            print("Chosen percentage: " + percentage + " is not a valid numeric value. Please choose a numeric value between [0, 1]")
            exit(1)

    else: 
        percentage="0.01"
        is_percent=False

        
    software = args.software
    
    if software=="seacr": 
            
        if args.bedgraph is not None:
            
            bedgraph = args.bedgraph
            if os.path.exists(bedgraph):
                 b_files =glob.glob(os.path.join(bedgraph, "*.bedgraph"))
                 if len(b_files) <1:
                     print("Chosen bedgraph directory: "+bedgraph+" is empty or does not contain any bedgraph files. \nPlease check your directory or select another one.")
                     exit(1)
            else:
                 print("Chosen bedgraph directory: "+bedgraph+" does not exist. \nPlease select check your path or select another one")
                 exit(1)
        else: 
            print("No bedgraph directory selected. \nPlease provide the path to a directory containing bedgraph files in the -bg parameter")
            exit(1)
    
        
    
    
    elif software=="macs2":
        
        if args.bam is not None:
            bam_dir = args.bam
            if os.path.exists(bam_dir):
                b_files = glob.glob(os.path.join(bam_dir, "*.mapped.bam"))
                if len(b_files) <1:
                    print("Chosen bam directory: "+bam_dir+" is empty or does not contain any mapped.bam files. \n Please check your directory or select another one.")
                    exit(1)
            else:
                print("Chosen bam directory: "+bam_dir+" does not exist. \nPlease select check your path or select another one")
                exit(1)
        else: 
            print("No bam directory selected. \nPlease provide the path to a directory containing bedgraph files in the -b parameter")
            exit(1)
            
        if args.genome_size is not None:
            genome_size= args.genome_size
            #test jbw
            #try:
            #      float(genome_size)
            #      
            #except ValueError:
            #    print("Chosen genome size: "+ genome_size+" is not numeric\nPlease choose provide the effective genome size of the organisme (i.e 1.87e9 for mus musculus")
            #    exit(1)
            #end test
                
    else:
        print("Software not recognized.\nPlease choose between the peak calling software 'seacr' and 'macs2' in the -soft parameter.")
        exit(1)
    
    fragments = args.fragments
    if os.path.exists(fragments):
        files =glob.glob(os.path.join(fragments, "*fragments_sorted.bed"))
        if len(files) <1:
            print("Chosen fragments directory: "+fragments+" is empty or does not contain any fragments_sorted.bed files. \nPlease check your directory or select another one.")
            exit(1)
    else:
        print("Chosen fragments directory: "+fragments+" does not exist. \nPlease select check your path or select another one")
        exit(1)
    
    control = args.control_index
    if isinstance(control, str):
        if not any(control in s for s in b_files):
            print("Chosen control index: "+control+" is not a part of the files names./nPlease check your spelling or inspect your bedgraph file names.")
            
    if args.fragment_table:
        tbl = args.fragment_table
        if os.path.exists(tbl) and os.path.isfile(tbl):
            if os.path.getsize(tbl) == 0:
                print("Chosen summary table file: "+ tbl+" is empty.\nPlease check your file or provide another one.")
                exit(1)
            df = pd.read_csv(tbl,nrows=1)
            column_names = df.columns.tolist()
            if not all( item in  column_names for item in ["Sample",	"Replication", "SequencingDepth", "MappedFragments", "AlignmentRate"]):
                print('Chosen summary table does not contain the following columns:\n ["Sample",	"Replication", "SequencingDepth", "MappedFragments", "AlignmentRate"]. Please update your table´s column names or chose another table')
                exit(1)
            
    else:
        tbl = os.path.join(summary_tables,"bowtie2_alignment_ref.csv")
        
        if not os.path.exists(tbl):
            while True:
                x = input("No summary table provided and non avalible in :"+ summary_tables+". \n:Do you wish to continue peak calling without peak reproducibility report and plot generation? y/n")
                if "y" in x or "yes" in x:
                    
                    print("Continuing without peak reproducibility report and plot generation")
                    skip_plot = True
                    break
                elif "n" in x or "no" in x:
                    print("Exiting.")
                    exit(0)
                else:
                    print("Invalid input. Please type 'y' or 'n'.")  
    return tbl, percentage, fragments, control, skip_plot, software, is_percent
       
def run(args):
    """
    
    Creates new directories, if not present, in the output directory:
        
        - Epimapper
        
        - Epimapper/peakCalling
        
        - Epimapper/peakCalling/seacr
        
        - Epimapper/summary_tables
        
        
        
        
        Function input:
            
            * args, class, containing the input from shell script or command line
            
            
            
    
    Runs functions created above with input data from args:
        
        * seacr_run()
        
        * seacr_summary()
        
        * peak_width()
        
        * bedtools()
        
        * plot()
        
        """

    if args.out_dir:
        
        projPath = args.out_dir
        
        if os.path.exists(args.out_dir) and os.path.isdir(args.out_dir):
            print("Chosen output directory: " +args.out_dir)
            projPath = args.out_dir
        else:
            print("Chosen output directory: " + args.out_dir + " does not exist or is not a directory. \n Creating directroy.")
            os.makedirs(args.out_dir)
            projPath = args.out_dir
    else:
        projPath = os.getcwd()
        print("Current output directory is: " + projPath)
    
        
    path=os.path.join(projPath,"Epimapper")
    if not os.path.exists(path):
        os.makedirs(path)
    
    summary_tables=os.path.join(path,"summary_tables")

    if not os.path.exists(summary_tables):
        os.makedirs(summary_tables) 
    
    
    peakCalling = os.path.join(path,"peakCalling")
    
    if not os.path.exists(peakCalling):
        os.mkdir(peakCalling)
    
        
          
    tbl, percentage, fragments, control, skip_plot, software,is_percent =check_input(args,summary_tables)
    
    


    if args.software == "macs2":
        print("Starting peak calling with macs2")
        
        bam_dir = args.bam
        
        genome_size = args.genome_size
        
        if  not skip_plot:
            sum_tbl = pd.read_csv(tbl)
        else:
            sum_tbl =False
        
        peakcall_macs2(peakCalling, bam_dir,control,percentage,summary_tables, fragments,sum_tbl,skip_plot,genome_size,is_percent)
    else:
        print("Starting peak calling with SEACR")
        
        if not skip_plot:
            sum_tbl = pd.read_csv(tbl)
        else: sum_tbl = False
        
        
        seacr_dir = pkg_resources.resource_filename('epimapper', 'seacr_tool')
        seacr_path= os.path.join(seacr_dir,"SEACR_1.3.sh")
        bedgraph = args.bedgraph
   
        #test jbw 16.06
        norm=args.seacr_norm
 
        peakcall_seacr(seacr_path,peakCalling, summary_tables,sum_tbl,bedgraph,control, fragments, percentage,skip_plot,norm)
        #end test   
    
if __name__=='__main__':
    
    """
    Exceutes the script if filename is called in terminal/ shell script. 
    
    
    Runs functions:
        
        * run()
    
    """
    
    
    
    args = set_parser(argparse.ArgumentParser('python peak_calling.py')).parse_args()
    
    run(args)
    
    
    
                
  
    
    
    
    
    
    
    
    
