#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import glob
import pathlib as pl
import re
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mlp
import subprocess
mlp.use("agg")

def set_parser(parser):
    
    """
    
      
    Creates and returns variables from terminal inputs
    
    Function input from shell or terminal:
    
        - Required input:
            * --sam (-s), str, full path to sam directory containing the sam files being analyzed
            
            * --chromosome_sizes (-cs), str, full path to file containing chromosome size information about the genome
            
            * --blacklist (-bl), str,  full path to bed file with blacklisted regions of the genome
            
            
            
        - Optinal input:
            * --out_dir (-o), str, The full pathway to desired output directory, where Epimapper directory will be made
        
    
        Function output:
            
            * args, object, containing the input data mentioned above
            
    """
    
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False, description="Input sam files for fragment length analysis")
    
    required_name=parser.add_argument_group("required arguments")
    
    optional_name = parser.add_argument_group("optional")
    
    required_name.add_argument("-o", "--out_dir", help = "Output directory", required = False, type = str)
    
    optional_name.add_argument("-s", "--sam", help = "Input directory with sam files", required = True, type = str)
    
    optional_name.add_argument("-cs", "--chromosome_sizes", help = "Path to file with chomosome sizes", required = True, type = str)
    
    optional_name.add_argument("-bl", "--blacklist", help= "Path to genome blacklist", required  = True, type = str)
    
    optional_name.add_argument("-sn", "--spike_in_norm", help = "If the samples are being spike-in normalized or not",required = False, type = bool, default = False)
    
    
    
    
    
    return parser


def make_bg(spike_in_norm, bed, chromosome_sizes, new_name, bedgraph):
      
      if not spike_in_norm:
          
          cmd_bg = "bedtools genomecov -bg -i " + bed +"/" + new_name + ".fragments_sorted.bed  -g "+ chromosome_sizes + " > " + bedgraph + "/" + new_name+ ".fragments.bedgraph"
          
          subprocess.run(cmd_bg,shell=True)
 
        
def check_command(out):
    if not out.returncode==0:
        print("Error in filtering. Exitig.")
        exit(1)
    
    
def do_filtering(LEN, files, minQualityScore, sam ,  sam_quality,
		bam, bed,chromosome_sizes,  align, genome_blacklist, genome, bedgraph, spike_in_norm ):
    
    """
    Uses samtools and bedtools to filter and convert sam files to quality and blacklist filtered bed files. 
    
    
    
    cmd1, cmd2 --> Takes blacklist bed file and chromosome sizes file to make a blacklist window bed file. 
    
    cmd3, cmd4 --> Filteres the sam files by removing the reads with a lower score than a set quilaity score and converts them to bam files.
    
    cmd5, cmd6 --> Filters and keeps the mapped read pairs, and sortes the bam files.
    
    cmd7 --> Converts bam files to bed files
    
    cmd8 --> Keeps the read pairs that are on the same chromosome and fragment length less than 1000bp.
    
    cmd9 --> Extracts the fragment related columns.
    
    cmd10 --> Maps the bed files to the Blacklist bed files and removes the overlaps.


    Function input:
        
        * LEN, str, string of a number defineing how many basepairs the  blacklist windows should be
        
        * files, list, list of stirng containing paths to the sam files being filtered and converted
        
        * minQualityScore, str, string of a number defining the lowest quality score for reads being filtered
        
        * sam, str, full path to the directory containing the sam files being filtered and converted
        
        * sam_quality, str, full path to the directory where quality score filtred sam files will be outputed
        
        * bam, str, full path to the directory where bam files will be outputed
        
        * bed, str, full path to the directory where bed files will be outputed
        
        * chromosome_sizes, str, full path to file containing chromosome size information about the genome
        
        * align, str, full path to where bed file containing blacklist window will be outputed
        
        * genome_blacklist, str, full path to bed file with blacklist areas of the genome
        
        * genome, str, genome name
        
    
    Function output:
        
        * tmp_files, list, list containing the path to quality score and blacklist filtered bed files
        
    
    
    
    """
    
    
    cmd1  = "bedtools makewindows -g "+ chromosome_sizes + " -w "+LEN+"  > "+ align+ "/" + genome +"."+ LEN+"b.windows.bed"
    out1 = subprocess.run(cmd1,shell=True)
    check_command(out1)
    
    cmd2 = "bedtools intersect -v -a "+align+ "/" + genome+ "."+ LEN+"b.windows.bed" + " -b "+ genome_blacklist+ " > "+ align + "/"+genome+ "."+ LEN+".b.windows.BlackListFiltered.bed"
    out2 = subprocess.run(cmd2,shell=True)
    check_command(out2)
        
        
    for file in files:
        
        new_name = re.split(r"[.]",pl.PurePath(file).name)[0]
       
        tmp_file_name = pl.PurePath(file).name
    
        cmd3 = "samtools view -Shq"+ minQualityScore+ " "+ sam + "/"+ tmp_file_name +"> "+ sam_quality + "/"+ new_name + "_QualityScore_" + minQualityScore+ ".sam"
        out3 =subprocess.run(cmd3,shell=True)
        check_command(out3)
        
        cmd4 = "samtools view -bS "+  sam_quality + "/"+ new_name + "_QualityScore_" + minQualityScore+ ".sam" +" >" + bam +"/"+new_name +".bam"
        out4 =subprocess.run(cmd4,shell=True)
        check_command(out4)
        
        cmd5 = "samtools view -b -F 0x04 "+ bam + "/" +new_name+".bam > " + bam +"/" + new_name+ ".mapped.bam"
        out5 = subprocess.run(cmd5,shell=True)
        check_command(out5)
        
        cmd6 = "samtools sort -n " + bam + "/" + new_name + ".mapped.bam -o " + bam +"/" + new_name+".mapped_sorted.bam "
        out6= subprocess.run(cmd6,shell=True)
        check_command(out6)
        
        cmd7 = "bedtools bamtobed -i " + bam +"/"+ new_name +".mapped_sorted.bam -bedpe > " + bed +"/" +new_name +".bed"
        out7=subprocess.run(cmd7,shell=True)
        check_command(out7)
        
        cmd8 = "awk '$1==$4 && $6-$2 < 1000 {print $0}' " + bed + "/" + new_name + ".bed" + "> "+ bed +"/" + new_name + ".clean.bed"
        out8=subprocess.run(cmd8,shell=True)
        check_command(out8)
        
        cmd9 = "cut -f 1,2,6 " + bed +"/" + new_name + ".clean.bed"+ " | sort -k1,1V -k2,2n -k3,3n  > " + bed +"/" + new_name + ".fragments.bed"
        out9=subprocess.run(cmd9,shell=True)
        check_command(out9)
        
   
        
        
        
        
        in_pd=pd.read_csv(bed +"/" + new_name + ".fragments.bed", sep='\t',header=None)

        in_pd.insert(3,3,1)
        
        in_pd.insert(4,'4',list(in_pd.iloc[:,2]-in_pd.iloc[:,1]),True)
    
        chrm = in_pd[in_pd[0].str.contains("M", na=False)]

        in_pd = in_pd.drop(list(chrm.index))

        new_pd = in_pd.append(chrm)

        new_pd[0] =  new_pd[0].apply(lambda x: 'chr' + str(x) if not str(x).startswith('chr') else str(x))
         
        new_pd.to_csv(bed +"/" + new_name + ".fragments_sorted.bed",sep='\t',index=False,header=None)
    
        cmd10= "bedtools map -g " +chromosome_sizes+" -a " +align + "/"+genome+ "."+ LEN+".b.windows.BlackListFiltered.bed -b "+ bed +"/" + new_name + ".fragments_sorted.bed  -c 4 -o sum > " \
            + bed+"/"+new_name+".fragments_sorted."+LEN+"b.windows.BlackListFiltered.bed"
    
        out10 =subprocess.run(cmd10,shell=True)
        check_command(out10)
   
        if not spike_in_norm:
            print("Creating bedgraph files")
            make_bg(spike_in_norm, bed, chromosome_sizes, new_name, bedgraph)
            
    cmd_rm_3 = "rm -r " +sam_quality
    subprocess.run(cmd_rm_3,shell=True)
     
    cmd_rm_8 = "rm " + os.path.join(bed,"*.clean.bed")
    subprocess.run(cmd_rm_8,shell=True)      
    tmp_file=os.path.join(bed,"*b.windows.BlackListFiltered.bed")
    
    tmp_files=glob.glob(tmp_file)
    
    print("Done with filering.")
    return tmp_files 
   
       
def calculation(filtered_files, count_cutoff):
    """
    
    
    Calculates log 2 transformed correlation matrixes between the different samples inputed. 
    
    1th Matrix is every count in the filtered bed files
    
    2nd Matrix is filtered, removing fragments with a count number lower than a input count_cutoff number
    
    
    
    
    Function input:
        
        * filtered_files, list, list of full paths to all filtred bed files being analyzed
        
        * count_cutoff, int, number determening the minimum number of counts of each fragment that will be filtered away before calculation the 2nd matrix
        
    
    Function output:
        
        * corr_df1, matrix, corrlation matrix between the different samples and their fragments
        
        * corr_df2 matrix, correlation matrix between the diffrent samples filtered on fragment counts
    
    
    
    
    """
    
    record_df=pd.DataFrame()
    for file in filtered_files:
        tmp_file_name = pl.PurePath(file).name
        sample = re.split(".fragments", tmp_file_name)[0]
        tmp_df=pd.read_csv(file,sep='\t',header=None)
        tmp_df.columns=['chrom','start','end',sample+"_counts"]
        tmp_df['ids']= tmp_df.chrom.astype(str)+[':']+tmp_df.start.astype(str)
        tmp_df['ids']=tmp_df.ids.apply(lambda x: x.strip() )
        tmp_df_clean = tmp_df[["ids",sample+"_counts"]]
        if record_df.shape[0] == 0:
            record_df = tmp_df_clean.copy()
        else:
            record_df=record_df.merge(tmp_df_clean,on=['ids'] ).copy()
    df=record_df
    df.index=df.ids.to_list()
    final_df=df.iloc[:,1:].copy()
    final_df=final_df.replace('.',np.nan)
    new_df=final_df.iloc[:,0:].dropna(axis=0, how='all').copy()
    new_df.fillna(0,inplace=True)
    new_df=np.log2(new_df.astype(float)+1)
    new_df = new_df.reindex(sorted(new_df.columns), axis=1)
    corr_df1 = new_df.corr()
    
    
    st_df=new_df.astype(float)<count_cutoff
    
    num_of_cols=st_df.shape[1] 
    
    filtered_new_df= new_df[~(st_df.sum(axis=1)==num_of_cols)].copy()

    st_df2= np.log2(filtered_new_df.astype(int)+1)
    
    st_df2=st_df2.reindex(sorted(st_df2.columns),axis=1)
    
    corr_df2=st_df2.corr()
    
    if corr_df1.shape[0] >0 and corr_df2.shape[0] >0:
        print("Done with calculating correlation matrixes")
    else:
        print("Error in calculating correlation matrixes")
        exit(1)
    
    
    
    
    return corr_df1, corr_df2
       
    
    
    
    
    
    

def plot_corr(corr_df1, corr_df2, summary_tables):
    
    """
    Plots heatmaps from correlation matrixes and saves them to directory summary_tables.
    
    
    Function input:
        
        * corr_df1, matrix, correlation matrix based in the  filtred bed files
        
        * corr_df2, matrix, count filtred correlation matrix based on the filtred bed files
        
        * summary_tables, str, full path to where the plots will be saved
        

    Function output:
        
        * Creates "corrcoef_heatmap4logCount_all.jpg" from corr_df1 and saves it to summary_tables
        
        * Creates "corrcoef_heatmap4logCount_filtered_gt8.jpg" from corr_df2 and saves it to summary_tables
    
    
    
    """
    sns.set(font_scale=1.4)
    plt.tight_layout()
    fig=plt.figure(figsize=(22,15))
    ch=sns.heatmap(corr_df1, annot=True, vmin =-1, vmax=1)
    ch.set_xticklabels(ch.get_xticklabels(),rotation=30, fontsize = 12)
    ch.set_yticklabels(ch.get_yticklabels(),rotation=45, fontsize= 12)

    out_fig1=os.path.join(summary_tables, "corrcoef_heatmap4logCount_all.jpg")
    fig.savefig(out_fig1)




    fig=plt.figure(figsize=(22,15))
    ch=sns.heatmap(corr_df2,annot=True, vmin=-1,vmax=1)
    ch.set_xticklabels(ch.get_xticklabels(),rotation=30)
    ch.set_yticklabels(ch.get_yticklabels(), rotation=45)
    out_fig2=os.path.join(summary_tables, "corrcoef_heatmap4logCount_filtered_gt8.jpg")

    fig.savefig(out_fig2)


    
    



def run(args):
    
    """
    
    Creates new directories, if not present, in the output directory:
        
        - Epimapper
        
        - Epimapper/alignment
        
        - Epimapper/alignment/sam_quality
        
        - Epimapper/alignment/bam
        
        - Epimapper/alignment/bed
        
        - Epimapper/summary_tables
        
        
    
    Function input:
        
        * args, object, containing the input from shell script or command line
        
        
    Runs functions created above with input data from args:
        
        * do_filtering()
        
        * calculation()
        
        * plot_corr()
        
    
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

    
    minQualityScore="2"
    
    path=os.path.join(projPath,"Epimapper")
    if not os.path.exists(path):
        os.mkdir(path)
        
    summary_tables = os.path.join(path, "summary_tables")
    if not os.path.exists(summary_tables):
        os.mkdir(summary_tables)
    
    align = os.path.join(path, "alignment")
    if not os.path.exists(align):
        os.mkdir(align)
        
    sam_quality = os.path.join(align,"sam_quality")
    if not os.path.exists(sam_quality):
        os.mkdir(sam_quality)
    
    bam = os.path.join(align,"bam")
    if not os.path.exists(bam):
        os.mkdir(bam)
        
    bed = os.path.join(align,"bed")
    if not os.path.exists(bed):
        os.mkdir(bed)
    
    spike_in_norm = args.spike_in_norm
    
    bedgraph = os.path.join(align,"bedgraph")
    if not spike_in_norm:
        if not os.path.exists(bedgraph):
            os.mkdir(bedgraph)
   
    
    sam = args.sam
    if  os.path.isdir(sam):
    
        files = glob.glob(os.path.join(sam,"*.sam"))
        if len(files) == 0:
            print("Chosen sam folder: " + sam +" is empty or does not contain any sam files \n Please check your directory or select another one")
            exit(1)
    else:
        print("Chosen sam folder: " + sam +" is not a directory or does not exist. \n Please select another directory.")
        exit(1)
    
    
    chromosome_sizes = args.chromosome_sizes
    if  os.path.exists(chromosome_sizes) and os.path.isfile(chromosome_sizes):
        if not os.path.getsize(chromosome_sizes)>0:
            print("Chosen chromosome sizes file: "+chromosome_sizes +" is empty.\n Please check your file or chose another one.")
            exit(1)
    else:
        print("Chosen chromosome sizes file: "+chromosome_sizes +" is not a file or does not exist. \n Please check your file or chose another one.")
        
    
    genome_blacklist = args.blacklist
    if os.path.exists(genome_blacklist) and os.path.isfile(genome_blacklist):
        if not os.path.getsize(genome_blacklist)>0:
            print("Chosen genome blacklist file: "+ genome_blacklist+" is empty.\n Please check your file or chose another one.")
            exit(1)
    else:
        print("Chosen genome blacklist file: "+ genome_blacklist+" is not a file or does not exist. \n Please check your file or chose another one.")
        exit(1)
    
    genome = pl.PurePath(chromosome_sizes).name
    
    
    LEN = "500"
    
    
    count_cutoff=8
    
    
    filtered_files2 = do_filtering(LEN, files, minQualityScore, sam, sam_quality, bam, bed, chromosome_sizes, align, genome_blacklist, genome,bedgraph, args.spike_in_norm)
    
    
    corr_df1,corr_df2 = calculation(filtered_files2, count_cutoff)
    
    cmd_rm = "rm " + os.path.join(align, "*.bed")
    subprocess.run(cmd_rm,shell=True)
    plot_corr(corr_df1, corr_df2, summary_tables)
    print("Done with plotting. Plots avalible at: "+summary_tables)
    
    
if __name__=='__main__':
    
    """
    
    Exceutes the script if filename is called in terminal/ shell script. 
    
    
    Runs functions:
        
        * run()
    
    """
    
    
    
        
    args = set_parser(argparse.ArgumentParser('python filtering.py')).parse_args()
    
    run(args)
 
   
            

    
    

        
        
        
        