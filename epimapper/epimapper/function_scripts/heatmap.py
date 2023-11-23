#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import os
import glob
import pathlib as pl
import pandas as pd
import subprocess










def set_parser(parser):
    """
    
      
    Creates and returns variables from terminal inputs
    
    Function input from shell or terminal:
    
    -Required input:
        
        * --bam_dir (-b), full pathway to directory with bam files being used to make bigwig files used to create matrixes
        
        * --ref (-r), full pathway to a refFlat file 
        
        * --peaks (-p), full pathway to directory with bed files containing peaks from SEACR peak calling
        
        * --blacklist (-bl), full pathway to bed file containing blacklisted parts of the genome
        
        
	-Optional input:
        
        * --out_dir(-o): The full pathway to desired output directory, where Epimapper directory will be made

        * --cores (-c), number of cores being used in matrix calculation, default = 8
    
    Function output:
        
        * args, object, containing the input data mentioned above
    


    """
    
    
    
    
    
    
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False, description="Input files for alignment")
    
    required_name=parser.add_argument_group("required arguments")
    
    optional_name = parser.add_argument_group("optinal arguments")
    
    required_name.add_argument("-b","--bam_dir", help="Input directory name", required=True, type=str)
    
    optional_name.add_argument("-out_dir", "--out_dir", help="Output directory", required=False, type=str)
    
    required_name.add_argument("-r","--ref", help="Referance BED file", required=True, type=str)
    
    required_name.add_argument("-p","--peaks", help="Peaks from SEACR", required=True, type=str)
    
    optional_name.add_argument("-c", "--cores", help="number of cores", required=False, type=str, default ="8")
    
    required_name.add_argument("-bl", "--blacklist", help= "Path to genome blacklist", required  = True, type = str)
    
    
    return parser



def create_bed(ref):

    refflat = pd.read_csv(ref, header=None, sep="\t")
    
    bed = refflat[[2,4,5,1,3]].copy()
    
    bed.to_csv(ref.replace(".txt", ".bed"), header=False, sep="\t", index=False) 
    
    cmd= "sort -k1,1V -k2,2n -k3,3 " + ref.replace(".txt", ".bed") +" > "+ ref.replace(".txt", "_sorted.bed")
    
    subprocess.run(cmd,shell=True)
    
    sorted_bed = pd.read_csv(ref.replace(".txt", "_sorted.bed"), header=None, sep="\t", names =["chr","start","end","gene","strand"])
    
    sorted_bed = sorted_bed[sorted_bed.chr.astype(str).str.contains('_') == False]
    
    chrm = sorted_bed[sorted_bed["chr"].astype(str).str.contains("M", na=False)]
    
    sorted_bed = sorted_bed.drop(list(chrm.index))

    sorted_bed = sorted_bed.append(chrm)
    
    rm_cmd = "rm " + ref.replace(".txt", "_sorted.bed")
    rm_cmd2= "rm " + ref.replace(".txt", ".bed")
    subprocess.run(rm_cmd,shell=True)
    subprocess.run(rm_cmd2,shell=True)
    
    sorted_bed.to_csv(ref.replace(".txt", "_cleaned.bed"),header=False, sep="\t", index=False)
    
    return ref.replace(".txt", "_cleaned.bed")
    
def coverage(tmp_files, bam_dir,bigwig, sample_names):   
    """
    Uses samtools and bamCoverage to sort bam files and create bigwig files
    
    Function input:
        
        * tmp_files, list, list of paths to bam files being analyzed
        
        * bam_dir, str, full path to directory with bam files being analyzed
        
        * bigwig, str, full path to directory where bigwig files will be outputed
    
    """
    for file in tmp_files:
        
        
        file_name  = pl.Path(file).name
        
        name = file_name.split("_")[0]
        
        if name in sample_names:
            
            sample =  file_name.split(".")[0]
            cmd = " \
                samtools sort " +file +" -o " + bam_dir +"/" + sample + ".chr_sorted.bam " + "\n" \
                    "samtools index " + bam_dir +"/" + sample + ".chr_sorted.bam " + "\n" \
                        "bamCoverage --bam " + bam_dir +"/" + sample + ".chr_sorted.bam" + " -o " + bigwig + "/" + sample + ".bw"
           
            out_code = subprocess.run(cmd,shell=True)
            if out_code.returncode  ==0:
                print("Done with - Creating bigwig files")
            else: 
                print("Error in - Creating bigwig files\n Please check your bam files")
                exit(1)
        
     



        
def compute_one_matrix(bigwig,ref,cores,bl): 
    """
    Uses deeptools to create a matrix over transcription units
    
    Function input:
        
        * bigwig, str, full pathway to directory containing bigwig files used to calculate the matrix
            The matrix will be outputed in bigwig folder
        
        * ref, str, full pathway to sorted reference bed file of genome 
        
        * cores, str, number of cores being used in the calculation
        
        * bl, str, full pathway to bed file containing blacklisted parts of the genome
        
        

    """
    
    bw_files = glob.glob(os.path.join(bigwig, "*.bw"))
    bw_files.sort()
    txt = " ".join(str(e) for e in bw_files)
    
        
    cmd2 =" computeMatrix scale-regions -S " + txt + "  -R " + ref + "  --beforeRegionStartLength 3000 --regionBodyLength 5000  --afterRegionStartLength 3000 --skipZeros -p " + cores \
        + " -bl " + bl + " -o " +bigwig+"/matrix_gene.mat.gz"
    
    
    out_code2 = subprocess.run(cmd2, shell=True)
    
    if out_code2.returncode  ==0:
        print("Done wtih - Computing region matrix")
    else:
        print("Error in - Computing region matrix")
        exit(1)
 
    
def plot_one_heatmap(bigwig, summary_tables):
    """
    
    Uses deeptools to create heatmap over transcriprion units
    
    Function input:
        
        * bigwig, str, full pathway to directory containing the matrix being used to create the heatmaps
        
        * summary_tables, str, full pathway to where the heatmap will be outputed
        
        
    
    """
    cmd3 = "plotHeatmap -m " + bigwig +"/matrix_gene.mat.gz -out "+ summary_tables+ "/matrix_heatmap.png --sortUsing sum"
    
    out_code3 =subprocess.run(cmd3, shell=True)
    
    if out_code3.returncode  ==0:
        print("Done with - Plotting region heatmap")
    else:
        print("Error in - Plotting region heatmap")
        exit(1)
    
    

def summitRegion(peaks_files):
    """
    Generates new bed files containing the midpoint information in column 6 of SEACR files to find the midpoint of signal block to align signals in heatmaps


    Function input:
        
        * peaks_files, list, list of paths to bed files from SEACR peak calling being used to generate heatmaps.
        
        
    Function output:
        
        * sample_names, list, list of sample names from the input files. 

    """
    
    sample_names = []
    for file in peaks_files:
        file_name  = pl.Path(file).name
        sample = file_name.split("_")[0]+"_"+file_name.split("_")[1]
        sample_names.append(sample)     
        
        txt2 = """awk  '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' """
        cmd4 = txt2 + file + " > " + file.replace(".bed", ".summitRegion.bed")
        
        out_code4 =subprocess.run(cmd4, shell=True)
        
        if out_code4.returncode ==0:
            print("Done with - Finding summit regions for "+sample)
        else: 
            print("Error in - Finding summit regions for "+ sample)
            exit(1)
       
                       
    sample_names = list(dict.fromkeys(sample_names))           
    sample_names.sort()

    return sample_names




def compute_many_matrixes(sample_names,bigwig,cores, peaks):
    
    
    """
    
    Uses deeptools to create matrixes on cut&tag peaks, based on bigwig and peak files.
    
    Function input:
        
        * sample_names, str, list of sample names from the input files
        
        * bigwig, str, full pathway to directory containing bigwig files used to calculate the matrix
         
        * cores, str, number of cores being used in the calculation
        
        * peaks, str, full pathway to directory containing peak bed files from SEACR analyzis
            The matrixes will be outputed in this directory 
        
        

    """
    
    for sample in sample_names:
        peak_file= glob.glob(os.path.join(peaks,sample+"*"))[0]
        cmd6 = "computeMatrix reference-point -S " +bigwig + "/"+ sample + ".bw -R " +peak_file +" --skipZeros -a 3000 -b 3000 --referencePoint center " \
            + "-p " + cores + " -o " + peaks + "/" + sample+"_SEACR.mat.gz"

        out_code6 = subprocess.run(cmd6, shell=True)
        
        if out_code6.returncode  == 0:
            print("Done with - Computing peak matrix for "+ sample)
        else:
            print("Error in - Computing peak matrix for " + sample)
            exit(1)
            
        

def plot_many_heatmaps(sample_names,summary_tables,peaks):
    """
    Uses deeptools to create heatmapes of cut&tag peaks
    
    
    
    Function input:
        
        * sample_names, str, list of sample names from the input files
        
        *  summary_tables, str, full pathway to where the heatmap will be outputed
        
        * peaks, str, full pathway to directory containing peak bed files from SEACR analyzis
        
    
    
    """
    for sample in sample_names:
        cmd7 = 'plotHeatmap --matrixFile ' + peaks  + "/"+ sample+'_SEACR.mat.gz --outFileName ' + summary_tables + '/' + sample +'SEACR_heatmap.png --sortUsing sum --startLabel "Peak start" --endLabel "Peak end" ' \
        + '--xAxisLabel "" --regionsLabel "Peaks" --samplesLabel ' + sample
            
        
        out_code7 =subprocess.run(cmd7, shell=True)
        if out_code7.returncode == 0:
            print("Done with - Plotting peaks heatmap for "+sample)
        else:
            print("Error in - Plotting peaks heatmap for "+sample)
            exit(1)

def check_input(args):
    
    bam_dir = args.bam_dir
    if os.path.exists(bam_dir):
        b_files = glob.glob(os.path.join(bam_dir, "*.mapped_sorted.bam"))
        if len(b_files) <1:
            print("Chosen bam directory: "+bam_dir+" is empty or does not contain any mapped_sorted.bam files. \n Please check your directory or select another one.")
            exit(1)
    else:
        print("Chosen bam directory: "+bam_dir+" does not exist. \nPlease select check your path or select another one")
        exit(1)
    
    peaks = args.peaks
    if os.path.exists(peaks):
         p_files =glob.glob(os.path.join(peaks, "*peaks*.bed"))
         if len(p_files) <1:
             print("Chosen  peak directory: "+peaks+" is empty or does not contain any peak bed files. \nPlease check your directory or select another one.")
             exit(1)

                 
    else:
         print("Chosen peak directory: "+peaks+" does not exist. \nPlease select check your path or select another one")
         exit(1)

    bl = args.blacklist
    if os.path.exists(bl) and os.path.isfile(bl):
        if not os.path.getsize(bl)>0:
           print("Chosen genome blacklist file: "+ bl+" is empty.\n Please check your file or chose another one.")
           exit(1)
    else:
        print("Chosen genome blacklist file: "+ bl+" is not a file or does not exist. \n Please check your file or chose another one.")
        exit(1)
    
    ref = args.ref
    if  os.path.exists(ref) and os.path.isfile(ref):
        if not os.path.getsize(ref)>0:
            print("Chosen reference refFlat genome file: "+ref +" is empty.\n Please check your file or chose another one.")
            exit(1)
    else:
        print("Chosen reference file: "+ ref+ " is not a file or does not exist. \n Please check your file or chose another one.")
        exit(1)
        
        
    if args.cores is not None:
        cores = args.cores
        try:
            float(cores)
        except ValueError:
            print("Number of cores is not numeric. /n Please provide a number of how many cores to be used in the calcualtion in the -c parameter, default = 8.")
            exit(1)
        
        
    return bam_dir, peaks,bl,ref, cores
    
    


        
def run(args):
    """
    Creates new directories, if not present, in the output directory:
        
        
        - Epimapper
        
        - Epimapper/alignment
        
        - Epimapper/alignment/bigwig
        
        - Epimapper/summary_tables
        
        
        
        
        Function input:
            
            * args, class, containing the input from shell script or command line
            
            
            
    
    Runs functions created above with input data from args:
        
        * coverage()
        
        * compute_one_matrix()
        
        * summitRegion()
        
        * plot_one_heatmap()
        
        * compute_many_matrixes()
        
        * plot_many_heatmaps()
    
    """
    
    
    
    
    
    if args.out_dir:
        
        projPath = args.out_dir
    else:
        projPath = os.getcwd()
    
    path=os.path.join(projPath,"Epimapper")

    if not os.path.exists(path):
        os.mkdir(path)
        
        
    summary_tables=os.path.join(path,"summary_tables")

    if not os.path.exists(summary_tables):
        os.mkdir(summary_tables)    
    
    align = os.path.join(path, "alignment")

    if not os.path.exists(align):
        os.mkdir(align)
    
    bigwig = os.path.join(align, "bigwig")

    if not os.path.exists(bigwig):
        os.mkdir(bigwig)
        
    bam_dir, peaks,bl,ref_txt, cores = check_input(args)
    
    ref = create_bed(ref_txt)
    
    files = os.path.join(bam_dir,"*.mapped_sorted.bam")
    
    tmp_files = glob.glob(files)
    
    


        
    peaks_files = os.path.join(peaks, "*peaks*.bed")
    p_files =peaks_files 
        
    peaks_files = glob.glob(p_files)
    
    
    a =[os.path.basename(file) for file in peaks_files]
    
    sample_names =[name.split("_")[0] for name in a]
    
    coverage(tmp_files, bam_dir, bigwig, sample_names)
    
    compute_one_matrix(bigwig, ref, cores, bl)
    
    sample_names = summitRegion(peaks_files)
    
    plot_one_heatmap(bigwig,summary_tables)
    
    compute_many_matrixes(sample_names, bigwig, cores, peaks)
    
    plot_many_heatmaps(sample_names, summary_tables, peaks)
        
if __name__=='__main__':
    """
    
    
    Exceutes the script if filename is called in terminal/ shell script. 
    
    
    Runs functions:
        
        * run()
    
    """
    
    
    
    args = set_parser(argparse.ArgumentParser('python heatmap.py')).parse_args()
        
    run(args)
    
    
    
    
    
    
    
    
    