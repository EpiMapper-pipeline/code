import os
import glob
import pathlib as pl
import re
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import shutil
import subprocess
import matplotlib as mlp
mlp.use("agg")




def set_parser(parser):
    """
    
    Creates and returns variables from terminal inputs
    
    Function input from shell or terminal:
    
        - Required input:
            
            * --sam_spike_in(-ss), str, full path to sam directory containing the sam files bowtie2 alignment to spike in genome
            
            * --chromosome_sizes (-cs), str, full path to file containing chromosome size information about the genome
            
            * --bed (-b), str, full path to directory with bed files being analyzes
            
            
        - Optinal input:
            
            * --out_dir(-o), str, The full pathway to desired output directory, where Epimapper directory will be made
            
            *--fragment_table, (-tbl), str, full path to table with information about number of mapped fragments, column names = ["Sample",	"Replication", "SequencingDepth", "MappedFragments", "AlignmentRate", "MappedFragments_SpikeIn",	"AlignmentRate_SpikeIn"]
                             Default: will be collected from summary_tables, made in python script "bowtie2_alignment.py"
        
    
        Function output:
            
            * args, object, containing the input data mentioned above
            
            
    
    
    
    """
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False, description="Input sam files for fragment length analysis")
    
    required_name=parser.add_argument_group("required arguments")
    
    optional_name = parser.add_argument_group("optional")
    
    required_name.add_argument("-ss", "--sam_spike_in", required = True, type = str)
    
    required_name.add_argument("-cs", "--chromosome_sizes", required = True, type = str)
    
    required_name.add_argument("-b","--bed", required = True, type = str)
    
    optional_name.add_argument("-o", "--out_dir", required = False, type= str)
    
    optional_name.add_argument("-tbl", "--fragment_table", required = False, type= str)
    
   
    
    
    return parser
    
    
    
def calibration(tmp_files, depth_files,bed_files,bedgraph, chromosome_sizes):
    
    """
    Calculates a scaling factor for each sample by:
    
        Scaling factor = C / (fragments mapped to spike-in genome)
    
    Where C is a constant, here 10000, to avoid small fractions in the normalized data.
    
    
    
    Uses bedtools genomecov to normalize and scale with the scaling factor and output bedgraph files.
    
    
    
    Function input:
        
        * tmp_files, list, containing the full path to sam files from spike in alignment with bowtie2
        
        * depth_files, str, full path to where the sequencing depth for each file will be stored
        
        * bed_files, str, full path to directory where bed files being analyzed are
        
        * bedgraph, str, full path to directory where bedgraph files will be saved
        
        * chromosome_sizes, str, full path to file file containing chromosome size information about the genome

    
    Funtion output:
        
        * Function outputs spike in calibrated and normalized bedgraph files for each sample to bedgraph directory
        
        

    """
    
    
    for file in tmp_files:
        new_name = re.split("_spike_in",pl.PurePath(file).name)[0]
       
        cmd1 = 'seqDepthDouble=$(samtools view -F 0x04 '+ file + ' | wc -l) \n seqDepth=$((seqDepthDouble/2)) \
            \n echo ${seqDepth} \n echo $seqDepth >' + depth_files+ '/'+new_name +'.spikeIn.seqDepth'
            
        out1 =subprocess.run(cmd1,shell=True)
        if not out1.returncode==0:
            print("Error in spike in calibration. Exiting.")
            exit(1)
        seqdepth = pd.read_csv(depth_files+ '/'+new_name +'.spikeIn.seqDepth', header=None)
        
        scale_factor = round(10000/seqdepth.iloc[0, 0],3)
        cmd2 = 'bedtools genomecov -bg -scale '+str(scale_factor)+' -i '+bed_files+'/'+new_name+ '.fragments_sorted.bed -g '+ chromosome_sizes + ' > ' + bedgraph + '/' + new_name+ '.fragments.normalized.bedgraph'
        
        out2 = subprocess.run(cmd2,shell=True)
        if not out2.returncode==0:
            print("Error in spike in calibration. Exiting.")
            exit(1)
        





def summary(depth_files, summary_tables, sum_tbl):
    
    """
    Creates summary table from spike in calibration and saves them to directory summary_tables.
    
    
    Calculates normalized fragments by muliplying the Mapped fragments with the scaling factor.
    
    
    Function input:
        
        * depth_files, str, full path to directory containing the spike in sequencing depth files
        
        * summary_tables, str, full path to directory where summary table will be saved
        
        * sum_tbl, str, full path to table with information about number of mapped fragments, column names = ["Sample",	"Replication", "SequencingDepth", "MappedFragments", "AlignmentRate", "MappedFragments_SpikeIn",	"AlignmentRate_SpikeIn"]
        
    
    Function output:
        
        * final_summary, pandas DataFrame, summary dataframe containting columns : [Sample,	Replication	, SequencingDepth, MappedFragments, AlignmentRate, MappedFragments_SpikeIn, AlignmentRate_SpikeIn, Scale_factor, NormalizedFragments]
    
    """
    
    
    scaleFactor = pd.DataFrame(columns =["Sample", "Replication" ,"Scale_factor"] )
    
    multiplier = 10000
    
    files = os.path.join(depth_files,"*.spikeIn.seqDepth")
    tmp_files = glob.glob(files)
    tmp_files.sort()
    for file in tmp_files:
        file_name = pl.PurePath(file).name
        name = file_name.split("_")[0]
        rep = "".join(re.findall("\Brep\d+", file_name))
        with open(file, "r") as t:
            seqDepth = t.readline()
            #test jbw 14.06
            if int(seqDepth)==0:
              scaleFactor.loc[len(scaleFactor)] = [name,rep,0 ]  
            else:
              scaleFactor.loc[len(scaleFactor)] = [name,rep,round(multiplier/int(seqDepth),2) ]
            #end test
            
        t.close()

    
    main_tbl = pd.read_csv(sum_tbl, index_col=0)
    
    
    final_summary = pd.merge(main_tbl, scaleFactor, on = ["Sample","Replication"])
    
    
    final_summary["NormalizedFragments"] = final_summary["MappedFragments"]* final_summary["Scale_factor"]
    
    
    if final_summary.shape[0]>0:
        
        print("Done with creating spike-in normalization summary. Summary avalible at: "+ os.path.join(summary_tables,"spike_in_calibration_summary.csv"))
    
        final_summary.to_csv(os.path.join(summary_tables,"spike_in_calibration_summary.csv"))
    else:
        print("Error in creating spike-in normalization summary. Exiting.")
        exit(1)
        
    return final_summary






def plot(final_summary, summary_tables):
    """
    Plots boxplots for sacling factor and normalized fragment counts from input summary table.
    
    Function input:
        
        * final_summary, pandas DataFrame, DataFrame with summmary information from spike in calibration and columns:  [Sample,	Replication	, SequencingDepth, MappedFragments, AlignmentRate, MappedFragments_SpikeIn, AlignmentRate_SpikeIn, Scale_factor, NormalizedFragments]
    
        * summary_tables, str, full path to where the plot will be saved
    
        
    Function outputs:
        
        * spike_in_calibration.png, png file, one png file containing 2 plots: spike in scaling factor boxplot and normalized fragment count boxplot
                                  Saved to directory summary_tables.
        
    
    """
    fig, axes= plt.subplots(nrows=1, ncols=2)
    sns.set_style("whitegrid")
    sns.set_palette("husl")
    fig, axes= plt.subplots(nrows=1, ncols=2)
    plt1 = sns.boxplot( x = final_summary.Sample, y = final_summary.Scale_factor, hue = final_summary.Sample, dodge=False, ax=axes[0])
    plt1.set_ylabel("Spike in Scaling Factor")
    plt.setp(plt1.xaxis.get_majorticklabels(), rotation=90)
    plt1.legend([],[], frameon=False)
    
    
    
    plt2 = sns.boxplot(x=final_summary.Sample, y = final_summary.NormalizedFragments, hue = final_summary.Sample, dodge=False, ax= axes[1])
    plt2.set_ylabel("Normalization Fragment Count")
    plt.setp(plt2.xaxis.get_majorticklabels(), rotation=90)
    
    fig.tight_layout(pad=3.0)  

    fig.savefig(os.path.join(summary_tables, "spike_in_calibration.png"))          


def run(args):
    
    """
  
    
    Creates new directories, if not present, in the output directory:
        
        - Epimapper
        
        - Epimapper/alignment
        
        - Epimapper/alignment/bedgraph
        
        - args.spike_in_sam/seqDepth
         
        - Epimapper/summary_tables
        
    
        
        
        
        
    Function input:
            
        * args, class, containing the input from shell script or command line
            
            
            
    
    Runs functions created above with input data from args:
        
        * calibration()
        
        * summary()
        
        * plot()
        
    
    
    
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
    
    summary_tables = os.path.join(path,"summary_tables")
    if not os.path.exists(summary_tables):
        os.mkdir(summary_tables)
    
    align = os.path.join(path, "alignment")
    if not os.path.exists(align):
        os.mkdir(align)
    
    bedgraph = os.path.join(align, "bedgraph")
    if not os.path.exists(bedgraph):
        os.mkdir(bedgraph)
    
    
    depth_files = os.path.join(args.sam_spike_in, "seqDepth")
    if  os.path.exists(depth_files):
            shutil.rmtree(depth_files)
            os.mkdir(depth_files)
    else:
            os.mkdir(depth_files)
                
    chromosome_sizes=args.chromosome_sizes         
    if  os.path.exists(chromosome_sizes) and os.path.isfile(chromosome_sizes):
        if not os.path.getsize(chromosome_sizes)>0:
            print("Chosen chromosome sizes file: "+chromosome_sizes +" is empty.\n Please check your file or chose another one.")
            exit(1)
    else:
        print("Chosen chromosome sizes file: "+chromosome_sizes +" is not a file or does not exist. \n Please check your file or chose another one.")
   
    sam = args.sam_spike_in
    if  os.path.isdir(sam):
   
       tmp_files = glob.glob(os.path.join(sam,"*.sam"))
       if len(tmp_files) < 1:
           print("Chosen sam spike in folder: " + sam +" is empty or does not contain any sam files. \n Please check your directory or select another one")
           exit(1)
    else:
       print("Chosen sam spike in folder: " + sam +" is not a directory or does not exist. \n Please select another directory.")
       exit(1)
   
    
    bed_files = args.bed
    if os.path.exists(bed_files):
        if len(glob.glob(bed_files+"/"+"*.bed")) < 1:
            print("Chosen bed directory: "+bed_files+" is empty or does not contain any bed files. \n Please check your directory or select another one")
    
    
    if args.fragment_table:
        sum_tbl = args.fragment_table
        if os.path.exists(sum_tbl) and os.path.isfile(sum_tbl):
            if os.path.getsize(sum_tbl) ==0:
                print("Chosen summary table file: "+ sum_tbl+" is empty.\n Please check your file or provide another one.")
                exit(1)
            df = pd.read_csv(sum_tbl,nrows=1)
            column_names = df.columns.tolist()
            if not ["Sample",	"Replication", "SequencingDepth", "MappedFragments", "AlignmentRate", "MappedFragments_SpikeIn",	"AlignmentRate_SpikeIn"] in column_names:
                print('Chosen summary table does not contain the following columns:\n ["Sample",	"Replication", "SequencingDepth", "MappedFragments", "AlignmentRate", "MappedFragments_SpikeIn",	"AlignmentRate_SpikeIn"]. Please update your table´s column names or chose another table')
                exit(1)
            
    else:
        sum_tbl = os.path.join(summary_tables,"bowtie2_alignment_ref_and_spike_in.csv")
        
        if not os.path.exists(sum_tbl):
            while True:
                x = input("No summary table provided and non avalible in :"+ summary_tables+". \n: Do you wish to continue spike in normalization of files without plot generation? y/n")
                if "y" in x or "yes" in x:
                    calibration(tmp_files, depth_files, bed_files, bedgraph, chromosome_sizes)
                    print("Done with spike in calibration.\n Spike in normalized files may be found in: "+ bedgraph+"\n Exitig.")
                    exit(0)
                    
                elif "n" in x or "no" in x:
                    print("Exiting.")
                    exit(0)
                else:
                    print("Invalid input. Please type 'y' or 'n'.")
        
            
                
    
    
    calibration(tmp_files, depth_files, bed_files, bedgraph, chromosome_sizes)
    
    final_summary = summary(depth_files, summary_tables, sum_tbl)
    
    plot(final_summary,summary_tables)
    print("Done with plotting. Plots avalible at: "+summary_tables)
    shutil.rmtree(depth_files)




if __name__=='__main__':
    
    
    """
    Exceutes the script if filename is called in terminal/ shell script. 
    
    
    Runs functions:
        
        * run()
    
    """
    
    
    
    
    args = set_parser(argparse.ArgumentParser('python spike_in_calibration.py')).parse_args()

    run(args)     

    
        
        
