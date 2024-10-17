import pandas as pd
import glob
import os
import numpy as np
import subprocess 
from scipy.stats import zscore
import pathlib as pl


def bin_count_macs2(peak_file,chromosome_sizes, blacklist_removed_bin):
    
    rm_y = "grep -v ^chrY " +peak_file + " > " +peak_file.replace(".bed", "_cleaned.bed")
    subprocess.run(rm_y,shell=True)
    
    rm_y = "grep -v ^chrM " + peak_file.replace(".bed", "_cleaned.bed") + " > " +peak_file.replace(".bed", "_2cleaned.bed")
    subprocess.run(rm_y,shell=True)
    
    #jbw 2024
    bed_sorted_file=peak_file.replace(".bed","_sorted.bed")
    sort_cmd = "sort -k1,1V -k2,2n -k3,3n " + peak_file.replace(".bed", "_2cleaned.bed") + " > " + bed_sorted_file
    subprocess.run(sort_cmd,shell=True)
    if (pl.Path(bed_sorted_file).stat().st_size<1):
       print("Empty in : ", bed_sorted_file, " is skipped !" )  
       fold_file = peak_file.replace(".bed", "_FolderEnrich.bed") 
       out_fold_file=None
    else:
       bin_cmd = "bedtools intersect -g " + chromosome_sizes + " -a " + blacklist_removed_bin + " -b " +peak_file.replace(".bed", "_sorted.bed")+ " -wao -sorted > "+ peak_file.replace(".bed", "_blacklistremoved_100bp.bed")  
       subprocess.run(bin_cmd, shell=True)
       fold_file =folderEnrichment(peak_file.replace(".bed", "_blacklistremoved_100bp.bed"))
       sort_fold = "sort -k1,1V -k2,2n -k3,3n " + fold_file + " > " + fold_file.replace(".bed", "_sorted.bed")
       subprocess.run(sort_fold,shell=True)
       out_fold_file= fold_file.replace(".bed", "_sorted.bed")

    cmd_rm = "rm -f " + peak_file.replace(".bed", "_cleaned.bed") + " " + peak_file.replace(".bed", "_blacklistremoved_100bp.bed") + " " + fold_file + " " + peak_file.replace(".bed", "_2cleaned.bed") +" " + peak_file.replace(".bed","_sorted.bed")  
    subprocess.run(cmd_rm,shell=True)

    #jbw 2024
    return out_fold_file
        
def folderEnrichment(macsFile):   
    out_file = macsFile.replace(".bed", "_FolderEnrich.bed")
    with open(out_file, "w") as out:
       num_of_bsite = 0
       uniq_data = {}
       all_data = {}
       #jbw 2024
       print(macsFile)
       print(out_file)
       with open(macsFile, "r") as data_file:
           for line in data_file:
               line = line.strip()
               lines = line.split("\t")
               #print(lines)
               tmp_line = "\t".join([lines[0], lines[1], lines[2], lines[6], lines[9]])
               all_data[num_of_bsite] = tmp_line
               num_of_bsite += 1
    
               identifier = ":".join([lines[0], lines[1], lines[2]])
               uniq_data.setdefault(identifier, []).append(tmp_line)
    
       uniq_ids = len(uniq_data)
      
    
       # Export unique id data with Maximum folderEnrichment used in export data
       num_of_exported = 0
       for identifier, tmp_data in uniq_data.items():
           if len(tmp_data) > 1:
               tmp_line = []
               max_data = 0
               tmp_id = ""
               for item in tmp_data:
                   tmp_line = item.split("\t")
                   if max_data < float(tmp_line[4]):
                       max_data = float(tmp_line[4])
                   tmp_id = ":".join([tmp_id, tmp_line[3]])
    
               tmp_out = "\t".join([tmp_line[0], tmp_line[1], tmp_line[2], tmp_id, str(max_data)])
               out.write(f"{tmp_out}\n")
           else:
               out.write(f"{tmp_data[0]}\n")
           num_of_exported += 1
    return out_file






def str2num(x, str1):
    if x==str1:
        x=0
    else:
        x=float(x)
    return x






def quantile_normalize(df):
    """
    input: dataframe with numerical columns
    output: dataframe with quantile normalized values
    """
    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis=0), 
                             index=df.index, 
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    return(df_qn)







def combined_dfs(files):
    sample_nub = len(files)
    record_col_name=[]
    all_df_list=[]
    for fi in files:
        
         tmp_df=pd.read_csv(fi, sep='\t',header=None)
         tmp_col_name=os.path.basename(fi).split('_')[0] + "_"+os.path.basename(fi).split('_')[1]
         tmp_df.columns=['chrom','starts','ends','name',tmp_col_name ]
         #replace miss value as zero
         tmp_df[tmp_col_name]=tmp_df[tmp_col_name].apply(lambda x: str2num(x,'.'))
         #jbw 2024
         all_df_list.append(tmp_df.copy())
         #print(fi)
         #print(all_df)
         #print(tmp_df)
         #all_df = pd.concat([all_df, tmp_df], ignore_index=True).copy()

         record_col_name.append(tmp_col_name)


    for i in range(1,len(all_df_list)):

        if i==1:
            combined_df=pd.merge(all_df_list[i-1],all_df_list[i],on=['chrom','starts','ends'], how='outer')
        else:
            combined_df=pd.merge(combined_df,all_df_list[i],on=['chrom','starts','ends'] , how='outer')
    filtered_df=combined_df[combined_df[record_col_name].sum(axis=1)>0].copy()
    
    
    new_record_col_name=[ i + '_qnt' for i in record_col_name]
    filtered_df[new_record_col_name ]=quantile_normalize(filtered_df[record_col_name]).copy()
    
    log2qnt_df= np.log2(filtered_df[new_record_col_name]+1)
    new_record_col_name2=[ i + '_Log2zscore' for i in new_record_col_name]
    filtered_df[new_record_col_name2]=log2qnt_df.apply(zscore).copy()
    
    return filtered_df






def select_sample2export(selected_sample, filtered_df,out_path):
  #selected_sample='TIME0-DNas'
  selected_columns=['chrom','starts','ends',selected_sample+'_qnt_Log2zscore',selected_sample+'_qnt',selected_sample]
  out_df1=filtered_df[selected_columns].copy()
  #set format
  out_df1[selected_sample+'_qnt_Log2zscore'] = out_df1[selected_sample+'_qnt_Log2zscore'].map('{:,.3f}'.format)

  out_file=os.path.join(out_path,selected_sample+"_100b.bed")
  out_df1.to_csv(out_file,sep='\t',index=False, header=None)
  sort_cmd = "sort -k1,1V -k2,2n -k3,3n " + out_file +" > "+ out_file.replace(".bed", "_sorted.bed")
  subprocess.run(sort_cmd, shell=True)
  rm = "rm -f " + out_file 
  subprocess.run(rm, shell=True)







def combine_files(in_files,out_dir):
    num_files = len(in_files)
    all_rh_exp = []
    all_ra_head = []
    ids = []

    for i in range(num_files):
        all_ra_head.append(in_files[i])
        with open(in_files[i], 'r') as data_file:
            counts = 0
            for line in data_file:
                line = line.strip().split('\t')
                chr_pos = ':'.join([line[0], line[1], line[2]])
                chr_pos = chr_pos.replace(' ', '').lower()

                if i == 0:
                    ids.append(chr_pos)

                if i == 0:
                    all_rh_exp.append(line[4])
                else:
                    temp_line = all_rh_exp[counts]
                    temp_line = '\t'.join([temp_line, line[4]])
                    all_rh_exp[counts] = temp_line

                counts += 1



    out_filename = f"{out_dir}/SimpleCombined_data.csv"
    with open(out_filename, 'w') as output_data:

        headline = "ID"
        num_head = len(all_ra_head)

        for idx in range(num_head):
            headline += f"\t{all_ra_head[idx]}"

        output_data.write(f"{headline}\n")

        for i in range(len(all_rh_exp)):
            lines = all_rh_exp[i]
            output_data.write(f"{ids[i]}\t{lines}\n")
    
    return f"{out_dir}/SimpleCombined_data.csv"






def filter_combined(combined, out_dir):
    
    df = pd.read_csv(combined, sep="\t")
    
    df.replace('.', 0, inplace=True)
    
    data_df2=df.iloc[:,1:].copy()
    

    data_df2 = data_df2.loc[(data_df2 != 0).any(axis=1)]
    
    percengate_of_non_zeros=data_df2.sum(axis=1)/data_df2.shape[1]
    
    cutoff=1/data_df2.shape[1]
    filtered_data_df=df.loc[percengate_of_non_zeros>cutoff,:].copy()
    
    if df.shape[1] > 5:
        num_occurrences = (data_df2.iloc[:, 1:] > 0).sum(axis=1)
        filtered_data_df=filtered_data_df.loc[num_occurrences > 2]
        
        

    return filtered_data_df
 




   
def make_combined_signal(peaks_files, chromosome_sizes,blacklist_removed_bin,peaks_folder, searchStr1,searchStr2, out_path):
    #jbw 2024
    tmp_out_file=[]
    for file in peaks_files:
        #jbw 2024
        print(file)
        tmp_out_file.append(bin_count_macs2(file,chromosome_sizes,blacklist_removed_bin))

    sorted_enrich = glob.glob(os.path.join(peaks_folder,"*FolderEnrich_sorted.bed"))
    
    filtered_df = combined_dfs(sorted_enrich)
    #jbw 2024
    tmp_out_files=[ os.path.basename(ti) for ti in tmp_out_file if ti is not None] 
    #check selected samples in tmp_out_file
    for si in searchStr1:
        if len([ti for ti in tmp_out_files if si in ti])>0:
           select_sample2export(si, filtered_df,out_path)
        else:
           print("Not find: ", si)

    for si in searchStr2:
        if len([ti for ti in tmp_out_files if si in ti])>0:
           select_sample2export(si, filtered_df,out_path)
        else:
           print("Not find: ",si)
    
    in_files = glob.glob(os.path.join(out_path,"*sorted.bed"))
    
    combined = combine_files(in_files, out_path)
    
    combined_df = filter_combined(combined, out_path)
    
    combined_df.insert(0,'id', combined_df.ID)
    new_df=combined_df.id.str.split(':',expand=True)
    new_df.columns=['chrs','pos_start','pos_end']
    out_df= pd.concat([new_df,combined_df],axis=1).copy()
    new_out_df=out_df.drop("ID", axis=1).copy()
    
    
    return new_out_df


def seacr_enrichement(peaks_folder,blacklist_removed_bin, chromosome_sizes, out_combined_files, searchStr1,searchStr2):
    #jbw 2024 
    tmp_file_list= glob.glob(os.path.join(peaks_folder,"*peaks.*.bed"))
    if len(tmp_file_list)==0:
       tmp_file_list=glob.glob(os.path.join(peaks_folder,"*peaks_sorted.bed"))
    #jbw 2024
    peak_files = tmp_file_list.copy()
    #jbw 2024
    #peaks_files = glob.glob(os.path.join(peaks_folder,"*peaks*.bed"))
    
    out_df = make_combined_signal(peaks_files, chromosome_sizes, blacklist_removed_bin, peaks_folder, searchStr1,searchStr2, out_combined_files)
    
    out_head_file=os.path.join(out_combined_files,'combined_signals_100b.head')
    out_data_file=os.path.join(out_combined_files,'combined_signals_100b.bed')
    
    out_head_df=out_df.iloc[0:1,:].copy()
    out_head_df.to_csv(out_head_file,sep='\t',index=False)
    out_df.to_csv(out_data_file,sep='\t',header=None,index=False)
    cmd = "gzip " + os.path.join(out_combined_files,'combined_signals_100b.bed')
    subprocess.run(cmd, shell=True)
    return



def enrichment_main(peaks_folder,blacklist_removed_bin, chromosome_sizes, out_combined_files, searchStr1,searchStr2):
    #jbw 2024 
    tmp_file_list= glob.glob(os.path.join(peaks_folder,"*peaks.*.bed"))
    if len(tmp_file_list)==0:
       tmp_file_list=glob.glob(os.path.join(peaks_folder,"*peaks_sorted.bed"))
    peak_files_all = tmp_file_list.copy()
    #jbw 2024
    #peak_files_all = glob.glob(os.path.join(peaks_folder,"*peaks*.bed"))
    
    peaks_files = [x for x in peak_files_all if "summitRegion" not in x]
    
    #jbw 2024
    out_df = make_combined_signal(peaks_files, chromosome_sizes, blacklist_removed_bin, peaks_folder, searchStr1,searchStr2, out_combined_files)
    
    out_head_file=os.path.join(out_combined_files,'combined_signals_100b.head')
    out_data_file=os.path.join(out_combined_files,'combined_signals_100b.bed')
    
    out_head_df=out_df.iloc[0:1,:].copy()
    out_head_df.to_csv(out_head_file,sep='\t',index=False)
    out_df.to_csv(out_data_file,sep='\t',header=None,index=False)
    cmd = "gzip " + os.path.join(out_combined_files,'combined_signals_100b.bed')
    subprocess.run(cmd, shell=True)
    return

