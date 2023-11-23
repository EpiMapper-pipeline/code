

import glob
import os
import pandas as pd
import argparse


def main(out_file, out_dmr_files,count_data_df, is_less_or_great_pval,logProb_cutoff):

  logReg_proba_cutoff=float(logProb_cutoff)

  for fil in out_dmr_files:
    
    #print(fil)
    in_file_df=pd.read_csv(fil,header=None, sep='\t')
    in_file_df.columns=['chrs','start_pos','end_pos','genome_info','mr_chrs','mr_start_pos','mr_end_pos','mr_info','mr_logReg_proba']

    #count mr and dmr in the region
    #genome_name=fil.split('_')[-4].lower()
    #print(fil)
    if '_gene_u' in fil.lower():
      genome_name='gene'
    elif '_tss_u' in fil.lower():
      genome_name='tss'
    elif '_tes_u' in fil.lower():
      genome_name='tes'
    elif '_5dist_u' in fil.lower():
      genome_name='5distUp'
    elif '_intergenic_u' in fil.lower():
      genome_name='intergenic'
    elif '_enhancer_' in fil.lower() or '_enhancers_' in fil.lower():
      genome_name='enhancer'
    elif '_5dist_d' in fil.lower():
      genome_name='5distDown'

    #input('Click ')
    uq_genome_info=in_file_df.genome_info.unique().shape
    uq_mr_info=in_file_df.mr_info.unique().shape

    uq_genome_info_df, uq_mr_info_df, uq_dmr_info_df=[],[],[]
    uq_genome_info_df=in_file_df.drop_duplicates(subset=['genome_info'])
    uq_mr_info_df=in_file_df.drop_duplicates(subset=['mr_info'])
    uq_dmr_info_df=uq_mr_info_df.copy()
    if is_less_or_great_pval==1:
       uq_dmr_info_df=uq_dmr_info_df[uq_dmr_info_df.mr_logReg_proba>=logReg_proba_cutoff]
    else:
       uq_dmr_info_df=uq_dmr_info_df[uq_dmr_info_df.mr_logReg_proba<logReg_proba_cutoff]

    
    total_dmr=uq_dmr_info_df.shape[0]
    total_hyper=uq_dmr_info_df[uq_dmr_info_df.mr_info.str.contains('up')].shape[0]
    total_hypo=uq_dmr_info_df[uq_dmr_info_df.mr_info.str.contains('down')].shape[0]
    

    count_data_df.at[genome_name,'Total']=total_dmr
    count_data_df.at[genome_name,'up']=total_hyper
    count_data_df.at[genome_name,'down']=total_hypo

  #print('Export at: ', out_file)
  count_data_df.to_csv(out_file,sep='\t')

def dmr_cal2genome_percent_main(in_outFile_folder,in_fileName_string,in_outFile_name, in_Ls_or_Gt_pval, in_LogReg_proba):
  out_folder=in_outFile_folder
  file_string=in_fileName_string
  out_file_name=in_outFile_name
  out_file=os.path.join(out_folder,out_file_name)

  is_less_or_great_pval=in_Ls_or_Gt_pval
  logProb_cutoff=in_LogReg_proba

  out_dmr_file_path=os.path.join(out_folder,file_string+'_*.bed')
  out_dmr_files=glob.glob(out_dmr_file_path)
  #print(out_dmr_file_path)

  count_data_df=pd.DataFrame(columns=['Total','up','down'],index=['gene','tss','tes','5dist','intergenic','enhancer'])
  main(out_file, out_dmr_files,count_data_df, is_less_or_great_pval,logProb_cutoff)

 

 


