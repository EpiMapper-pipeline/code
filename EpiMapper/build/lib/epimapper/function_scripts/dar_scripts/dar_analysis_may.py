import pandas as pd
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from random import randint
import argparse
from sklearn.decomposition import PCA
from .pca_analysis import do_PCA_analysis
from .pca_analysis import plot_3D_PCA
from .dmr_map2genome import dmr_map2_genome_main
from .dmr_cal2genome_percent import dmr_cal2genome_percent_main
from .dmr_percent2plot import dmr_percent2plot_main
import logging
import warnings
warnings.filterwarnings("ignore")
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib as mlt
from .find_dar import find_col_index4sample, do_dar_test2
from scipy.stats import ttest_ind, zscore
#for write to file
mlt.use('Agg')
#for show plot
#mlt.use('TkAgg')
#exec(open('dar_analysis.py').read())

def make_pvalue_files(in_dar_file,in_dar_pval_file,pval_cutoff):
  '''
  this script is used to make dmr_analysis format file for genome map
  
  '''

  dar_df=pd.read_csv(in_dar_file,sep='\t', header=None)
  
  dar_df.columns=['chrs','pos_start','pos_end','sample_info']



  dar_pval_df=pd.read_csv(in_dar_pval_file,sep='\t')
  tmp_cols=dar_pval_df.columns.to_list()
  tmp_cols[0]='ID'
  dar_pval_df.columns=tmp_cols
  
  dar_df['ID']=dar_df.chrs+ ':'+dar_df.pos_start.astype(str)+':'+dar_df.pos_end.astype(str)


  #merge two dataframe
  merged_df=pd.merge(dar_df, dar_pval_df, how='left',on='ID')

  merged_df['out_info1']= merged_df['sample_info'].apply(lambda x: str(len(x.split(','))))  
  
  merged_df['out_info2']=  [i for i in range(0,merged_df.shape[0])]
  merged_df['out_info3']= merged_df.chrs+':' + merged_df.pos_start.astype(str) +':'+  merged_df.pos_end.astype(str) + ':mr' + merged_df['out_info2'].astype(str)
  merged_df['out_sample_info']=merged_df.out_info3+ ':' + merged_df.out_info1

  merged_df['delta']=merged_df.tval
  merged_df['isDAR']='U'
  merged_df.loc[merged_df['pval']<pval_cutoff,['isDAR']]='D'
  merged_df['up_or_down']='up'
  merged_df.loc[merged_df.delta<0,['up_or_down']]='down'

  merged_df['out_sample_info2']=merged_df['out_sample_info']+':'+ merged_df['up_or_down']+':'+ merged_df['isDAR']
  out_df=merged_df[['chrs','pos_start','pos_end','out_sample_info2','pval']].copy()

  out_file=in_dar_file.replace('.bed','_pval_min_'+str(pval_cutoff)+'.bed')

  out_df.to_csv(out_file,sep='\t', header=None,index=None)
  return out_file, merged_df, out_df


def make_genome_annotation_files_by_hmst(out_gfolder, in_file, 
           out_file_name, data, sample_names,cutoff):
  ''' use hmst-seq analyzer to map AR or DAR to genomic regions such as tss, tes, gene, 5dist, enhancers'''


  out_folder2genome=os.path.join(out_gfolder,'genome')
  if not os.path.exists(out_folder2genome):
     print('Create , ' , out_folder2genome)
     os.mkdir(out_folder2genome)

  #test jbw
  pval_cutoff=cutoff
  in_sortedDMR_file =os.path.join(out_gfolder,in_file+'.bed')
  in_geneRegion_file = glob.glob(os.path.join(data,'*region*.txt'))[0]
  in_refFlat_file = glob.glob(os.path.join(data,'*refFlat*.bed'))[0]
  in_minimum_overlap4bedtools = 1e-9
  in_outFile_folder =out_folder2genome
  dmr_min_cutoff = pval_cutoff
  dmr_map2_genome_main(in_sortedDMR_file, in_geneRegion_file, in_refFlat_file, in_minimum_overlap4bedtools, in_outFile_folder, dmr_min_cutoff)
  #print ("dmr_map2_genome_main")  


  ##calculate percentage of DMR in annotated genomic regions
  
  in_outFile_folder = os.path.join(out_gfolder,'genome')
  in_fileName_string = in_file 
  in_outFile_name = out_file_name+'_DAR_ttest_pval_'+str(pval_cutoff)+'.csv'
  in_Ls_or_Gt_pval = 0
  in_LogReg_proba = str(pval_cutoff)
  
  dmr_cal2genome_percent_main(in_outFile_folder, in_fileName_string, in_outFile_name, in_Ls_or_Gt_pval, in_LogReg_proba)
  #print ("dmr_cal2genome_percent_main")



  in_countFile_folder = os.path.join(out_gfolder,'genome')
  in_countFile_name = out_file_name+ '_DAR_ttest_pval_'+str(pval_cutoff)+'.csv' 
  names = sample_names
  dmr_percent2plot_main(in_countFile_folder, in_countFile_name, names)
  #print ("dmr_percent2plot_main")



def make_color_array(num_of_color):
  #make color list
  colors = []
  for i in range(num_of_color):
    colors.append('#%06X' % randint(0, 0xFFFFFF))
  return colors



def pie_plot_of_dars(in_folder,file_name_string, regions,
           column_idx2pval,pval_cutoff,out_fig_name):
  ''' column_idx2pval=8, pval_cutoff=0.001'''  
  #in_files=glob.glob(os.path.join('out_signal/genome','combined*.bed'))
  #regions=['tss','tes','gene','5dist','enhancers','intergenic']
  in_files= glob.glob(os.path.join(in_folder,file_name_string))
  
  #plot pie chart for all AR
  pie_plot_tbl = pd.DataFrame(columns=["region","number_of_total_dars","number_of_dars_0.001"])
  pie_plot_tbl["region"] = regions
  pie_plot_tbl = pie_plot_tbl.set_index("region")
  for fi in in_files:
    base_name = list(map(str.lower,fi.split("/")[-1].split("min")[1].split("_")))
    
    for item in pie_plot_tbl.index:
        if item.lower() in base_name:
            tmp_region = item.lower()
        

    tmp_df=pd.read_csv(fi,sep='\t',header=None)
    
    tmp_df_clean=  tmp_df.drop_duplicates(subset=[7])
    
    tmp_num=tmp_df_clean.shape[0]
    if not pd.isna(pie_plot_tbl.at[tmp_region,"number_of_total_dars"]):
        pie_plot_tbl.loc[tmp_region,"number_of_total_dars"] += tmp_num
    else:
        pie_plot_tbl.loc[tmp_region,"number_of_total_dars"] = tmp_num

    tmp_df_filtered=tmp_df_clean[tmp_df_clean.iloc[:,column_idx2pval]<pval_cutoff]
  
    if not pd.isna(pie_plot_tbl.at[tmp_region,"number_of_dars_0.001"]):
        pie_plot_tbl.loc[tmp_region, "number_of_dars_0.001"]= tmp_df_filtered.shape[0]
    else:
        pie_plot_tbl.loc[tmp_region,"number_of_dars_0.001"] = tmp_df_filtered.shape[0]
    

  colors=make_color_array(len(regions))
  
  plt.figure(figsize=(14, 7))
  plt.subplot(1,2,1)
  plt.pie(pie_plot_tbl["number_of_total_dars"],colors=colors,labels=pie_plot_tbl.index, textprops={'fontsize': 15}) #fontsize
  plt.title("Number of total DARs", fontsize=18, fontweight='bold')
  
  #test jbw
  plt.subplot(1,2,2)
  plt.pie(pie_plot_tbl["number_of_dars_0.001"],colors=colors,labels=pie_plot_tbl.index, textprops={'fontsize': 15}) #fontsize
  plt.title("Number of DARs (p-value < "+ str(pval_cutoff) +")", fontsize=18, fontweight='bold')
  plt.savefig(out_fig_name,dpi=150)

  # save png
  out_fig_name2 = out_fig_name.replace('.pdf','.png')
  plt.savefig(out_fig_name2,format='png',dpi=300)
  print("finished pie_plot_of_dars.png")
  
  return 


def plot_pca4samples(head_file,data_file, normal_group_str,tumor_group_str,out_fig_name):
  '''PCA plot for samples based on overall signals'''
  #file1='out_signal/200bp/combined_signals.head'
  #file2='out_signal/200bp/combined_signals.bed.gz'
  file1=head_file
  file2=data_file

  file_head_df=pd.read_csv(file1,sep='\t')
  file_data_df=pd.read_csv(file2,sep='\t',compression='gzip', header=None)

  #find vat and sat column index
  #str1=['N_rep']
  #str2=['TL_rep']
  str1=normal_group_str
  str2=tumor_group_str
  sat_idx, vat_idx=find_col_index4sample(file_head_df,str1,str2) 
  out_df=pd.concat([file_data_df.iloc[:,3], file_data_df.iloc[:,sat_idx], file_data_df.iloc[:,vat_idx]],axis=1)

  #extract data
  out_idx=np.concatenate((np.array([3]),sat_idx,vat_idx))
  out_df.columns=file_head_df.columns[out_idx]

  #remove the whole row with 0 elements
  print(out_df)
  out_df2=out_df[out_df.sum(axis=1)!=0].copy()
  out_data=out_df2.iloc[:,1:].to_numpy()

  #prepare for PCA analysis
  len_gcb=len(sat_idx)
  len_tumor=len(vat_idx)
  gcb_col_idx=(np.array([i for i in range(0,len_gcb)]),)
  tumor_col_idx=(np.array([i for i in range(len_gcb,len_gcb+len_tumor)]),)
  tmp_column_names=out_df2.columns.values[1:]
  tmp_column_names2=np.array([os.path.basename(i).split("100")[0] for i in tmp_column_names ],dtype=str)
  #print(out_data.T, gcb_col_idx,tumor_col_idx, tmp_column_names2)
  pca_df, col_name, pca_X =do_PCA_analysis(out_data.T,gcb_col_idx,tumor_col_idx, tmp_column_names2)

  #plot PCA
  fig=plt.figure(figsize=(8,7))
  sub_plot_num=111
  wildType_fileString="Group A"  #'H3K27me3'
  other_group_str = "Group B" #'H3K4me3'

  #mlt.use('TkAgg')

  #in PCA plot, first start wildtype group then is the tumor group
  tmp_ax=plot_3D_PCA(fig,sub_plot_num, pca_df,wildType_fileString,len_gcb,len_tumor, other_group_str)
  for idx, row in pca_df.iterrows():
    tmp_ax.text(row['PCA 1'], row['PCA 2'], row['PCA 3'], tmp_column_names2[idx])
  #plt.show()
  
  plt.savefig(out_fig_name,dpi=150)

  # save png
  out_fig_name2 = out_fig_name.replace('.pdf','.png')
  plt.savefig(out_fig_name2,format='png',dpi=300)
  print("finished plot_3D_PCA.png")
  return out_fig_name, pca_df, col_name





