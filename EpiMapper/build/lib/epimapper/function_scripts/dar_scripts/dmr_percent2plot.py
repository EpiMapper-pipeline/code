#bar plot of percentage of DMR in different genomic regions
import pandas as pd
import os
import numpy as np
import matplotlib as mlt
import matplotlib.pyplot as plt

mlt.use('Agg')
#mlt.use('TkAgg')


def main(count_file_folder, count_file_name,names):
  count_file=os.path.join(count_file_folder,count_file_name)
  if not os.path.exists(count_file):
     print('Input count table file is not found please check input file name or path, I stop !', count_file)
     exit(1)
  count_data_df=pd.read_csv(count_file,index_col=0,sep='\t')
  count_data_matrix=count_data_df.to_numpy()
  percent_up=count_data_matrix[:,1]/count_data_matrix[:,0]*100
  percent_down=count_data_matrix[:,2]/count_data_matrix[:,0] *100
  data_size=count_data_matrix[:,1].shape[0]
  title_name =names[0] + " vs " +names[1]
  percent_data_matrix=np.concatenate((percent_up.reshape(1,data_size),percent_down.reshape(1,data_size)),axis=0)
  percent_data_df=pd.DataFrame(data=percent_data_matrix.T,columns=['up_percent','down_percent'],index=count_data_df.index.to_list())
  percent_data_df['genome']=percent_data_df.index.to_list()
  ax=percent_data_df.plot.bar(x='genome',y=['up_percent','down_percent'], title=title_name ,figsize=(18,24)) # font size
  ax.set_title(title_name, fontsize=42, fontweight='bold')
  ax.tick_params(axis='x', labelsize=32) #axis font
  ax.tick_params(axis='y', labelsize=32)#axis font
  ax.set_xlabel('',fontdict={'fontsize': 37, 'weight': 'bold'}) # ylabel font size
  ax.set_ylabel('Percentage',fontdict={'fontsize': 37, 'weight': 'bold'}) # ylabel font size
  ax.legend(fontsize=37) # legend font size
  out_fig_file=count_file.replace('.csv','.pdf')
  out_fig_file=out_fig_file.replace("differential_analysis/out_combined_files/genome","summary_tables")

  #out_fig_file=out_fig_file.replace('counts','percent')
  #print(out_fig_file)

  plt.savefig(out_fig_file, format='pdf')

  # save png
  out_fig_file2=out_fig_file.replace('.pdf','.png')
  plt.savefig(out_fig_file2, format='png', dpi=300)
  print("finished Epimapper_DAR_ttest_pval_0.01.png")
  #plt.show()

def dmr_percent2plot_main(in_countFile_folder,in_countFile_name, names):
  count_file_folder=in_countFile_folder
  count_file_name= in_countFile_name
  main(count_file_folder, count_file_name, names)



