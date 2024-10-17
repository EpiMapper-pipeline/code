 #this script is used to find genome position for DAR peaks in TSS, gene, 5distance up and enhancers
import os
import pandas as pd
import multiprocessing as mp
import glob
#exec(open("dar_peaks2genome.py").read())

def initpool(dar_df):
  global global_dar_df
  global_dar_df=dar_df




def find_genome4peaks(args):
  global global_dar_df
  
  in_genome_file=args
  
  genome_df=pd.read_csv(in_genome_file,sep='\t',header=None)
  
  record_genome=[]
  
  for idx, rows in global_dar_df.iterrows():
      
      tmp_df=genome_df[ genome_df[7].apply(lambda x: rows.ids in x )].copy()
      
      if tmp_df.shape[0]>0:
          
        tmp_str='~'.join(tmp_df[3].to_list())
        
      else:
          
        tmp_str=''
        
      record_genome.append(tmp_str)
      
  column_name='_'.join(in_genome_file.split('_0.01_')[-1].split('_')[0:3]) 
  
  print(column_name +  ' - Done' )
  
  out_df=pd.DataFrame(data=record_genome,columns=[column_name],index=global_dar_df.ids) 
  
  return out_df.copy()


def annotation(diff_dir, DAR, in_genome_folder, in_genome_files):
  out_folder=DAR


  in_dar_file= glob.glob(os.path.join(DAR,'combined_peaks_merged_pvals*.csv'))[0]
  in_dar_df=pd.read_csv(in_dar_file,sep='\t',index_col=0)
  in_dar_df.insert(0,'ids',in_dar_df.index)
  #in_dar_df['ids']=in_dar_df.index

  
  all_genome_files=[os.path.join(in_genome_folder, os.path.basename(fi)) for fi in in_genome_files]
  genome_files = [item for item in all_genome_files if "intergenic" not in item]

  num_of_process=len(genome_files)
  #parallel search
  pool= mp.Pool(processes=num_of_process, initializer=initpool,initargs=(in_dar_df ,))
  out_dfs= pool.map(find_genome4peaks,[(genome_files[loop]) for loop in range(0,num_of_process)],1)
  pool.close()

  #merge results from multiple processes
  for oi in out_dfs:
      in_dar_df=in_dar_df.merge(oi, left_index=True, right_index=True,how='left').copy()

  #export
  out_file=in_dar_file.replace('.csv','_genome.csv')
  out_file=out_file.replace("out_combined_files", "DAR")
  out_folder, out_name=os.path.split(out_file)
  if not os.path.exists(out_folder):
      os.makedirs(out_folder)
      
  
  
 
  in_dar_df.to_csv(out_file,sep='\t',index=False)




