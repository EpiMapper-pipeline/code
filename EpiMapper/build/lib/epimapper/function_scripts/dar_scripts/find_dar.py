#this script is used to extract eack called peak region from merged peaks and show or do DAR analysis
import pandas as pd
from scipy.stats import norm, ttest_ind, ks_2samp, mannwhitneyu, ranksums 
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
#exec(open("find_dar_nature.py").read())

def rratios(tmp_x,tmp_y):
  '''compute relative rations between two samples'''
  if tmp_x>0 and tmp_y>0:
    rr=(tmp_x-tmp_y)/(tmp_x+tmp_y)*2
  else:
    rr=(np.exp(tmp_x)-np.exp(tmp_y))/(np.exp(tmp_x)+np.exp(tmp_y))*2
  #scale=0.2223 representative relative ratio=0.2223 and real ratios=1.25
  #print(tmp_x,tmp_y,rr)
  pval=norm.pdf(abs(rr),loc=0,scale=0.222223)*2
  return rr, pval





def do_dar_test2(test_methods,tmp_xx,tmp_yy):
  '''input two numpy array for test'''
  #take mean of all samples
  tmp_x=tmp_xx.astype(float).mean(axis=1).to_numpy()
  tmp_y=tmp_yy.astype(float).mean(axis=1).to_numpy()
  
  if len(tmp_x)>1:
     #peaks with multiple data points
     if test_methods == 'kstest':
      tmp_p=ks_2samp(tmp_x,tmp_y)
     elif test_methods == 'mannwhitneyu':
      tmp_p=mannwhitneyu(tmp_x,tmp_y)
     elif test_methods == 'ranksumtest':
      tmp_p=ranksums(tmp_x,tmp_y)
     else: # ttest
      tmp_p=ttest_ind(tmp_x,tmp_y)
     max_x=np.max(tmp_x)
     max_y=np.max(tmp_y)
     num_x=len(tmp_x)
     num_y=len(tmp_y)
  elif tmp_x.size==1 :
     #peaks with single data points
     tmp_x=tmp_xx.astype(float).to_numpy()[0]
     tmp_y=tmp_yy.astype(float).to_numpy()[0]
     if test_methods == 'kstest':
      tmp_p=ks_2samp(tmp_x,tmp_y)
     elif test_methods == 'mannwhitneyu':
      tmp_p=mannwhitneyu(tmp_x,tmp_y)
     elif test_methods == 'ranksumtest':
      tmp_p=ranksums(tmp_x,tmp_y)
     else:
      tmp_p=ttest_ind(tmp_x,tmp_y)
     max_x=np.max(tmp_x)
     max_y=np.max(tmp_y)
     #print(tmp_x, tmp_y)
     #tmp_p=[0,0]
     #tmp_p[0],tmp_p[1]= rratios(tmp_x[0],tmp_y[0])
     #max_x=tmp_x[0]
     #max_y=tmp_y[0]
     num_x=1
     num_y=1
  else:
     #peaks with empty array
     tmp_p=[1,1]
     max_x=0
     max_y=0
     num_x=0
     num_y=0
  #tmp include both tval and pval 
  return tmp_p,max_x,max_y,num_x, num_y 

def initpool(peakSignal_df):
  global global_peakSignal_df
  global_peakSignal_df=peakSignal_df


def do_DAR_test(args):
 ''' do ttest sat_idx/X vs. vat_idx/Y'''
 global global_peakSignal_df
 in_tmp_ar_ids,sat_idx,vat_idx,test_methods= args
 all_p={}
 for ii in in_tmp_ar_ids:
     
   tmp_df=global_peakSignal_df[global_peakSignal_df.id==ii].copy()
   tmp_x=tmp_df.iloc[:,sat_idx-3+7]
   tmp_y=tmp_df.iloc[:,vat_idx-3+7]
   tmp_p,max_x,max_y, num_x,num_y=do_dar_test2(test_methods,tmp_x,tmp_y)
   all_p[ii]=[tmp_p[0],tmp_p[1],max_x,max_y,num_x,num_y]
 return all_p




def parallel_do_DAR(args):
  ''' sat_idx for group1 , vat_idx for group2
      do ttest for group1 vs. group2
  '''
  #do parallel computation for AR
  all_p=[]
  #num_of_process=15
  num_of_process, in_ar_df, in_peakSignal_df, sat_idx, vat_idx, test_methods= args

  #divide data to multiple processes
  len_of_data=len(in_ar_df.id.to_list())
  interval=int(len_of_data/num_of_process)
  index2data=[]
  tmp_index=[]
  for i in range(0,len_of_data,interval):
     tmp_index.append(i)

  if tmp_index[-1]<=len_of_data:
     tmp_index[-1]=len_of_data+1
  if len(tmp_index)==num_of_process:
     num_of_process -=1

  #do parallel calculation
  ar_ids=in_ar_df.id.to_numpy()
  pool= mp.Pool(processes=num_of_process, initializer=initpool,initargs=(in_peakSignal_df ,))
  out_pvals= pool.map(do_DAR_test,[(ar_ids[tmp_index[loop]:tmp_index[loop+1]], sat_idx, vat_idx, test_methods) for loop in range(0,num_of_process)],1)
  pool.close()
  out_pvals= list(filter(None, out_pvals))

  #collect returned data from all processes
  all_p={}
  for dict_p in out_pvals:
     all_p.update(dict_p)

  #find DARs with simple filtering
  test_p=all_p.copy()
  for dd in test_p.keys():
    if (test_p[dd][1]>=0.05 or ( test_p[dd][2] <1 and test_p[dd][3]<1 )):
       del all_p[dd]
  print(len(all_p))

  return test_p, all_p



def find_indices(search_strings,tmp_cols):
    indices = []
    for search_string in search_strings:
        search_indices = np.where(np.char.find(tmp_cols, search_string) >= 0)[0]
        indices.extend(search_indices)
    return np.unique(indices)


def find_col_index4sample(in_ar_head_df, searchStr1, searchStr2):
    '''sample1_idx for searchStr1 and sample2_idx for searchStr2 '''
    tmp_cols = in_ar_head_df.columns.values.astype(str)
    
    
    sample1_idx = find_indices(searchStr1,tmp_cols)
    sample2_idx = find_indices(searchStr2,tmp_cols)
    
    return sample1_idx, sample2_idx



