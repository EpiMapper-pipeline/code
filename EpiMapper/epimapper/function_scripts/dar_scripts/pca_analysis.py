import numpy as np
import pandas as pd
import logging
import warnings
warnings.filterwarnings("ignore")
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def do_PCA_analysis(all_smooth,gcb_col_idx,tumor_col_idx, tmp_column_names) :
   ''' Perform PCA analysis on all_smooth, gcb_col_idx and tumor_col_idx are column index of gcb and tumor group
       Output: PCA dataframe and column names
       ALL input data are numpy array format
   '''
   #do PCA analysis on data 
   #the first 4 is gcb, and the rest is tumor 
   X= StandardScaler().fit_transform(all_smooth)
   col_name= np.concatenate([tmp_column_names[gcb_col_idx[0]],tmp_column_names[tumor_col_idx[0]]],axis=0)
   #feat_cols= pd.DataFrame(X )

   #PCA analysis
   #print(all_smooth,col_name)
   pca_X=PCA(n_components=3)
   principalComponents_X=pca_X.fit_transform(X)
   pca_df=pd.DataFrame(data=principalComponents_X, columns=['PCA 1','PCA 2','PCA 3'])
   logging.info('Explained variation per principal component: {}'.format(pca_X.explained_variance_ratio_))
   return pca_df, col_name, pca_X



def plot_3D_PCA(fig,sub_plot_num, pca_df,wildType_fileString,num_of_gcb,num_of_tumor, other_group_str ):
   ''' 3-D plot of the first 3 PCAs
   '''
   #first is GCB or wildtype group , then comes Tumor group
   #plot 3D PCA
   ax=fig.add_subplot(sub_plot_num,projection='3d')
   my_color=[1]*num_of_gcb + [2]*num_of_tumor
   ax1= ax.scatter(pca_df['PCA 1'][0:num_of_gcb], pca_df['PCA 2'][0:num_of_gcb], pca_df['PCA 3'][0:num_of_gcb], c='green', s=100,label=wildType_fileString )
   ax2= ax.scatter(pca_df['PCA 1'][num_of_gcb:num_of_tumor+num_of_gcb], pca_df['PCA 2'][num_of_gcb:num_of_tumor+num_of_gcb], pca_df['PCA 3'][num_of_gcb:num_of_tumor+num_of_gcb], c='red', s=100,label=other_group_str)

   # make simple, bare axis lines through space:
   xAxisLine = ((min(pca_df['PCA 1']), max(pca_df['PCA 1'])), (0, 0), (0,0))
   ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'orange')
   yAxisLine = ((0, 0), (min(pca_df['PCA 2']), max(pca_df['PCA 2'])), (0,0))
   ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'orange')
   zAxisLine = ((0, 0), (0,0), (min(pca_df['PCA 3']), max(pca_df['PCA 3'])))
   ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'orange')

   #add label for plotted data points
   if False:
     for idx, row in pca_df.iterrows():
          ax.text(row['PCA 1'], row['PCA 2'], row['PCA 3'], idx)

   # label the axes
   ax.set_xlabel("PC1",fontdict={'fontsize': 16, 'weight': 'bold'})
   ax.set_ylabel("PC2",fontdict={'fontsize': 16, 'weight': 'bold'})
   ax.set_zlabel("PC3",fontdict={'fontsize': 16, 'weight': 'bold'})
   ax.tick_params(labelsize=13) #axis font
   ax.set_title("PCA ",fontsize=21, fontweight='bold')
   ax.legend(fontsize=14)



   return ax
