#map DMR to genome reference file
import os
import pandas as pd
import argparse
import shutil



def count_dmrs_not_mapped2genome(in_sorted_dmr_file,record_out_files,dmr_min_cutoff, is_greater=True):
  #count how many DMRs are not mapped to annotated geneomic regions
  #coumt MR or DMR in genomic files
  all_data_df=[]
  for fil in record_out_files:
     tmp_data_df=pd.read_csv(fil,header=None, sep='\t')
     #jbw 2024
     #all_data_df.append(tmp_data_df.copy())
     all_data_df = pd.concat([all_data_df, tmp_data_df], ignore_index=True).copy()

  lines=0
  for i in all_data_df:
     lines+= len(i)

  all_indata_df=pd.concat(all_data_df)
  #if all_indata_df.shape[1]==9:
  all_indata_df.columns=['chrs','start_pos','end_pos','genome_info','mr_chrs','mr_start_pos','mr_end_pos','mr_sites','mr_logReg_proba']
  #elif all_indata_df.shape[1]=8:
  #   all_indata_df.columns=['chrs','start_pos','end_pos','genome_info','mr_chrs','mr_start_pos','mr_end_pos','mr_sites']

  uq_indata_df=all_indata_df.drop_duplicates(subset=['mr_sites'])
  total_uq_mrs=uq_indata_df
  min_cutoff=dmr_min_cutoff
  
  if is_greater:
    print('DMR selection based on >= '+ str(min_cutoff))
    total_uq_dmrs=uq_indata_df[uq_indata_df['mr_logReg_proba']>=min_cutoff]
  else:
    print('DMR selection based on <= ' + str(min_cutoff))
    total_uq_dmrs=uq_indata_df[uq_indata_df['mr_logReg_proba']<=min_cutoff]

  in_dmrs_df=pd.read_csv(in_sorted_dmr_file,header=None,sep='\t')
  in_dmrs_df.columns=['mr_chrs','mr_start_pos','mr_end_pos','mr_sites','mr_logReg_proba']
  total_in_mrs=in_dmrs_df
  
  if is_greater:
    total_in_dmrs=in_dmrs_df[in_dmrs_df['mr_logReg_proba']>=min_cutoff]
  else:
    total_in_dmrs=in_dmrs_df[in_dmrs_df['mr_logReg_proba']<=min_cutoff]

  print('Number of MR or DMR do not find mapped genome information')
  print(total_uq_mrs.shape[0]-total_in_mrs.shape[0])
  print(total_uq_dmrs.shape[0]-total_in_dmrs.shape[0])

  print('Perentage of MR or DMR mapped to genome info')
  print(total_uq_mrs.shape[0]/total_in_mrs.shape[0])
  print(total_uq_dmrs.shape[0]/total_in_dmrs.shape[0])

  diff_in_mr=set(total_in_mrs.mr_sites.to_list())- set(total_uq_mrs.mr_sites.to_list()) 
  diff_in_dmr=set(total_in_dmrs.mr_sites.to_list())- set(total_uq_dmrs.mr_sites.to_list())

  #check the distribution of unmapped MRs in chromes
  chrs=[]
  for i in range(1,25):
    if i<23:
      chrs.append('chr'+str(i))
    elif i==23:
      chrs.append('chrX')
    elif i==24:
      chrs.append('chrY')

  diff_in_dmr_df=pd.DataFrame(data=list(diff_in_dmr),columns=['dmr_sites'])
  dict_chr={}
  total=0
  for ii in chrs:
    dict_chr[ii]=diff_in_dmr_df[diff_in_dmr_df.dmr_sites.str.contains(ii+':')]
    total += dict_chr[ii].shape[0]
  return total, dict_chr






def prepare_result_folder(res_folder):
    if os.path.exists(res_folder):
       print("The result folder, ", res_folder, " , exits ")
    else:
       print("Create result folder, ", res_folder)
       os.mkdir(res_folder)

   
def clear_folder(folder):
    if os.path.exists(folder):
        shutil.rmtree(folder)

    os.mkdir(folder)


def main(region_files, methylation_file, reference_file, dmr_min_cutoff,out_folder,min_overlap,is_greater=True):
  #for each region to find its DMR or MRs
  #there is a bug in recent bedtools but bedtools2.2 version works, now use *_range instead of *_all file for regions overlapping, which works in both old and new version of bedtools
  record_out_files=[]
  for fil in region_files:
    region_file=fil
    region_name=os.path.basename(fil).split('_')[0].lower()
    
    out_methylation_name = out_folder + '/' + os.path.basename(methylation_file)[:-4]
    out_region_name = os.path.basename(region_file)[:-4]


    dist5_methylation_file = out_methylation_name + '_' + 'noGenes.bed'
    

    #print(out_methylation_name,out_region_name)
    if region_name == '5dist':
        # For 5distance we first remove genes(TSS, geneBody and TES) from the two methylation files
        cmd='bedtools intersect -a ' + methylation_file + ' -b ' + reference_file + \
                      ' -v > ' + dist5_methylation_file

        results=os.system(cmd)

        if results !=0:
           print('Error in 5dist bedtools, ', cmd )
           exit(1)

        out = dist5_methylation_file[:-4] + '_' + out_region_name + '_overlap' + str(min_overlap) + '.bed'

        cmd='bedtools intersect -a ' + region_file + ' -b ' + \
                  dist5_methylation_file + ' -wa -wb -f ' + str(min_overlap) + ' > ' + out
        results=os.system(cmd)
        if results !=0:
           print('Error in bedtools intersect, ', cmd)
           exit(1)

        os.system('rm ' + dist5_methylation_file)  # removes temporary file
    else:
        out = out_methylation_name + '_' + out_region_name + '_overlap' + str(min_overlap) + '.bed' 
        cmd= 'bedtools intersect -a ' + region_file + ' -b ' + methylation_file + \
                  ' -wa -wb -f ' + str(min_overlap) + ' > ' + out
        results=os.system(cmd)
        if results != 0:
           print('Error in bedtools intersect, ', cmd )
           exit(1)


  
    record_out_files.append(out)


  return 






def dmr_map2_genome_main(in_sortedDMR_file,in_geneRegion_file, in_refFlat_file, in_minimum_overlap4bedtools, in_outFile_folder, dmr_min_cutoff ):
  in_sorted_dmr_file=in_sortedDMR_file
  in_region_files=in_geneRegion_file
  reference_file=in_refFlat_file
  min_overlap=in_minimum_overlap4bedtools
  out_folder=in_outFile_folder

  prepare_result_folder(out_folder)

  region_files=pd.read_csv(in_region_files,header=None)
  region_files=region_files.loc[:,0].to_list()
  methylation_file=in_sorted_dmr_file

  if dmr_min_cutoff ==0 :
    
    parameter_file=in_sorted_dmr_file.replace('.bed','_parameter.bed')
    print(parameter_file)
    parameter_df=pd.read_csv(parameter_file,sep='\t',header=None)
    tmp_str=parameter_df[0].to_list()[-1]
    print(tmp_str)
    tmp_str=tmp_str.split('_')[-1]
    dmr_min_cutoff=float(tmp_str)
    main(region_files, methylation_file, reference_file, dmr_min_cutoff,out_folder,min_overlap)
    
  else:
    #use <= cutoff value for selecting DMRs
    dmr_min_cutoff=dmr_min_cutoff
    is_greater=False
    main(region_files, methylation_file, reference_file, dmr_min_cutoff,out_folder,min_overlap,is_greater)







