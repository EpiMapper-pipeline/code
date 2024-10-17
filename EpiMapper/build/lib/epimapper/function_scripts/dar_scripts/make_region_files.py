#!/usr/bin/env python3
# -*- coding: utf-8 -*-




import pandas as pd
import os
import subprocess
def make_names(chr, start, end, rna_name, x, y, reg, gene_name, strand,
               orig_start, orig_end):
    """Creates the string names and returns it. The format of names is:
    chr:start_pos:end_pos:rna_name||region_name:X:Y||gene_name:strand:gene_start:gene_end

    Keyword arguments:
    row -- row of dataframe passed to make_TSS_TES_gene_5dist()
    start -- start position
    end -- end position
    x -- X
    y -- Y
    name -- region name"""

    names = ('%s:%s:%s:%s||%s:%d:%d||%s:%c:%s:%s' % (
             chr,                     #chr
             start,                   #start_pos depending on type(TSS, TES, gene, 5dist)
             end,                     #end_pos of type(TSS, TES, gene, 5dist)
             rna_name,                #rna_name
             reg,                     #TSS, TES, gene or 5dist
             x,                       #X
             y,                       #Y
             gene_name,               #gene_name
             strand,                  #strand
             orig_start,              #start_pos
             orig_end))               #end_pos
    return names


def make_BED_of_ref(df):
    """Creates BED formatted dataframe from input file and returns the dataframe

    Keyword arguments:
    df -- dataframe of reference"""

    df['gene_rna_name'] = df['gene_name'].astype(str) + ':' + df['rna_name'].astype(str)
    df['dot'] = '.'
    df_BED = pd.DataFrame(df[['chr', 'start_pos', 'end_pos', 'gene_rna_name',
                              'dot', 'strand']])
    return df_BED



def find_uniq_genes(df, out_folder):
    """Filters and orders data and writes to file

    Keyword arguments:
    df -- pandas dataframe
    num_chr_name -- numerical chromosome name?
    human -- data from human? if False -> mouse or rat
    remove_mir -- remove genes with names starting with mir? (default True)
    """
    df_copy = df.copy()

    

    #removes rows where chr contains '_':
    df_copy = df_copy[df_copy.chr.astype(str).str.contains('_') == False]

    #removes duplicates
    df_copy = df_copy.drop_duplicates(subset=['chr', 'strand', 'TSS', 'TES'])
    
    df_copy = df_copy[['chr', 'TSS', 'TES','strand', 'gene_name', 'rna_name']]

    df_copy.to_csv(os.path.join(out_folder,"tmp_ref.csv"), header=False, index=False, sep="\t")
    
    cmd_sort= "sort -k1,1V -k2,2n -k3,3 " + os.path.join(out_folder,"tmp_ref.csv")+ " > " +os.path.join(out_folder,"tmp_ref_sorted.csv")
    
    subprocess.run(cmd_sort, shell=True)
    
    in_pd = pd.read_csv(os.path.join(out_folder,"tmp_ref_sorted.csv"),sep = "\t", names=['chr', 'TSS', 'TES','strand', 'gene_name', 'rna_name'])
    
    chrm = in_pd[in_pd["chr"].astype(str).str.contains("M", na=False)]

    in_pd = in_pd.drop(list(chrm.index))

    #jbw 2024 
    #new_pd = in_pd.append(chrm)
    new_pd = pd.concat([in_pd, chrm], ignore_index=True).copy()


    new_pd["chr"] =  new_pd["chr"].apply(lambda x: 'chr' + str(x) if not str(x).startswith('chr') else str(x))
    
    rm_cmd = "rm "+ os.path.join(out_folder,"tmp*.csv")

    subprocess.run(rm_cmd, shell=True)
    
    return new_pd


def combine_non_unique(df):
    """Combines similar genomic regions starting and ending at the same site, but
    having different rna_names"""
    pd.options.display.max_colwidth = 100
    #print(df)
    groups = df.groupby(['chr', 'start_pos', 'end_pos', 'x', 'y', 'gene_name', 'reg'], sort=False)
    all_rows = []
    for g in groups:
        group_df = g[1].reset_index()

        chr = group_df.chr[0]
        start_pos = group_df.start_pos[0]
        end_pos = group_df.end_pos[0]
        X = group_df.x[0]
        Y = group_df.y[0]
        gene_name = group_df.gene_name[0]
        reg = group_df.reg[0]
        strand = group_df.strand[0]

        if len(group_df) > 1:
            rna_name = '&'.join(group_df.rna_name.tolist())
            orig_start = group_df.orig_start.drop_duplicates().astype(str).tolist()
            orig_end = group_df.orig_end.drop_duplicates().astype(str).tolist()
            if len(orig_start) > 1:
                orig_start = '&'.join(orig_start)
            else:
                orig_start = orig_start[0]
            if len(orig_end) > 1:
                orig_end = '&'.join(orig_end)
            else:
                orig_end = orig_end[0]
        else:
            rna_name = group_df.rna_name[0]
            orig_start = group_df.orig_start[0].astype(str)
            orig_end = group_df.orig_end[0].astype(str)

        names = make_names(chr, start_pos, end_pos, rna_name, X, Y, reg, gene_name,
                           strand, orig_start, orig_end)
        new_row = [chr, start_pos, end_pos, names]
        all_rows.append(new_row)

    new_df = pd.DataFrame(all_rows, columns=['chr', 'start_pos', 'end_pos', 'names'])
    return new_df



def make_TSS_TES_gene_5dist(df, X, Y, M, N, rem):
    """Creates four BED files: TSS, TES, gene and 5dist from the reference

    df -- sorted dataframe, columns=['chr', 'start_pos', 'end_pos', 'gene_name:rna_name',
                                     'dot', 'strand']
    X -- number of upstream bp, TSS, TES, gene
    Y -- number of downstream bp, TSS, TES, gene
    M -- number of bp from gene start site, 5dist
    N -- number of bp from gene start site, 5dist
    rem -- True if regions should be removed, False if TSS, TES and 5dist should be kept
    """
    removed = []
    #added jbw
    TSS = []; TES = []; gene = []; dist5 = []; dist5D=[] 
    size = len(df.values)
    
    #test jbw 06.23
    #print(df)
    #print(X,Y,M,N, rem)
    #end test
    for i, row in enumerate(df.values):
        gene_start_pos = int(row[1])+Y
        gene_end_pos = int(row[2])-Y
        rna_gene = row[3].split(':')
        rna_name = rna_gene[1]
        gene_name = rna_gene[0]
        if gene_start_pos < gene_end_pos:
            gene.append({'chr': row[0],
                    'start_pos': gene_start_pos,
                    'end_pos': gene_end_pos,
                    'x': X,
                    'y': Y,
                    'gene_name': gene_name,
                    'rna_name': rna_name,
                    'orig_start': row[1],
                    'orig_end': row[2],
                    'strand': row[5],
                    'reg':'gene'})
                    #'names': make_names(row, gene_start_pos, gene_end_pos,
                    #                    X, Y, 'gene')})

            TSS_start_pos = row[1]-X if row[5]=='+' else row[2]-Y
            TSS_end_pos = row[1]+Y if row[5]=='+' else row[2]+X
            TSS.append({'chr': row[0],
                        'start_pos': TSS_start_pos,
                        'end_pos': TSS_end_pos,
                        'x': X,
                        'y': Y,
                        'gene_name': gene_name,
                        'rna_name': rna_name,
                        'orig_start': row[1],
                        'orig_end': row[2],
                        'strand': row[5],
                        'reg':'TSS'})
                        #'names': make_names(row, TSS_start_pos, TSS_end_pos, X, Y,
                        #                    'TSS')})

            #added wang
            if TSS_start_pos<1:
               TSS_start_pos=1

            TES_start_pos = row[1]-X if row[5]=='-' else row[2]-Y
            TES_end_pos = row[1]+Y if row[5]=='-' else row[2]+X
            TES.append({'chr': row[0],
                        'start_pos': TES_start_pos,
                        'end_pos': TES_end_pos,
                        'x': X,
                        'y': Y,
                        'gene_name': gene_name,
                        'rna_name': rna_name,
                        'orig_start': row[1],
                        'orig_end': row[2],
                        'strand': row[5],
                        'reg':'TES'})
                        #'names': make_names(row, TES_start_pos, TES_end_pos, X, Y,
                        #                    'TES')})
            #added wang
            if TES_start_pos<1:
               TES_start_pos=1

            #5'distance upstream
            dist5_start_pos = row[1]-N if row[5]=='+' else row[2]+M
            dist5_end_pos = row[1]-M if row[5]=='+' else row[2]+N
            if dist5_start_pos < 1:
                dist5_start_pos = 1

            #added wang
            if dist5_end_pos <1 :
                dist5_end_pos =1

            if not (dist5_end_pos==1 & dist5_start_pos==1) :
               dist5.append({'chr': row[0],
                        'start_pos': dist5_start_pos,
                        'end_pos': dist5_end_pos,
                        'x': M,
                        'y': N,
                        'gene_name': gene_name,
                        'rna_name': rna_name,
                        'orig_start': row[1],
                        'orig_end': row[2],
                        'strand': row[5],
                        'reg':'5dist'})
                        #'names': make_names(row, dist5_start_pos, dist5_end_pos, M,
                        #                    N, '5dist')})

            #added jbw for 5distD Downstream
            dist5D_start_pos=row[1] +M  if row[5]=='+' else row[2] -N
            dist5D_end_pos=row[1] + N if row[5]=='+' else row[2] -M
            if dist5D_start_pos <1:
               dist5D_start_pos =1

            if dist5D_end_pos <1:
               dist5D_end_pos=1

            if not (dist5D_end_pos==1 & dist5D_start_pos==1):
               dist5D.append( {'chr': row[0],
                        'start_pos': dist5D_start_pos,
                        'end_pos': dist5D_end_pos,
                        'x': M,
                        'y': N,
                        'gene_name': gene_name,
                        'rna_name': rna_name,
                        'orig_start': row[1],
                        'orig_end': row[2],
                        'strand': row[5],
                        'reg':'5distD'})
        else:
            if rem: # remove TSS, TES, 5dist regions even though geneBody too small
                removed.append({'chr': row[0],
                        'start_pos': row[1],
                        'end_pos': row[2],
                        'x': 0,
                        'y': 0,
                        'gene_name': gene_name,
                        'rna_name': rna_name,
                        'orig_start': row[1],
                        'orig_end': row[2],
                        'strand': row[5],
                        'reg':'removed_from_geneBody_TSS_TES_5dist'})
                        #'names': make_names(row, row[1], row[2], 0, 0,
                        #                    'removed_from_geneBody_TSS_TES_5dist')})
            else:
                removed.append({'chr': row[0],
                        'start_pos': row[1],
                        'end_pos': row[2],
                        'x': 0,
                        'y': 0,
                        'gene_name': gene_name,
                        'rna_name': rna_name,
                        'orig_start': row[1],
                        'orig_end': row[2],
                        'strand': row[5],
                        'reg':'removed_only_from_geneBody'})
                        #'names': make_names(row, row[1], row[2], 0, 0,
                        #                    'removed_only_from_geneBody')})

                TSS_start_pos = row[1]-X if row[5]=='+' else row[2]-Y
                TSS_end_pos = row[1]+Y if row[5]=='+' else row[2]+X
                TSS.append({'chr': row[0],
                            'start_pos': TSS_start_pos,
                            'end_pos': TSS_end_pos,
                            'x': X,
                            'y': Y,
                            'gene_name': gene_name,
                            'rna_name': rna_name,
                            'orig_start': row[1],
                            'orig_end': row[2],
                            'strand': row[5],
                            'reg':'TSS'})
                            #'names': make_names(row, TSS_start_pos, TSS_end_pos, X, Y,
                            #                    'TSS')})

                #added wang
                if TSS_start_pos<1:
                   TSS_start_pos=1

                TES_start_pos = row[1]-X if row[5]=='-' else row[2]-Y
                TES_end_pos = row[1]+Y if row[5]=='-' else row[2]+X
                TES.append({'chr': row[0],
                            'start_pos': TES_start_pos,
                            'end_pos': TES_end_pos,
                            'x': X,
                            'y': Y,
                            'gene_name': gene_name,
                            'rna_name': rna_name,
                            'orig_start': row[1],
                            'orig_end': row[2],
                            'strand': row[5],
                            'reg':'TES'})
                            #'names': make_names(row, TES_start_pos, TES_end_pos, X, Y,
                            #                    'TES')})
 
                #added wang
                if TES_start_pos<1:
                   TES_start_pos=1

                #5'distance upstream
                dist5_start_pos = row[1]-N if row[5]=='+' else row[2]+M
                dist5_end_pos = row[1]-M if row[5]=='+' else row[2]+N
                if dist5_start_pos < 1:
                    dist5_start_pos = 1
                #added wang
                if dist5_end_pos <1:
                    dist5_end_pos =1

                if not (dist5_start_pos==1 & dist5_end_pos==1):
                   dist5.append({'chr': row[0],
                            'start_pos': dist5_start_pos,
                            'end_pos': dist5_end_pos,
                            'x': M,
                            'y': N,
                            'row': row,
                            'gene_name': gene_name,
                            'rna_name': rna_name,
                            'orig_start': row[1],
                            'orig_end': row[2],
                            'strand': row[5],
                            'reg':'5dist'})
                            #'names': make_names(row, dist5_start_pos, dist5_end_pos, M,
                            #                    N, '5dist')})

                #added jbw 5distance downstream
                dist5D_start_pos= row[1]+M  if row[5]=='+' else row[2] -N
                dist5D_end_pos = row[1] +N  if row[5]=='+' else row[2] -M
                if dist5D_start_pos<1:
                   dist5D_start_pos=1

                if dist5D_end_pos<1:
                   dist5D_end_pos=1

                if not (dist5D_start_pos==1 & dist5D_end_pos==1):
                   dist5D.append({'chr': row[0],
                            'start_pos': dist5D_start_pos,
                            'end_pos': dist5D_end_pos,
                            'x': M,
                            'y': N,
                            'row': row,
                            'gene_name': gene_name,
                            'rna_name': rna_name,
                            'orig_start': row[1],
                            'orig_end': row[2],
                            'strand': row[5],
                            'reg':'5distD'})
 

    TSS_df = pd.DataFrame(TSS)
    TSS_df = combine_non_unique(TSS_df)

    TES_df = pd.DataFrame(TES)
    TES_df = combine_non_unique(TES_df)

    gene_df = pd.DataFrame(gene)
    gene_df = combine_non_unique(gene_df)

    dist5_df = pd.DataFrame(dist5)
    dist5_df = combine_non_unique(dist5_df)

    removed_df = pd.DataFrame(removed)
    if not removed_df.empty:
       removed_df = combine_non_unique(removed_df)
   
    #added jbw
    dist5D_df=pd.DataFrame(dist5D)
    dist5D_df=combine_non_unique(dist5D_df)
  
    return TSS_df, TES_df, gene_df, dist5_df, removed_df, dist5D_df








def make_intergenic_region(genome_file, reference, bed_files, min_intergenic_reg_len, max_intergenic_reg_len, X, Y,
                           out_folder, intergenic_between_genes):
    """Finds intergenic regions with regions bigger than min_intergenic_reg_len,
    and returns dataframe.

    genome_file -- genome file name
    bed_files -- BED filenames: [TSS_filename, TES_filename, gene_filename]
    min_intergenic_reg_len -- minimum size of intergenic region
    out_folder -- folder name where output should go
    intergenic_between_genes -- complement on gene only or TSS, TES and gene? find intergenic
                  regions between gene regions, or include TSS and TES as well?"""

    if intergenic_between_genes:
        # Runs the command line which creates intergenic regions BED file
        subprocess.run('bedtools complement -i '+reference+' -g '
                  +genome_file+' > '+out_folder+'/'+'intergenic_regions.bed',shell=True)
    else:
        subprocess.run('cat '+bed_files[0]+' '+ bed_files[1]+' '+bed_files[2]+\
                  ' > '+out_folder+'/'+'combined.bed',shell=True)

        #added jbw to check whether negative values in bed file and replace them as 1
        tmp_file=out_folder+'/'+'combined.bed'
        tmp_df=pd.read_csv(tmp_file,sep='\t',header=None)
        tmp_df.loc[tmp_df[1].apply(lambda x: x<0),1]=1
        tmp_df.to_csv(tmp_file,sep='\t',index=False, header=None)
        #end change jbw

        subprocess.run('bedtools sort -g '+genome_file+' -i '+out_folder+'/'+\
                  'combined.bed > '+out_folder+'/'+'combined_sorted.bed',shell=True)

        subprocess.run('bedtools merge -i '+out_folder+'/'+'combined_sorted.bed > \
                  '+out_folder+'/'+'combined_sorted_merged.bed',shell=True)

        subprocess.run('bedtools complement -i '+out_folder+'/'+'combined_sorted_merged.bed -g ' +
                  genome_file + ' > '+out_folder+'/'+'intergenic_regions.bed',shell=True)
        #jbw test
        subprocess.run('rm '+out_folder+'/'+'combined*.bed', shell=True)

    inter_gen = pd.read_csv(out_folder+'/'+'intergenic_regions.bed', sep='\t',
                names=['chr', 'start_pos', 'end_pos'])
    #jbw test
    subprocess.run('rm '+out_folder+'/'+'intergenic_regions.bed', shell=True)

    # Removes intergenetic regions shorter than parameter
    inter_gen = inter_gen[((inter_gen['end_pos'] - inter_gen['start_pos']) >
                min_intergenic_reg_len) & ((inter_gen['end_pos'] - inter_gen['start_pos']) <
                max_intergenic_reg_len)]

    inter_gen['names'] = inter_gen.chr.astype(str) + ':' +\
                         inter_gen.start_pos.astype(str) + ':' +\
                         inter_gen.end_pos.astype(str) + '||' +\
                         'intergenic:' +\
                         ('None' if intergenic_between_genes else str(X)) + ':' +\
                         ('None' if intergenic_between_genes else str(Y))
    return inter_gen





def make_region_files(reference, genomeFile, X, Y, M, N, intergenic_between_genes, min_intergenic_len, max_intergenic_len,
         remove_short, out_folder):

    
    
    #test jbw 06.23
    #input refFlat file which has to be cleaned before using it !!
    ref_df = pd.read_csv(reference, index_col=False, header=0, sep='\t')
    #df = pd.read_csv(reference, sep='\t',header=0, usecols=[0, 1, 2, 3, 4, 5],
    #            names=['gene_name', 'rna_name', 'chr', 'strand', 'TSS', 'TES'])
    #new_ref = find_uniq_genes(df, numerical_chr_name, human, remove_mir)
 

    reference_name = reference.split('/')[-1][:-4]

    dfs = make_TSS_TES_gene_5dist(ref_df, X, Y, M, N, remove_short)
    #added jbw
    TSS_df, TES_df, gene_df, dist5_df, removed_df , dist5D_df= dfs


    if remove_short:
        TSS_filename = out_folder+'/TSS_Up'+str(X)+'_Down'+str(Y)+'_removedShort'+'.bed'
        TES_filename = out_folder+'/TES_Up'+str(X)+'_Down'+str(Y)+'removedShort'+'.bed'
        gene_filename = out_folder+'/gene_Up'+str(X)+'_Down'+str(Y)+'removedShort'+'.bed'
        dist5_filename = out_folder+'/5dist_Up'+str(N)+'_Up'+str(M)+'removedShort'+'.bed'
        removed_filename = out_folder+'/removed_regions_all_TSS_TES_5dist_geneBodyLessThan0.bed'
        #added jbw
        dist5D_filename= out_folder+ '/5dist_Down'+str(N)+'_Down'+str(M)+'removedShort'+'.bed'



    order = ['chr', 'start_pos', 'end_pos', 'names']
    TSS_df = TSS_df[order]
    TES_df = TES_df[order]
    gene_df = gene_df[order]
    dist5_df = dist5_df[order]
    if not removed_df.empty :
       removed_df = removed_df[order]
    #added jbw
    dist5D_df= dist5D_df[order]

    #added wang check negative postion in exported files
    TSS_df.loc[TSS_df['start_pos'].apply(lambda x: x<0),'start_pos']=1
    TES_df.loc[TES_df['start_pos'].apply(lambda x: x<0),'start_pos']=1

    TSS_df.to_csv(TSS_filename, sep='\t', index=False, header=None)
    TES_df.to_csv(TES_filename, sep='\t', index=False, header=None)
    gene_df.to_csv(gene_filename, sep='\t', index=False, header=None)
    dist5_df.to_csv(dist5_filename, sep='\t', index=False, header=None)
    if not removed_df.empty:
       removed_df.to_csv(removed_filename, sep='\t', index=False, header=None)
    #added jbw
    dist5D_df.to_csv(dist5D_filename, sep='\t', index=False, header=None)

    # Creates filtered intergenic region BED file
    bed_files = [TSS_filename, TES_filename, gene_filename]
    inter_gen = make_intergenic_region(genomeFile, reference, bed_files, min_intergenic_len, max_intergenic_len,
                                       X, Y, out_folder, intergenic_between_genes)
    inter_gen = inter_gen[order]

    if intergenic_between_genes:
        intergen_filename = out_folder+'/'+'intergenic_uniqueSorted_betweenGenes_minLen'+str(min_intergenic_len)+'.bed'
    else:
        intergen_filename = out_folder+'/'+'intergenic_uniqueSorted_betweenTSS_TES_genes_minLen'+str(min_intergenic_len) + '.bed'
    inter_gen.to_csv(intergen_filename, sep='\t', index=False, header=None)
    
    #added jbw
    return [TSS_filename, TES_filename, gene_filename, dist5_filename, intergen_filename, removed_filename, dist5D_filename]





def make_reference(reference,out_folder):
    df = pd.read_csv(reference, sep='\t', usecols=[0, 1, 2, 3, 4, 5],
                names=['gene_name', 'rna_name', 'chr', 'strand', 'TSS', 'TES'])
    
    outfilename = reference.split('/')[-1][:-4]
    
    new_ref = find_uniq_genes(df, out_folder)
    
    new_ref.rename(columns={'TSS': 'start_pos', 'TES': 'end_pos'}, inplace=True)
    
    ref_BED_df = make_BED_of_ref(new_ref)

    # Makes sure the columns are in this exact order
    new_ref = new_ref[['rna_name', 'chr', 'strand', 'start_pos', 'end_pos', 'gene_name']]
    ref_BED_df = ref_BED_df[['chr', 'start_pos', 'end_pos', 'gene_rna_name', 'dot', 'strand']]

    # Keeps header
    new_ref.to_csv(out_folder+'/'+outfilename+'_clean_sorted.txt', sep='\t', index=False)

    out_bed = out_folder+'/'+outfilename+'_clean_sorted.bed'
    # No header
    ref_BED_df.to_csv(out_bed, sep='\t', index=False, header=None)
    
    return ref_BED_df, out_bed



def main(reference, genomeFile, X, Y, M, N, intergenic_between_genes, min_intergenic_len, max_intergenic_len,
         remove_short, out_folder, enchancer):
    
    remove_short = True
    intergenic_between_genes = True
    
    ref_BED_df,reference_bed = make_reference(reference, out_folder)
    ref_BED_df.columns = ['chr', 'start_pos', 'end_pos', 'gene_name:rna_name',
                                     'dot', 'strand']
    
    dfs = make_TSS_TES_gene_5dist(ref_BED_df, X, Y, M, N, remove_short)
    
   #added jbw
    TSS_df, TES_df, gene_df, dist5_df, removed_df , dist5D_df= dfs
   
    TSS_filename = out_folder+'/TSS_Up'+str(X)+'_Down'+str(Y)+'_removedShort'+'.bed'
    TES_filename = out_folder+'/TES_Up'+str(X)+'_Down'+str(Y)+'removedShort'+'.bed'
    gene_filename = out_folder+'/gene_Up'+str(X)+'_Down'+str(Y)+'removedShort'+'.bed'
    dist5_filename = out_folder+'/5dist_Up'+str(N)+'_Up'+str(M)+'removedShort'+'.bed'
    removed_filename = out_folder+'/removed_regions_all_TSS_TES_5dist_geneBodyLessThan0.bed'

    dist5D_filename= out_folder+ '/5dist_Down'+str(N)+'_Down'+str(M)+'removedShort'+'.bed'
    
    order = ['chr', 'start_pos', 'end_pos', 'names']
    TSS_df = TSS_df[order]
    TES_df = TES_df[order]
    gene_df = gene_df[order]
    dist5_df = dist5_df[order]
    if not removed_df.empty :
       removed_df = removed_df[order]
    #added jbw
    dist5D_df= dist5D_df[order]
    
    TSS_df.loc[TSS_df['start_pos'].apply(lambda x: x<0),'start_pos']=1
    TES_df.loc[TES_df['start_pos'].apply(lambda x: x<0),'start_pos']=1

    TSS_df.to_csv(TSS_filename, sep='\t', index=False, header=None)
    TES_df.to_csv(TES_filename, sep='\t', index=False, header=None)
    gene_df.to_csv(gene_filename, sep='\t', index=False, header=None)
    dist5_df.to_csv(dist5_filename, sep='\t', index=False, header=None)
    if not removed_df.empty:
      removed_df.to_csv(removed_filename, sep='\t', index=False, header=None)

    dist5D_df.to_csv(dist5D_filename, sep='\t', index=False, header=None)
    
    
    
    # Creates filtered intergenic region BED file
    bed_files = [TSS_filename, TES_filename, gene_filename]
    inter_gen = make_intergenic_region(genomeFile, reference_bed, bed_files, min_intergenic_len, max_intergenic_len,
                                       X, Y, out_folder, intergenic_between_genes)
    inter_gen = inter_gen[order]

    if intergenic_between_genes:
        intergen_filename = out_folder+'/'+'intergenic_uniqueSorted_betweenGenes_minLen'+str(min_intergenic_len)+'.bed'
    else:
        intergen_filename = out_folder+'/'+'intergenic_uniqueSorted_betweenTSS_TES_genes_minLen'+str(min_intergenic_len) + '.bed'
    inter_gen.to_csv(intergen_filename, sep='\t', index=False, header=None)
    
   
    
    out_list_file = out_folder+'/list_region_files.txt'
    
    outfiles = open(out_list_file, 'w')
    outfiles.write(TSS_filename+'\n')
    outfiles.write(gene_filename+'\n')
    outfiles.write(TES_filename+'\n')
    outfiles.write(dist5_filename+'\n')
    outfiles.write(intergen_filename+'\n')
    #added jbw
    outfiles.write(dist5D_filename +'\n')
    if isinstance(enchancer, str):
        outfiles.write(enchancer +'\n')
    outfiles.close()
    
    return out_list_file
    
    
    
    
