'''
Copyright (c) 2020, 2021 Zhenjia Wang <zhenjia@virginia.edu>, Chongzhi Zang <zang@virginia.edu>

This code is free software; you can redistribute it and/or modify it 
under the terms of the BSD License.

@author: Zhenjia Wang, Chongzhi Zang
@contact: zhenjia@virginia.edu, zang@virginia.edu
'''

import sys,argparse
import os,glob,shutil
import numpy as np
import pandas as pd
import random,string
import scipy
from scipy import stats



def linear_regression(x,y):
    xmean = np.mean(x)
    ymean = np.mean(y)
    assert len(x)==len(y)
    sum1,sum2=0,0
    x_values,y_values = x.values,y.values
    for i in np.arange(len(x)):
        sum1 += (x_values[i]-xmean)*(y_values[i]-ymean)
        sum2 += np.power((x_values[i]-xmean),2)
    a = sum1/sum2
    b = ymean-a*xmean
    return a,b,xmean,ymean


def return_residual_and_std(df,x,y,a,b):
    # residual of each region/gene
    y_reg = a*x+b
    y_residual = y-y_reg
    miu = np.mean(y_residual)
    std = np.std(y_residual)
    # p-value of each residual
    df['residual']=y_residual
    df['residual_std']=std
    df['residual_pvalue'] = scipy.stats.norm(miu,std).sf(y_residual)
#     df['residual_pvalue'] = pd.concat([df['residual_pvalue'],(1-df['residual_pvalue'])],axis=1).min(axis=1)
    return df


def read_avg_coverage_from_matrices(matrices):
    # read the average of normalized coverage from the matrices
    df = pd.DataFrame()
    for matrix_ii in np.arange(len(matrices)):
        # for each matrix, read the local coverages (diagonal values)
        matrix_data = matrices[matrix_ii]
        df_coverage = pd.DataFrame()
        coverage_col='coverage_{}'.format(matrix_ii)
        for chrom in matrix_data.keys():
            df_coverage_chr = matrix_data[chrom][[0]] +1 ## pseudo-count =1
            df_coverage_chr.columns=[coverage_col] 
            df_coverage_chr['chr']=chrom
            df_coverage_chr['pos']=df_coverage_chr.index
            df_coverage = pd.concat([df_coverage,df_coverage_chr])
        ## normalize the reads in each sample
        # df_coverage = df_coverage[df_coverage[coverage_col]>0]
        avg_reads = df_coverage[coverage_col].mean()
        df_coverage[coverage_col]=df_coverage[coverage_col]/avg_reads
        if df.shape[0]==0:
            df = df_coverage.copy()
        else:
            df = pd.merge(df,df_coverage,how='inner',left_on=['chr','pos'],right_on=['chr','pos'])      
    # get the averaged coverage
    df['avg_coverage'] = df[df.columns.difference(['chr','pos'])].mean(axis=1)
    # df = df[df['avg_coverage']>0]
    return df



def profile_normalization(args,treat_matrix_data,control_matrix_data):

    dci = pd.read_csv(args.diff_file,sep='\t',header=None)
    dci.columns=['chr','start','end','DCI']
    
    # ==== read averaged coverage from matrix
    treat_matrix_df = read_avg_coverage_from_matrices(treat_matrix_data)
    control_matrix_df = read_avg_coverage_from_matrices(control_matrix_data)
    coverage_df = pd.merge(treat_matrix_df,control_matrix_df,how='inner',left_on=['chr','pos'],right_on=['chr','pos'])
    coverage_df['avg_coverage_log2FC'] = np.log2(coverage_df['avg_coverage_x']/(coverage_df['avg_coverage_y']))

    # ==== merge coverage log2FC and DCI
    compr_df = pd.merge(coverage_df,dci,how='inner',left_on=['chr','pos'],right_on=['chr','start'])
#     compr_df.to_csv(outdir+os.sep+'{}_coverage_and_DCI.csv'.format(outpur_prename))

    # ==== linear regression
    x = compr_df['avg_coverage_log2FC']
    y = compr_df['DCI']
    a,b,xmean,ymean = linear_regression(x,y)
    # ==== residual
    compr_df = return_residual_and_std(compr_df,x,y,a,b)      
    normalized_dci = compr_df[['chr','start','end','residual']]
    # ==== PCC 
    s,p = stats.pearsonr(x,y)
    print('Pearson correlation coefficient between DCI and coverage changes BEFORE normalization: {:.4f}'.format(s))
    y = compr_df['residual']
    s,p = stats.pearsonr(x,y)
    print('Pearson correlation coefficient between DCI and coverage changes AFTER normalization: {:.4f}'.format(s))
    
    # ==== save the normalized DCI file
    diff_file_name = '{}/{}_differential_score_after_coverageNormalization.bed'.format(args.outdir,args.outFileName)
    normalized_dci.to_csv(diff_file_name,sep='\t',index=False,header=None)  
       
    # ==== generate the bigwig file    
    try:
        os.system('LC_COLLATE=C sort -k1,1 -k2,2n {}> {}.sorted'.format(diff_file_name,diff_file_name))
        os.system('{} {}.sorted {} {}.sorted.bw'.format(args.BgToBigWig,diff_file_name,args.chromsize,diff_file_name))
        os.system('rm {}.sorted'.format(diff_file_name))
    except:
        print('{} {}.sorted {} {}.sorted.bw'.format(args.BgToBigWig,diff_file_name,args.chromsize,diff_file_name))
        print('Cannot generate the bigwig file. Please make sure the bedGraphToBigWig file is executable.')
        print('It may help to try: \n $ chmod 755 {}'.format(args.BgToBigWig))
    
    return diff_file_name
    




