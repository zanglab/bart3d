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
from scipy import stats
from bart3d import local_interaction_from_hicpro, local_interaction_from_hic, local_interaction_from_cool


def normalized_mirror_concat(df):
    # ==== normalize the contact matrix by dividing averaged reads of bins with the same distance 
    df = df/df.mean()  # median =0 in most cases
    init_cols = df.columns
    # ==== get up-stream chromatin interactions
    for column_pos in np.arange(1,len(init_cols)):
        column = init_cols[column_pos]
        mir_column = '-{}'.format(column)
        new_values = np.append([0]*column_pos,df[column].values)[:len(df.index)]
        df.insert(0,mir_column,new_values)
    return df


        
def differential_score(args):
    
    species=args.species
    genomic_distance = args.genomicDistance
    tmpdir = args.outdir
    os.makedirs(tmpdir,exist_ok=True)
    
    # choose the right chrom and chrom_len
    chrom_len = pd.read_csv(args.chromsize,header=None,index_col=0,sep='\t')[1].to_dict()
    chroms = list(chrom_len.keys())
    
    # == read data from each of the treatment files
    treat_matrix_files = args.treatment.split(',')
    treat_matrix_data = {}
    for treat_i in np.arange(len(treat_matrix_files)):
        treat_matrix_file = treat_matrix_files[treat_i]
        # ==== read local chromatin interaction 
        if args.fileFormat=='hicpro':
            if not args.bedFileHicpro:
                sys.stderr.write('Please provide the HiC-Pro abs/ord .bed file through --bedFileHicpro! \n')
                sys.exit(1)
            # ==== read the ord/matrix files    
            print('== Reading treatment contact matrix: {}'.format(treat_i+1))
            treat_ord_file = args.bedFileHicpro
            treat_interaction_dfs,treat_resolution = local_interaction_from_hicpro.get_local_interaction_from_matrix(treat_ord_file,treat_matrix_file,tmpdir,chrom_len,genomic_distance)
            
        elif args.fileFormat=='hic':
            print('== Reading treatment contact matrix: {}'.format(treat_i+1))
            treat_interaction_dfs,treat_resolution = local_interaction_from_hic.get_local_interaction_from_matrix(args.treatment,args.resolution,tmpdir,chrom_len,genomic_distance)
        
        elif args.fileFormat=='cool':
            print('== Reading treatment contact matrix: {}'.format(treat_i+1))
            treat_interaction_dfs,treat_resolution = local_interaction_from_cool.get_local_interaction_from_matrix(args.treatment,tmpdir,chrom_len,genomic_distance)
        
        # save the treatment contact matrix
        treat_matrix_data[treat_i] = treat_interaction_dfs
 
    # == read data from each of the control files
    control_matrix_files = args.control.split(',')
    control_matrix_data = {}     
    for control_j in np.arange(len(control_matrix_files)):
        control_matrix_file = control_matrix_files[control_j]
        # ==== read local chromatin interaction 
        if args.fileFormat=='hicpro':
            print('== Reading control contact matrix: {}'.format(control_j+1))
            control_ord_file = args.bedFileHicpro
            control_interaction_dfs,control_resolution = local_interaction_from_hicpro.get_local_interaction_from_matrix(control_ord_file,control_matrix_file,tmpdir,chrom_len,genomic_distance)
            
        elif args.fileFormat=='hic':
            print('== Reading control contact matrix: {}'.format(control_j+1))
            control_interaction_dfs,control_resolution = local_interaction_from_hic.get_local_interaction_from_matrix(args.control,args.resolution,tmpdir,chrom_len,genomic_distance)
        
        elif args.fileFormat=='cool':
            print('== Reading control contact matrix: {}'.format(control_j+1))
            control_interaction_dfs,control_resolution = local_interaction_from_cool.get_local_interaction_from_matrix(args.control,tmpdir,chrom_len,genomic_distance)
        
        # save the control contact matrix
        control_matrix_data[control_j] = control_interaction_dfs    
 

    # ==== treat and control file should with the same resolution
    try:
        assert treat_resolution == control_resolution
    except AssertionError:
        print('Input matrices are generated at different resolutions.')
        exit(1)


    # ==== write out genome wide differential scores    
    check_bins = int(genomic_distance/treat_resolution)
    diff_dict={}
    # ==== paired t-test between each treatment and each control
    for treat_i in np.arange(len(treat_matrix_files)):
        for control_j in np.arange(len(control_matrix_files)):
            loop_key='treatment{}_vs_control{}'.format(treat_i,control_j)
            diff_dict[loop_key] ={}
            
            treat_interaction_dfs = treat_matrix_data[treat_i]
            control_interaction_dfs = control_matrix_data[control_j]  
            # print('paired ttest, {} {}'.format(treat_i,control_j))
            for chr in chroms:
                if (chr in treat_interaction_dfs.keys()) and (chr in control_interaction_dfs.keys()):
                    treat_interaction = treat_interaction_dfs[chr].astype(float)
                    control_interaction = control_interaction_dfs[chr].astype(float)
                    # ==== in case of un-matched index at the end of each chromosome
                    shared_index = [i for i in treat_interaction.index if i in control_interaction.index]
                    treat_interaction = treat_interaction.loc[shared_index]
                    control_interaction = control_interaction.loc[shared_index]
                    # ==== get up-stream interactions
                    treat_interaction = normalized_mirror_concat(treat_interaction)
                    control_interaction = normalized_mirror_concat(control_interaction)
                    # ==== differential score,fill nan statistics with 0
                    s,p = stats.ttest_rel(treat_interaction,control_interaction,axis=1,nan_policy='omit')
                    #s[np.isnan(s)]=0
            
                    # ==== bins with #left/right interaction < check_bins
                    for ii in np.arange(0,treat_interaction.shape[0]):
                        # ==== start position of current bin
                        bin_start_pos = treat_interaction.index[ii]#;print('b',ii)
                        # ==== remove added 0 values 
                        if ii<check_bins or treat_interaction.shape[0]-1-ii<check_bins:
                            local_interaction_left_bin = max(check_bins-1-ii,0)
                            local_interaction_right_bin = min(treat_interaction.shape[1], check_bins+treat_interaction.shape[0]-1-int(ii))
                            # print('====',local_interaction_left_bin,local_interaction_right_bin)
                            treat_val = treat_interaction.iloc[ii][local_interaction_left_bin:local_interaction_right_bin]
                            control_val = control_interaction.iloc[ii][local_interaction_left_bin:local_interaction_right_bin]
                            s_current,p_current = stats.ttest_rel(treat_val,control_val)
                            # s_current=0 if np.isnan(s_current) else s_current 
                            statistic_score = -2*np.log(p_current)*np.sign(s_current)
                        # ==== otherwise use values in array(s)
                        else:
                            # s_current = s[ii]
                            statistic_score = -2*np.log(p[ii])*np.sign(s[ii])
                            # print('====',treat_interaction.shape[1])
                        # ==== check if is the last bin
                        if ii < treat_interaction.shape[0]-1:
                            diff_index='{}_{}_{}'.format(chr,bin_start_pos,bin_start_pos+treat_resolution)
                        else:
                            diff_index='{}_{}_{}'.format(chr,bin_start_pos,chrom_len[chr])
                        # ==== save the quantified score
                        diff_dict[loop_key][diff_index] = statistic_score
    
    # keep the abs max to quantify if the treatment is increased or decreased than control
    diff_df = pd.DataFrame.from_dict(diff_dict)
    diff_df_cp = diff_df.copy()
    
    # p-value of fisher's combined score for increased interactions
    diff_df['pos_fisher_score'] = diff_df_cp[diff_df_cp>0].sum(axis=1)
    diff_df['pos_fisher_freedom'] = 2*diff_df_cp[diff_df_cp>0].notnull().sum(axis=1)
    diff_df['pos_fisher_pvalue'] = stats.chi2.sf(diff_df['pos_fisher_score'],diff_df['pos_fisher_freedom'])
    diff_df['pos_fisher_neglog_pvalue'] = -1*np.log(diff_df['pos_fisher_pvalue']).fillna(0)

    # p-value of fisher's combined score for decreased interactions
    diff_df['neg_fisher_score'] = -1*diff_df_cp[diff_df_cp<0].sum(axis=1)
    diff_df['neg_fisher_freedom'] = 2*diff_df_cp[diff_df_cp<0].notnull().sum(axis=1)
    diff_df['neg_fisher_pvalue'] = stats.chi2.sf(diff_df['neg_fisher_score'],diff_df['neg_fisher_freedom'])
    diff_df['neg_fisher_neglog_pvalue'] = -1*np.log(diff_df['neg_fisher_pvalue']).fillna(0)
    
    # keep the neglog fisher pvalue with more degrees of freedom (more supporter), if equal, keep the bigger one
    flag_abs = diff_df['pos_fisher_neglog_pvalue']-diff_df['neg_fisher_neglog_pvalue']
    diff_df['fisher_neglog_pvalue_abs_max']=(flag_abs>=0)*diff_df['pos_fisher_neglog_pvalue']+(flag_abs<0)*(-1)*diff_df['neg_fisher_neglog_pvalue']
    # compare the degrees of freedom
    flag = np.sign(diff_df['pos_fisher_freedom']-diff_df['neg_fisher_freedom']) 
    diff_df['final_score']=(flag==1)*diff_df['pos_fisher_neglog_pvalue']+(flag==-1)*(-1)*diff_df['neg_fisher_neglog_pvalue']+(flag==0)*diff_df['fisher_neglog_pvalue_abs_max']
    # diff_df.to_csv('{}/{}_paired_ttest_results.csv'.format(args.outdir,args.outFileName))
    
    # ==== write the bed format DCI profile
    diff_file_name = '{}/{}_differential_score.bed'.format(args.outdir,args.outFileName)
    diff_file = open(diff_file_name,'w')
    for ii in diff_df.index:
        chr,start,end = ii.split('_')   
        score = diff_df.at[ii,'final_score']     
        diff_file.write('{}\t{}\t{}\t{}\n'.format(chr,start,end,score))
    diff_file.close()
    
    # ==== generate the bigwig file    
    try:
        os.system('LC_COLLATE=C sort -k1,1 -k2,2n {}> {}.sorted'.format(diff_file_name,diff_file_name))
        os.system('{} {}.sorted {} {}.sorted.bw'.format(args.BgToBigWig,diff_file_name,args.chromsize,diff_file_name))
        os.system('rm {}.sorted'.format(diff_file_name))
    except:
        print('{} {}.sorted {} {}.sorted.bw'.format(args.BgToBigWig,diff_file_name,args.chromsize,diff_file_name))
        print('Cannot generate the bigwig file. Please make sure the bedGraphToBigWig file is executable.')
        print('It may help to try: \n $ chmod 755 {}'.format(args.BgToBigWig))
    
    return diff_file_name,treat_matrix_data,control_matrix_data
    





# def main(args):
# 
#     # ==== test demo
#     args.species='hg38'
#     matrix_format='hicpro'
#     args.genomic_distance = 200000
#     args.ofilename=None
#     args.chromsize='/nv/vol190/zanglab/zw5j/data/Genome/ucsc/{}/{}_clean.chrom.sizes'.format(args.species,args.species)
#     
#     indir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/11_HiC_analysis/f1_preprocess/HiC_Pro/hicPro_raw_fq/panos_hic/hic_results/matrix/'
#     args.treat1 = indir+os.sep+'Jurkat/raw/5000/Jurkat_5000_ord.bed'
#     args.treat2= indir+os.sep+'Jurkat/raw/5000/Jurkat_5000.matrix'
#     args.control1 = indir+os.sep+'A6010/raw/5000/A6010_5000_ord.bed'
#     args.control2 = indir+os.sep+'A6010/raw/5000/A6010_5000.matrix'
#    
#     # tmp dir to save local interaction for each chrom
#     args.outdir = 'tmp_'+''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(20))
#     differential_score(args)
#     
# 
# if __name__ == '__main__':
# 
#     parser = argparse.ArgumentParser()
# #     parser.add_argument('-v', '--genomic_distance', action = 'store', type = int,dest = 'genomic_distance', help = 'input file of', metavar = '<int>')
# #     parser.add_argument('-n', '--normalization', action = 'store', type = str,dest = 'normalization', help = 'input file of', metavar = '<str>')
# #     parser.add_argument('-c', '--chrom', action = 'store', type = str,dest = 'chrom', help = 'input file of', metavar = '<str>')
#     #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
#     #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
#     #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
#     #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
#     
#     args = parser.parse_args()
#     if(len(sys.argv))<0:
#         parser.print_help()
#         sys.exit(1)
#   
#     main(args)
