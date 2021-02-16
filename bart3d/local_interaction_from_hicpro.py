'''
Copyright (c) 2020, 2021 Zhenjia Wang <zhenjia@virginia.edu>, Chongzhi Zang <zang@virginia.edu>

This code is free software; you can redistribute it and/or modify it 
under the terms of the BSD License.

@author: Zhenjia Wang, Chongzhi Zang
@contact: zhenjia@virginia.edu, zang@virginia.edu
'''

import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect
import random,time,string


def get_lines(infile):
    with open(infile,'rb') as f:
        lines = 0
        buf_size = 1024*1024
        buf = f.raw.read(buf_size)
        while buf:
            lines += buf.count(b'\n')
            buf = f.raw.read(buf_size)
    return lines


    
def get_local_interaction_from_matrix(ordfile,matrixfile,tmpdir,chrom_len,genomic_distance):
     
    # ==== for each chr, get the first/end ID in ord file 
    # ==== read out the genomic_distance limited interactions for each bin
    ord_df = pd.read_csv(ordfile,index_col=0,header=None,sep='\t')
    ord_df.columns = ['start','end','id']
    chroms = ord_df.index.drop_duplicates()
    # chroms = [i for i in chroms if i in chrom_len.keys()]
    resolution = ord_df.iloc[0,1] - ord_df.iloc[0,0]#;print(resolution)
    check_bins = int(genomic_distance/resolution)


    # ==== local chromatin interaction
    chr_ord_id = {} # ids of all bins in each chrom
    chr_local_interaction = {} # local interaction score for each chr
    for chr in chroms:
        # ==== first/end ID in ord file
        ord_df_tmp = ord_df[ord_df.index==chr]
        # print(chr,ord_df_tmp.shape)
        chr_start_ID = ord_df_tmp.iloc[0,-1]
        chr_start_pos = ord_df_tmp.iloc[0,0]
        chr_end_ID = ord_df_tmp.iloc[-1,-1]
        chr_end_pos = ord_df_tmp.iloc[-1,0]
        # == check genome assembly
        if chr in chrom_len.keys():
            try:
                assert chr_end_pos <= chrom_len[chr]
            except AssertionError:
                print('End coordinate {} bigger than {} size {}!'.format(chr_end_pos,chr,chrom_len[chr]))
                print('Make sure you are using the matched genome assembly to process the HiC data!')
                exit(1)
        # initiate the local interaction dict
        chr_ord_id[chr] = [chr_start_ID,chr_start_pos,chr_end_ID,chr_end_pos]
        chr_local_interaction[chr] = {}
        for bin in np.arange(chr_start_ID,chr_end_ID+1):
            chr_local_interaction[chr][int(bin)] = [0]*check_bins


    # ==== all start ID of all chroms
    all_chr_starts = [chr_ord_id[i][0] for i in chroms]
    assert sorted(all_chr_starts)==all_chr_starts

    # ==== read the matrix, and keep only those bins with pairwise dis < genomic_distance
    total_lines = get_lines(matrixfile)
    chunksize=2e7
    total_chunks = int(np.ceil(total_lines/chunksize))
    current_line=0;chunk=1
    
    matrix_inf = open(matrixfile) 
    line = matrix_inf.readline()
    while line:        
        # time1 = time.time()
        sline = line.strip().split()
        leftID = int(sline[0])
        rightID = int(sline[1])
        # ==== the inter-chr interactions might be included but will be removed later for differential score
        if 0<=rightID-leftID <check_bins:
            score = float(sline[2])
            current_chr = chroms[bisect.bisect_right(all_chr_starts,leftID)-1]#;time2 = time.time();print('a',time2-time1)
            # print(leftID,current_chr)
            chr_local_interaction[current_chr][leftID][rightID-leftID] = score#;time3 = time.time();print('b',time3-time2)
        line = matrix_inf.readline()
        # ==== keep record of reading
        if current_line%chunksize==0:
            print('== processing {}/{}'.format(chunk,total_chunks))
            chunk+=1
        current_line+=1
    matrix_inf.close()

    
    # ==== write out the local interaction, for double check
#     flag = os.path.basename(matrixfile).split('.matrix')[0]
#     flag = '{}_{}'.format(flag,''.join(random.choice(string.ascii_letters + string.digits) for _ in range(6)))
#     for chr in chroms:
#         if chr in chrom_len.keys():
#             chr_start_id = chr_ord_id[chr][0]
#             chr_start_pos = chr_ord_id[chr][1]
#             outfile = tmpdir+os.sep+'tmp_{}_{}_res_{}_view_region_{}.csv'.format(flag,chr,resolution,genomic_distance)
#             outf = open(outfile,'w')
#             outf.write('{}\t{}\n'.format('dis','\t'.join(map(str,np.arange(0,genomic_distance,resolution)))))
#             for view_id in chr_local_interaction[chr].keys():
#                 view_pos = chr_start_pos+(view_id-chr_start_id)*resolution
#                 outf.write('{}\t{}\n'.format(view_pos,'\t'.join(map(str,chr_local_interaction[chr][view_id]))));
#             outf.close()
    
    # save the local interaction
    interaction_dfs = {} # interaction dataframe
    for chr in chroms:
        if chr in chrom_len.keys():
            chr_start_id = chr_ord_id[chr][0]
            chr_start_pos = chr_ord_id[chr][1]
            interaction_dfs[chr] = pd.DataFrame.from_dict(chr_local_interaction[chr],orient='index')
            columns = np.arange(0,genomic_distance,resolution)
            index = [chr_start_pos+(view_id-chr_start_id)*resolution for view_id in interaction_dfs[chr].index]
            interaction_dfs[chr].columns = columns
            interaction_dfs[chr].index = index
            
            # ==== save the interaction matrix
            flag = os.path.basename(matrixfile).split('.matrix')[0]
            # interaction_dfs[chr].to_csv(tmpdir+os.sep+'tmp_{}_{}.csv'.format(flag,chr))
                
    return interaction_dfs,resolution



def main(ordfile,matrixfile,tmpdir,species,genomic_distance):
    
    # test demo
    chrom_len = GenomeData.species_chrom_lengths['hg38']

    indir='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/11_HiC_analysis/f1_preprocess/HiC_Pro/hicPro_raw_fq/panos_hic/hic_results/matrix/'
    ordfile=indir+'/Jurkat/raw/5000/Jurkat_5000_ord.bed'
    matrixfile=indir+'/Jurkat/raw/5000/Jurkat_5000.matrix'
    
    indir='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/patient_hic_compr_GSE115896/f0_process_panos/panos_patient_out/hic_results/matrix/PD9/raw/50000'
    ordfile=indir+'/PD9_50000_ord.bed'
    matrixfile=indir+'/PD9_50000.matrix'
    
    os.makedirs(args.outdir,exist_ok=True)
    
    # for input ord/matrix file, get the local interaction for each chrom  
    flag,resolution = get_local_interaction_from_matrix(ordfile,matrixfile,tmpdir,chrom_len,genomic_distance)




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--infile1', action = 'store', type = str,dest = 'infile1', help = 'input ord.bed file', metavar = '<dir>')
    parser.add_argument('-m', '--infile2', action = 'store', type = str,dest = 'infile2', help = 'input matrix file', metavar = '<dir>')
    parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    parser.add_argument('-v', '--genomicDistance', action = 'store', type = int,dest = 'genomicDistance', help = 'genomic distance for local interaction', metavar = '<int>',default=200000)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile1,args.infile2,args.outdir,args.species,args.genomic_distance)
