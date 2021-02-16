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
import random,time
import cooler
import string


def get_local_interaction_from_matrix(matrixfile,tmpdir,chrom_len,genomic_distance):
     
    # for each chr, extract data from .hic file
    coolfile = cooler.Cooler(matrixfile)
    # check the resolution
    bin_start = coolfile.bins()[:1].start
    bin_end = coolfile.bins()[:1].end
    resolution = int(bin_end - bin_start)
    check_bins = int(genomic_distance/resolution)
    
    # check genome assembly
    try:
        assert max(coolfile.chroms()["length"][:]) == chrom_len['chr1']
    except AssertionError:
        print('\n****\nERROR! the length of chr1 ({}) does not equal to {}.'.format(max(coolfile.chroms()["length"][:]),chrom_len['chr1']))
        print('Please make sure you are using the accurate genome assembly: hg38/mm10.\n****\n')
        raise
    

    # ==== read start and end pos for each chrom
    chroms = coolfile.chroms()["name"][:]
    # chroms = [i for i in chroms if 'chr'+i.lstrip('chr') in chrom_len.keys()]
    chr_local_interaction = {}
    chr_ord_id = {}
    for chr in chroms:
        bins_chr = coolfile.bins().fetch(chr)
        chr_start_ID = bins_chr.index[0]
        chr_start_pos = bins_chr.start[chr_start_ID]
        chr_end_ID = bins_chr.index[-1]
        chr_end_pos = bins_chr.start[chr_end_ID]
        # initiate the local interaction dict
        chr_ord_id[chr] = [chr_start_ID,chr_start_pos,chr_end_ID,chr_end_pos]
        chr_local_interaction[chr] = {}
        for bin in np.arange(chr_start_ID,chr_end_ID+1):
            chr_local_interaction[chr][int(bin)] = [0]*check_bins
        

    # == Split the contact matrix pixel records into equally sized chunks to save memory 
    pixels_len = coolfile.pixels().shape[0]
    chunksize=2e7
    total_chunks = int(np.ceil(pixels_len/chunksize))
    for chunk in np.arange(total_chunks):
        chunk_start = int(chunk*chunksize)
        chunk_end = int(min(chunk_start+chunksize,pixels_len))#;print(chunk_start,chunk_end)
        pixels_chunk = coolfile.pixels()[chunk_start:chunk_end]
        # ==== keep only check_bins-limited diagonal pixels
        bin_shift = pixels_chunk['bin2_id']-pixels_chunk['bin1_id']
        filtered_pixels = pixels_chunk[(bin_shift>=0)&(bin_shift<check_bins)]
        print('== processing {}/{} '.format(chunk+1,total_chunks))
        
        for chr in chroms:
            chr_start_ID = chr_ord_id[chr][0]
            chr_end_ID = chr_ord_id[chr][2]
            pixels_chr = filtered_pixels[(filtered_pixels['bin1_id']>=chr_start_ID)&(filtered_pixels['bin1_id']<=chr_end_ID)]
            for pixel_id in pixels_chr.index:
                leftID = pixels_chr.loc[pixel_id,'bin1_id'] 
                rightID = pixels_chr.loc[pixel_id,'bin2_id'] 
                score = pixels_chr.loc[pixel_id,'count'] 
                chr_local_interaction[chr][leftID][rightID-leftID] = score  
            

    # ==== write out the local interaction, for double check
#     flag= os.path.basename(matrixfile).split('.cool')[0]
#     flag = '{}_{}'.format(flag,''.join(random.choice(string.ascii_letters + string.digits) for _ in range(6)))
#     for chr in chroms:
#         if 'chr'+chr.lstrip('chr') in chrom_len.keys():
#             chr_start_id = chr_ord_id[chr][0]
#             chr_start_pos = chr_ord_id[chr][1]
#             outfile = tmpdir+os.sep+'tmp_{}_chr{}_res_{}_view_region_{}.csv'.format(flag,chr.lstrip('chr'),resolution,genomic_distance)
#             outf = open(outfile,'w')
#             outf.write('{}\t{}\n'.format('dis','\t'.join(map(str,np.arange(0,genomic_distance,resolution)))))
#             for view_id in chr_local_interaction[chr].keys():
#                 view_pos = chr_start_pos+(view_id-chr_start_id)*resolution
#                 outf.write('{}\t{}\n'.format(view_pos,'\t'.join(map(str,chr_local_interaction[chr][view_id]))));
#             outf.close()

    # save the local interaction
    interaction_dfs = {} # interaction dataframe
    for chr in chroms:
        if 'chr'+chr.lstrip('chr') in chrom_len.keys():
            chr_start_id = chr_ord_id[chr][0]
            chr_start_pos = chr_ord_id[chr][1]
            chr_key='chr{}'.format(chr.lstrip('chr'))
            interaction_dfs[chr_key] = pd.DataFrame.from_dict(chr_local_interaction[chr],orient='index')
            columns = np.arange(0,genomic_distance,resolution)
            index = [chr_start_pos+(view_id-chr_start_id)*resolution for view_id in interaction_dfs[chr_key].index]
            interaction_dfs[chr_key].columns = columns
            interaction_dfs[chr_key].index = index
                   
    return interaction_dfs,resolution




def main(args):
    
    args.infile='a6010_r5000.cool'
    args.infile = 'pd9_r50000.cool'
    os.makedirs(args.outdir,exist_ok=True)
    
    chrom_len = GenomeData.species_chrom_lengths[args.species] 
    matrixfile = args.infile
    tmpdir = args.outdir
    genomic_distance = args.genomicDistance
    resolution = 50000
    
    # for input ord/matrix file, get the local interaction for each chrom  
    flag,resolution = get_local_interaction_from_matrix(matrixfile,tmpdir,chrom_len,genomic_distance)




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--infile', action = 'store', type = str,dest = 'infile', help = 'input matrix file', metavar = '<dir>')
    parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    parser.add_argument('-v', '--genomicDistance', action = 'store', type = int,dest = 'genomicDistance', help = 'genomic distance for local interaction', metavar = '<int>',default=200000)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args)
