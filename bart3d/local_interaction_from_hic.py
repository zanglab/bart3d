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
from bart3d import straw
import struct,string


####### COPY START #######

"""
The following code was copied/revised from hic2cool, a stand-alone Python file written by 
Carl Vitzthum of the HMS DBMI Park lab. Originally published 1/26/17.

This code was originally based on the straw project by Neva C. Durand and 
Yue Wu (https://github.com/theaidenlab/straw).
"""

def readcstr(f):
    # buf = bytearray()
    buf = b""
    while True:
        b = f.read(1)
        if b is None or b == b"\0":
            # return buf.encode("utf-8", errors="ignore")
            return buf.decode("utf-8")
        elif b == "":
            raise EOFError("Buffer unexpectedly empty while trying to read null-terminated string")
        else:
            buf += b
            
def read_header(infile):
    """
    Takes in a .hic file and returns a dictionary containing information about
    the chromosome. Keys are chromosome index numbers (0 through # of chroms
    contained in file) and values are [chr idx (int), chr name (str), chrom
    length (str)]. Returns the masterindex used by the file as well as the open
    file object.
    """
    req = open(infile, 'rb')
    chrs = {}
    resolutions = []
    magic_string = struct.unpack(b'<3s', req.read(3))[0]
    req.read(1)
    if (magic_string != b"HIC"):
        error_string = ('... This does not appear to be a HiC file; '
                       'magic string is incorrect')
        req.close()
        print(error_string, file=sys.stderr)
        sys.exit(1)
        # force_exit(error_string, req)
    global version
    version = struct.unpack(b'<i', req.read(4))[0]
    masterindex = struct.unpack(b'<q', req.read(8))[0]
    genome = b""
    c = req.read(1)
    while (c != b'\0'):
        genome += c
        c = req.read(1)
    genome = genome.decode('ascii')
    # metadata extraction
    metadata = {}
    nattributes = struct.unpack(b'<i', req.read(4))[0]
    for x in range(nattributes):
        key = readcstr(req)
        value = readcstr(req)
        metadata[key] = value
    nChrs = struct.unpack(b'<i', req.read(4))[0]
    for i in range(0, nChrs):
        name = readcstr(req)
        length = struct.unpack(b'<i', req.read(4))[0]
        if name and length:
            chrs[i] = [i, name, length]
    nBpRes = struct.unpack(b'<i', req.read(4))[0]
    # find bp delimited resolutions supported by the hic file
    for x in range(0, nBpRes):
        res = struct.unpack(b'<i', req.read(4))[0]
        resolutions.append(res)
    return chrs, resolutions #, masterindex, genome, metadata
    

####### COPY END #######


def get_local_interaction_from_matrix(matrixfile,resolution,tmpdir,chrom_len,genomic_distance):
     
    # for each chr, extract data from .hic file
    
    # ==== for each chr, get the first/end ID in ord file 
    # ==== read out the genomic_distance limited interactions for each bin
    chroms_info,resolutions = read_header(matrixfile)
    # print(chroms_info,resolutions);exit()
    # == check if given resolution is supported, otherwise use the smallest one.
    if resolution not in resolutions:
        print('\n****\nERROR! The given resolution {} (in bp) is not supported in {}.\n'.format(resolution,matrixfile))
        print('Available resolutions are {}.\n****\n'.format(resolutions))
        exit(1)
    # resolution = min(resolutions) if resolution==0 else resolution 
    print('== Use resolution {} (in bp)'.format(resolution))
    
    # ==== check the genome assembly. Either '1' or 'chr1' is used in .hic file
    check_chr = 'chr{}'.format(chroms_info[1][1].lstrip('chr'))
    try:
        assert chroms_info[1][2]==chrom_len[check_chr]
    except AssertionError:
        print('\n****\nERROR! the length of {} ({}) does not equal to {}.'.format(check_chr,chroms_info[1][2],chrom_len[check_chr]))
        print('Please make sure you are using the accurate genome assembly: hg38/mm10.\n****\n')
        exit(1)


#     flag= os.path.basename(matrixfile).split('.hic')[0]
#     flag = '{}_{}'.format(flag,''.join(random.choice(string.ascii_letters + string.digits) for _ in range(6)))

    check_bins = int(genomic_distance/resolution)
    chroms = [chroms_info[i][1] for i in chroms_info.keys() if chroms_info[i][1].lower()!='all']
    # chroms = [i for i in chroms if 'chr'+i.lstrip('chr') in chrom_len.keys()]
    
    chr_local_interaction = {} # local interaction score for each chr
    interaction_dfs = {} # interaction dataframe
    for chr in chroms:
        if 'chr'+chr.lstrip('chr') not in chrom_len.keys():
            continue
        print('== Extrac data from chrom {}'.format(chr))
        results = straw.straw('NONE',matrixfile,chr,chr,'BP',resolution)
        # == the end position
        chr_end_pos = max(set(results[0]+results[1]))
        chr_local_interaction[chr] = {}
        # == initiate the data for local interaction
        for bin in np.arange(0,chr_end_pos+1,resolution):
            chr_local_interaction[chr][int(bin)] = [0]*check_bins

        # == read contact matrix
        for result_id in np.arange(len(results[0])):
            leftID = results[0][result_id]
            rightID = results[1][result_id]
            if 0<=rightID-leftID <genomic_distance:
                score = results[2][result_id]
                chr_local_interaction[chr][leftID][int((rightID-leftID)/resolution)] = score

        # ==== write out the local interaction, for double check
#         outfile = tmpdir+os.sep+'tmp_{}_chr{}_res_{}_view_region_{}.csv'.format(flag,chr.lstrip('chr'),resolution,genomic_distance)
#         outf = open(outfile,'w')
#         outf.write('{}\t{}\n'.format('dis','\t'.join(map(str,np.arange(0,genomic_distance,resolution)))))
#         for view_pos in chr_local_interaction[chr].keys():
#             outf.write('{}\t{}\n'.format(view_pos,'\t'.join(map(str,chr_local_interaction[chr][view_pos]))));
#         outf.close()
    
        # save the local interaction
        chr_key='chr{}'.format(chr.lstrip('chr'))
        interaction_dfs[chr_key] = pd.DataFrame.from_dict(chr_local_interaction[chr],orient='index')
        columns = np.arange(0,genomic_distance,resolution)
        interaction_dfs[chr_key].columns = columns
                   
    return interaction_dfs,resolution




def main(args):
    
    # ==== test demo
    args.infile= '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/f8_hic_TADs/f0_hic_data_hg38/view_hic/A6010_inter_30.hic'  
    args.infile='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/patient_hic_compr_GSE115896/f1_hicpro_to_juicer_new/f1_hicpro_to_juicer_input/PD9_allValidPairs.hic'
    # args.infile='GSE71831_Patski_paternal.hic'
    os.makedirs(args.outdir,exist_ok=True)
    
    chrom_len = GenomeData.species_chrom_lengths[args.species] 
    matrixfile = args.infile
    tmpdir = args.outdir
    genomic_distance = args.genomicDistance
    resolution = 50000
    
    # for input ord/matrix file, get the local interaction for each chrom  
    flag,resolution = get_local_interaction_from_matrix(matrixfile,resolution,tmpdir,chrom_len,genomic_distance)




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
