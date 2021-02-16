'''
Copyright (c) 2020, 2021 Zhenjia Wang <zhenjia@virginia.edu>, Chongzhi Zang <zang@virginia.edu>

This code is free software; you can redistribute it and/or modify it 
under the terms of the BSD License.

@author: Zhenjia Wang, Chongzhi Zang
@contact: zhenjia@virginia.edu, zang@virginia.edu
'''

import os,sys,time
import argparse
import pandas as pd
import numpy as np
from operator import itemgetter
import bisect



def get_overlapped_udhs(start,end,udhs_starts_chr,udhs_ends_chr,udhs_ids_chr):
    # for each interval/row in  the input bed file, return ids of all overlapped UDHS
    s = bisect.bisect_left(udhs_ends_chr,start)
    e = bisect.bisect_right(udhs_starts_chr,end)
    overlapped_ids = udhs_ids_chr[s:e]
    return overlapped_ids
    

def read_regions_from_bed_not_sep_strand(chromsize,dhsfile,species):
    
    chrom_len = pd.read_csv(chromsize,header=None,index_col=0,sep='\t')[1].to_dict()
    udhs_starts,udhs_ends,udhs_ids = {},{},{}

    with open(dhsfile,'r') as inf:
        lines = inf.readlines(100000)
        while lines:
            for line in lines:# and (if not re.match("#",lines)):
                sline = line.strip().split()
                chrom = sline[0]
                start = int(sline[1])
                end = int(sline[2])
                id = int(sline[3])
                if chrom in chrom_len.keys():
                    if chrom not in udhs_starts:
                        udhs_starts[chrom]=[]
                        udhs_ends[chrom]=[]
                        udhs_ids[chrom]=[]
                    udhs_starts[chrom].append(start)
                    udhs_ends[chrom].append(end)
                    udhs_ids[chrom].append(id)
                else:
                    pass   
            #print(regions);exit(0)
            lines = inf.readlines(100000)  
    return udhs_starts,udhs_ends,udhs_ids



def score_on_DHS(args):
    '''
    Assign the score on bed interval to the overlapping UDHS
    '''

    # read regions from UDHS
    udhs_starts,udhs_ends,udhs_ids = read_regions_from_bed_not_sep_strand(args.chromsize,args.dhsfile,args.species)
    #print(udhs_ids['chr1'][:10]);exit()
    
    #read interval from input bed file
    bedfile = args.diff_file
    bfile = open(bedfile,'r')
    line = bfile.readline()
    # initiate the score on each of UDHS
    counting = {}
    for chrom in udhs_ids:
        for udhs_id in udhs_ids[chrom]:
            counting[udhs_id]=[]
    #print(counting);exit()
    while line:
        # for each row of the input bed file, assign score to each of the overlapped DHS
        line = line.strip().split()
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        # score from specific column
        try:
            score = float(line[3])
        except:
            score = 0
        if chrom in udhs_starts:
            # get the score on this DHS, 0 if not overlapping with inpit bed file
            if start < end:
                overlapped_ids = get_overlapped_udhs(start,end,udhs_starts[chrom],udhs_ends[chrom],udhs_ids[chrom])
            else:
                sys.stdout.write("Warning: start positions should be smaller than end positions!\n")
            #print(overlapped_ids);exit()
            for overlapped_id in overlapped_ids:
                counting[overlapped_id].append(score)
            
        line = bfile.readline()
    bfile.close()
    #print(counting);exit()
    # get the averaged value 
    for count_id in counting:
        if len(counting[count_id])==0:
            counting[count_id]=0.0
        else:
            counting[count_id]=np.mean(counting[count_id])
    return counting

    # with open('test2.out','w') as outf:
    #     for count_id in counting:
    #         outf.write('{}\t{}\n'.format(count_id,counting[count_id]))



def main(args):


    args.dhsfile = '/nv/vol190/zanglab/zw5j/data/unionDHS/hg38_unionDHS_fc4_50merge.bed'
    #args.dhsfile = 'test.bed'
    score_on_DHS(args)




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input bed file', metavar = '<file>')
    parser.add_argument('-c', '--scorecol', action = 'store', type = int,dest = 'scorecol', help = 'specify the column for score', metavar = '<file>',default=3)
    parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args)
