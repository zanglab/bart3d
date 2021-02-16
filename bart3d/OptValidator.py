'''
Copyright (c) 2020, 2021 Zhenjia Wang <zhenjia@virginia.edu>, Chongzhi Zang <zang@virginia.edu>

This code is free software; you can redistribute it and/or modify it 
under the terms of the BSD License.

@author: Zhenjia Wang, Chongzhi Zang
@contact: zhenjia@virginia.edu, zang@virginia.edu
'''

import configparser
import os,sys,re

script_dir = os.path.dirname(os.path.realpath(__file__))
# script_dir = os.path.dirname(script_dir)

def conf_validate():
    '''
    Read user provided path from 'bart.conf' config file
    '''
    config = configparser.ConfigParser()
    config_path = os.path.join(script_dir, 'bart3d.conf')
    if not os.path.exists(config_path):
        sys.stderr.write("CRITICAL: bart3d.conf does not exist in {}!\n".format(script_dir))
        sys.exit(1)
    config.read(config_path)
    return config

def opt_validate(options):
    '''
    Validate input options and specify used data.
    '''
    config = conf_validate()
    
    if not options.outdir:
        options.outdir = os.path.join(os.getcwd(), 'bart3d_output') # create output directory at current working directory
    
    if not options.outFileName:
        # input only contains .bam/.bed/.txt
        first_treatment_file=options.treatment.split(',')[0]
        first_control_file=options.control.split(',')[0]
        treat_base = re.split('.matrix|.hic|.cool',os.path.basename(first_treatment_file))[0]
        control_base = re.split('.matrix|.hic|.cool',os.path.basename(first_control_file))[0]
        outfile_base = '{}_OVER_{}'.format(treat_base,control_base)
        options.outFileName = outfile_base

    
    # === hg38 ===
    if options.species == 'hg38':   
        data_dir = os.path.join(config['path']['hg38_library_dir'], 'hg38_library')

    # === mm10 ===
    elif options.species == 'mm10': 
        data_dir = os.path.join(config['path']['mm10_library_dir'], 'mm10_library')
        
    options.normfile = data_dir+os.sep+'bart2_{}_H3K27ac.dat'.format(options.species)
    options.dhsfile = data_dir+os.sep+'bart2_{}_UDHS.bed'.format(options.species)
    options.tffile = data_dir+os.sep+'bart2_{}_TF_file.json'.format(options.species)
    options.tfoverlap = data_dir+os.sep+'bart2_{}_TF_overlap.json'.format(options.species)
    options.chromsize = script_dir+os.sep+'utility'+os.sep+'{}_clean.chrom.sizes'.format(options.species)
    options.BgToBigWig = script_dir+os.sep+'utility'+os.sep+'bedGraphToBigWig'
    
    return options

