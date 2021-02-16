'''
Copyright (c) 2020, 2021 Zhenjia Wang <zhenjia@virginia.edu>, Chongzhi Zang <zang@virginia.edu>

This code is free software; you can redistribute it and/or modify it 
under the terms of the BSD License.

@author: Zhenjia Wang, Chongzhi Zang
@contact: zhenjia@virginia.edu, zang@virginia.edu
'''

import os,sys
from bart3d import OptValidator, get_differential_score_per_bin, differential_score_on_UDHS, AUCcalc, StatTest, DCI_normalization

def stat_by_position(args,positions,name_flag):

    # ==== Start calculating the AUC score for each TF ChIP-seq dataset...
    sys.stdout.write('== Start calculating the AUC score for each TR ChIP-seq dataset...\n')
    sys.stdout.flush()

    auc_file = args.outdir+os.sep+args.outFileName + '_{}_auc.txt'.format(name_flag)
    tf_aucs, tf_index = AUCcalc.cal_auc(args, positions,auc_file)
    

    # ==== Statistical tests start
    sys.stdout.write('== Statistical tests start.\n')
    stat_file = args.outdir + os.sep+args.outFileName + '_{}_bart3d_results.txt'.format(name_flag)
    StatTest.stat_test(tf_aucs, tf_index, stat_file, args.normfile)



def run_bart3d(options):

    args = OptValidator.opt_validate(options)
#     print(args);exit()
    # ==== check if the library exist
    if not os.path.isfile(args.tffile):
        sys.stderr.write('Oops! No data under the library directory: {}.\n'.format(os.path.dirname(os.path.dirname(args.tffile))))
        sys.stderr.write('Please download the {} libraries from https://faculty.virginia.edu/zanglab/bart/.\n'.format(args.species))
        sys.stderr.write('Please revise the bart3d.conf before installation!\n')
        sys.exit(1)  
    else:      
        sys.stderr.write("Library directory: {}\n".format(os.path.dirname(os.path.dirname(args.tffile))))
        
    # ==== create output directory
    try:
        os.makedirs(args.outdir, exist_ok=True)
    except:
        sys.stderr.write('Output directory: {} could not be created. \n'.format(args.outdir))
        sys.exit(1)
    sys.stdout.write("Output directory: {} \n".format(args.outdir))
    sys.stdout.write("Output file prefix: {} \n\n".format(os.path.basename(args.outFileName)))

    if args.species == 'hg38':
        sys.stdout.write("Start prediction on hg38...\n")
    elif args.species == 'mm10':
        sys.stdout.write("Start prediction on mm10...\n")
    sys.stdout.flush()


    # ==== file of differential score
    sys.stdout.write('Start calculating differential chromatin interactions...\n')
    sys.stdout.flush()
    args.diff_file,treat_matrix_data,control_matrix_data = get_differential_score_per_bin.differential_score(args)  
    # args.diff_file = '{}/{}_differential_score.bed'.format(args.outdir,args.outFileName)
    if args.coverageNormalization:
        args.diff_file=DCI_normalization.profile_normalization(args,treat_matrix_data,control_matrix_data)
        
    
    # ==== Start Sorting scored UDHS...
    sys.stdout.write('\nStart sorting scored UDHS...\n')
    sys.stdout.flush()
    counting = differential_score_on_UDHS.score_on_DHS(args)#;print(counting)
    
    # TRs associated with chromatin interaction increased regions
    sys.stdout.write('Predicting TRs associated with INCREASED chromatin interactions...\n')
    positions = sorted(counting.keys(),key=counting.get,reverse=True)
    stat_by_position(args,positions,'Interaction_Increased')
    
    sys.stdout.write('Predicting TRs associated with DECREASED chromatin interactions...\n')
    positions = sorted(counting.keys(),key=counting.get,reverse=False)
    stat_by_position(args,positions,'Interaction_Decreased')

    sys.stdout.write("Congratulations! BART3D job finished successfully!\n\n")




