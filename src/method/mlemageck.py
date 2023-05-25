'''
MAGeCK MLE main entry
'''

from __future__ import print_function
from mageck.mleclassdef import *
import numpy as np
import pandas as pd
import sys

def create_design_matrix(args, df):
    '''
    post-processing of argument parsing
    '''
    
    # read screen to sequence map file
    screen_sequence_map=pd.read_csv(args.screen_sequence_map)

    # create replicate to model dictionary
    replicate_model_dict={}

    for i in range(len(screen_sequence_map)):
        row=screen_sequence_map.iloc[i]

        if row['SequenceID'] not in replicate_model_dict.keys():
            if row['ModelID'] != 'pDNA': # exclude pDNA (control)
                replicate_model_dict[row['SequenceID']]=row['ModelID']

    # create design matrix
    sequence_id = np.unique(df.columns[2:])
    sequence_id = np.intersect1d(sequence_id, screen_sequence_map['SequenceID'])
    
    ## get model id from sequence id
    model_id = []

    for i in range(len(sequence_id)):
        if sequence_id[i] not in model_id:
            model_id.append(replicate_model_dict[sequence_id[i]])
    model_id = np.unique(model_id)
    model_id = np.append('baseline', model_id)

    ## initialize pandas dataframe
    desmat=pd.DataFrame(0, index=sequence_id,columns=model_id)

    ## fill in the design matrix
    for i in range(len(desmat.columns)):
        col_id = desmat.columns[i]
        
        if col_id == 'baseline':
            desmat.loc[:, col_id] = 1
        else:
            seq_id = screen_sequence_map.loc[screen_sequence_map['ModelID'] == col_id, 'SequenceID']
            seq_id = np.intersect1d(seq_id, desmat.index)

            if len(seq_id) == 0:
                print('No sequence id found for model id: ' + col_id)
                sys.exit(-1)
            
            desmat.loc[seq_id, col_id] = 1

    args.design_matrix=desmat.to_numpy()
    
    # parsing sample label
    args.include_samples=list(desmat.index)
    args.beta_labels=list(desmat.columns)
    
    return args


def read_gene_from_file(df,includesamples=None):
    '''
    Reading gene models 
    Parameters
    ----------
    filename
        file name of the read count table 
    includesamples
        If not None, only samples in the includesampels are included
    '''
    # first, read count table
    allgenedict={}
    nline=0
    nsamples=0
    ngene=0
    
    sampleids=[]
    sampleindex=[]
    sampleids_toindex={}
        
    for line in range(len(df)):
        field=df.iloc[line]

        if nline == 0:
            # The first line: check sample columns
            header=list(df.columns)
            nsamples=len(header)-2
            sampleids=list(header[2:])
            for i in range(nsamples):
                sampleids_toindex[sampleids[i]]=i
            if includesamples != None:
                for si in includesamples:
                    if si not in sampleids_toindex:
                        sys.exit(-1)
                sampleindex=[sampleids_toindex[si] for si in includesamples]
            else:
                sampleindex=[i for i in range(nsamples)]

        sgid=str(field[0])
        gid=str(field[1])

        if gid not in allgenedict:
            sks=SimCaseSimple()
            sks.prefix=gid
            sks.nb_count=[]
            sks.sgrnaid=[]
            ngene+=1
            for i in sampleindex:
                sks.nb_count+=[[]]
            allgenedict[gid]=sks
        else:
            sks=allgenedict[gid]
        sks.sgrnaid+=[sgid]
        for i in range(len(sampleindex)): 
            ni=sampleindex[i]
            try:
                nrt=float(field[ni+2])+1 # add 1 pseudocount
                sks.nb_count[i]+=[nrt]
            except ValueError:
                print('Error loading line '+str(nline))

        nline+=1

    # convert nb_count to matrix
    for (gene,ginst) in allgenedict.items():
        ginst.nb_count=np.matrix(ginst.nb_count)

    return allgenedict


def mageckmle_main(parsedargs=None,returndict=False):
    '''
    Main entry for MAGeCK MLE
    ----
    Parameters:

    parsedargs
        Arguments for parsing
    returndict
        If set true, will not try to run the whole prediction process, but will return after mean variance modeling
    '''
        
    import numpy as np
    from mageck.mleinstanceio import write_gene_to_file
    from mageck.mleem import iteratenbem
    from mageck.mlemeanvar import MeanVarModel
    from mageck.mageckCount import normalizeCounts
    from mageck.mlesgeff import read_sgrna_eff
    from mageck.mlemultiprocessing import runem_multiproc,iteratenbem_permutation,iteratenbem_permutation_by_nsg
    from mageck.cnv_normalization import read_CNVdata,betascore_piecewisenorm,betascore_piecewisenorm,highestCNVgenes
    from mageck.cnv_estimation import mageckmleCNVestimation

    # main process
    args=parsedargs
    maxfittinggene=args.genes_varmodeling
    maxgene=np.inf

    # reading sgRNA efficiency
    read_sgrna_eff(args)

    # reading read count table
    count_table=pd.read_csv(args.count_table)

    # parsing arguments
    args=create_design_matrix(args, count_table)

    # create gene dictionary
    allgenedict=read_gene_from_file(count_table,includesamples=args.include_samples)
    
    # check the consistency of negative control sgRNAs
    sgrna2genelist={}
    for (geneid,gk) in allgenedict.items():
        sgid=gk.sgrnaid
        for sg_i in sgid:
            sgrna2genelist[sg_i]=geneid

    # calculate the size factor
    cttab_sel={}
    for (geneid,gk) in allgenedict.items():
        sgid=gk.sgrnaid
        sgreadmat=gk.nb_count.getT().tolist()
        for i in range(len(sgid)):
            cttab_sel[sgid[i]]=sgreadmat[i]
    if hasattr(args,'norm_method'):
        if args.norm_method!='none':
            size_f=normalizeCounts(cttab_sel,method=args.norm_method,returnfactor=True,reversefactor=True,controlsgfile=args.control_sgrna)
        else:
            size_f=None
    else:
        size_f=normalizeCounts(cttab_sel,returnfactor=True,reversefactor=True)

    desmat=args.design_matrix
    ngene=0
    for (tgid,tginst) in allgenedict.items():
        tginst.design_mat=desmat
    
    # perform copy number estimation if option selected
    CN_arr = None
    CN_celldict = None
    CN_genedict = None
    genes2correct = False
    CN_celllabel=args.beta_labels[1:]
    if args.cnv_norm is not None or args.cnv_est is not None: 
        # check if --cell-line option is set
        if args.cell_line is not None:
            CN_celllabel=[args.cell_line]*len(args.beta_labels[1:]) # replace it with all cell lines provided
        if args.cnv_norm is not None: 
            # get copy number data from external copy number dataset
            # here is used just check the cnv files
            (CN_arr,CN_celldict,CN_genedict) = read_CNVdata(args.cnv_norm,CN_celllabel)
            genes2correct = False # do not select only subset of genes to correct (i.e. correct all genes)
        elif args.cnv_est is not None:
            # estimating CNVS
            # organize sgRNA-gene pairing into dictionary
            sgrna2genelist = {sgrna: gene for gene in allgenedict for sgrna in allgenedict[gene].sgrnaid}
            # estimate CNV and write results to file
            mageckmleCNVestimation(args.cnv_est,cttab_sel,desmat,sgrna2genelist,CN_celllabel,args.output_prefix)
            # read into the data structures
            (CN_arr,CN_celldict,CN_genedict) = read_CNVdata(str(args.output_prefix)+'.CNVestimates.txt',CN_celllabel)
            genes2correct = highestCNVgenes(CN_arr,CN_genedict,percentile=98)
        # if no match was found
    
    # run the EM for a few genes to perform gene fitting process
    meanvardict={}
    for (tgid,tginst) in allgenedict.items():
        #iteratenbem(tginst,debug=False,alpha_val=0.01,estimateeff=False,size_factor=size_f)
        ##sgrna_wide_dispersion_estimation_MAP_v2(tginst,tginst.design_mat)
        ngene+=1
        tginst.w_estimate=[]
        meanvardict[tgid]=tginst
        if ngene>maxfittinggene:
            break

    argsdict={'debug':False, 'alpha_val':0.01, 'estimateeff':False,'size_factor':size_f}
    runem_multiproc(meanvardict,args,nproc=args.threads,argsdict=argsdict)

    for (tgid,tginst) in meanvardict.items():
        allgenedict[tgid]=tginst

    # model the mean and variance
    if maxfittinggene>0:
        mrm=MeanVarModel()
        # old: linear model
        mrm.get_mean_var_residule(allgenedict)
        mrm.model_mean_var_by_lm()
        # new method: generalized linear model
        #mrm.model_mean_disp_by_glm(allgenedict,args.output_prefix,size_f)
    else:
        mrm=None
    
    if returndict:
        return (allgenedict,mrm,size_f)
    
    # run the test again...
    if hasattr(args,'threads') and args.threads>0:
        # multipel threads
        argsdict={'debug':False,'estimateeff':True,'meanvarmodel':mrm,'restart':False,'removeoutliers':args.remove_outliers,'size_factor':size_f,'updateeff':args.update_efficiency,'logem':False}
        runem_multiproc(allgenedict,args,nproc=args.threads,argsdict=argsdict)
    else:
        # only 1 thread
        # the following codes should be merged to the above code section
        ngene=0
        for (tgid,tginst) in allgenedict.items():
            #try:
            if hasattr(args,'debug_gene') and args.debug_gene!=None and tginst.prefix != args.debug_gene:
                continue
            iteratenbem(tginst,debug=False,estimateeff=True,meanvarmodel=mrm,restart=False,removeoutliers=args.remove_outliers,size_factor=size_f,updateeff=args.update_efficiency)
            # Tracer()()
            ngene+=1
            if ngene>maxgene:
                break

    # set up the w vector
    for (tgid,tginst) in allgenedict.items():
        if len(tginst.w_estimate)==0:
            tginst.w_estimate=np.ones(len(tginst.sgrnaid))
    #Tracer()()
    # permutation, either by group or together
    if args.no_permutation_by_group:
        iteratenbem_permutation(allgenedict,args,nround=args.permutation_round,removeoutliers=args.remove_outliers,size_factor=size_f)
    else:
        iteratenbem_permutation_by_nsg(allgenedict,args,size_f=size_f)
    # correct for FDR
    from mageck.mleclassdef import gene_fdr_correction;
    gene_fdr_correction(allgenedict,args.adjust_method);
    # correct for CNV
    if args.cnv_norm is not None or args.cnv_est is not None:
        betascore_piecewisenorm(allgenedict,CN_celllabel,CN_arr,CN_celldict,CN_genedict,selectGenes=genes2correct)

    # write to file
    genefile=args.output_prefix+'.gene_summary.txt'
    # sgrnafile=args.output_prefix+'.sgrna_summary.txt'
    write_gene_to_file(allgenedict,genefile,args,betalabels=args.beta_labels)
    # write_sgrna_to_file(allgenedict,sgrnafile)
    return (allgenedict,mrm)


if __name__ == '__main__':
    try:
        mageckmle_main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) Bye!\n")
        sys.exit(0)


