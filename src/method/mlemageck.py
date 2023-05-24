'''
MAGeCK MLE main entry
'''

from __future__ import print_function
import sys

def mageckmle_postargs(args):
    '''
    post-processing of argument parsing
    '''
 
    from mageck.mledesignmat import parse_designmat,parse_designmat_from_day0
    
    try:
        import scipy
        from scipy.stats import nbinom
    except ImportError:
        print('Cannot find scipy (required for mle approach). Please check your scipy installation.')
        sys.exit(-1)
    try:
        import numpy as np
        import numpy.linalg as linalg
    except ImportError:
        print('Cannot find numpy (required for mle approach). Please check your numpy installation.')
        sys.exit(-1)
    # parsing design matrix
    if args.design_matrix != None:
        (desmat,sampleid,betalabel)=parse_designmat(args.design_matrix)
    else:
        (desmat,sampleid,betalabel)=parse_designmat_from_day0(args)

    args.design_matrix=desmat
    
    # parsing sample label
    if sampleid ==None:
        # design matrix is provided as a string
        if args.include_samples !=None:
            args.include_samples=args.include_samples.split(',')
            if len(args.include_samples) != desmat.shape[0]:
                print('The number of samples in the --include-samples option do not match rows in design matrix.')
                sys.exit(-1)
        if args.beta_labels!=None:
            args.beta_labels=args.beta_labels.split(',')
            if len(args.beta_labels) != desmat.shape[1]:
                print('The number of labels in the --beta-labels option do not match columns in design matrix.')
                sys.exit(-1)
    else:
        # design matrix is provided as file
        if args.include_samples !=None:
            print('Sample labels are included in the design matrix file '+args.design_matrix+'. The --include-samples option should not be used.')
            sys.exit(0)
        if args.beta_labels!=None:
            print('Beta labels are included in the design matrix file '+args.design_matrix+'. The --beta-labels option should not be used.')
            sys.exit(0)

        args.include_samples=sampleid
        args.beta_labels=betalabel

        if len(args.include_samples) != desmat.shape[0]:
            print('The number of samples in the --include-samples option do not match rows in design matrix.')
            sys.exit(-1)
        if len(args.beta_labels) != desmat.shape[1]:
            print('The number of labels in the --beta-labels option do not match columns in design matrix.')
            sys.exit(-1)
    
    return args


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
    # parsing arguments
    args=mageckmle_postargs(parsedargs)
        
    import scipy
    from scipy.stats import nbinom
    import numpy as np
    import numpy.linalg as linalg
    from mageck.mleinstanceio import read_gene_from_file,write_gene_to_file,write_sgrna_to_file
    from mageck.mleem import iteratenbem
    from mageck.mlemeanvar import MeanVarModel
    from mageck.mageckCount import normalizeCounts
    from mageck.mlesgeff import read_sgrna_eff,sgrna_eff_initial_guess
    from mageck.dispersion_characterization import sgrna_wide_dispersion_estimation_MAP_v2
    from mageck.mlemultiprocessing import runem_multiproc,iteratenbem_permutation,iteratenbem_permutation_by_nsg
    from mageck.cnv_normalization import read_CNVdata,match_sgrnaCN,betascore_piecewisenorm,betascore_piecewisenorm,highestCNVgenes
    from mageck.cnv_estimation import mageckmleCNVestimation
    from mageck.mageckCount import mageckcount_checkcontrolsgrna

    # main process
    maxfittinggene=args.genes_varmodeling
    maxgene=np.inf
    # reading sgRNA efficiency
    read_sgrna_eff(args)
    # reading read count table
    allgenedict=read_gene_from_file(args.count_table,includesamples=args.include_samples)
    #
    # check the consistency of negative control sgRNAs
    sgrna2genelist={}
    for (geneid,gk) in allgenedict.items():
        sgid=gk.sgrnaid
        for sg_i in sgid:
            sgrna2genelist[sg_i]=geneid

    mageckcount_checkcontrolsgrna(args,sgrna2genelist)
    # 
    sgrna_eff_initial_guess(args,allgenedict)
    #
    #
    #
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

    # desmat=np.matrix([[1,1,1,1],[0,0,1,0],[0,0,0,1]]).getT()
    desmat=args.design_matrix
    ngene=0
    for (tgid,tginst) in allgenedict.items():
        tginst.design_mat=desmat
    #
    #
    #
    #
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

    #
    #
    #
    #
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
    #
    #
    #
    #
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
            #except:
            #    print('Error occurs while calculating beta values of gene '+tgid+'.')
            #    sys.exit(-1)
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


