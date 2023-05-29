"""
MAGeCK argument parser
"""

from __future__ import print_function
import sys
import argparse
import os
import pandas as pd
import tempfile
from mlemageck import *


def arg_count(subparser):
    """
    Adding parsers for count command
    """
    subp_count=subparser.add_parser('count',help='Collecting read counts from fastq files.')
    
    cnt_reqgroup=subp_count.add_argument_group(title='Required arguments',description='')
    cnt_reqgroup.add_argument('-l','--list-seq',required=False,help='A file containing the list of sgRNA names, their sequences and associated genes. Support file format: csv and txt. Provide an empty file for collecting all possible sgRNA counts.')
    subp_count_mg=cnt_reqgroup.add_mutually_exclusive_group(required=True)
    subp_count_mg.add_argument('--fastq',nargs='+',help='Sample fastq files (or fastq.gz files, or SAM/BAM files after v0.5.5), separated by space; use comma (,) to indicate technical replicates of the same sample. For example, "--fastq sample1_replicate1.fastq,sample1_replicate2.fastq sample2_replicate1.fastq,sample2_replicate2.fastq" indicates two samples with 2 technical replicates for each sample.')
    subp_count_mg.add_argument('-k', '--count-table',help='The read count table file. Only 1 file is accepted.')
    
    cnt_norm_group=subp_count.add_argument_group(title='Optional arguments for normalization',description='')
    cnt_norm_group.add_argument('--norm-method',choices=['none','median','total','control'],default='median',help='Method for normalization, including "none" (no normalization), "median" (median normalization, default), "total" (normalization by total read counts), "control" (normalization by control sgRNAs specified by the --control-sgrna option).')
    cnt_norm_group_ctrl=cnt_norm_group.add_mutually_exclusive_group()
    cnt_norm_group_ctrl.add_argument('--control-sgrna',help='A list of control sgRNAs for normalization and for generating the null distribution of RRA.')
    cnt_norm_group_ctrl.add_argument('--control-gene',help='A list of genes whose sgRNAs are used as control sgRNAs for normalization and for generating the null distribution of RRA.')
    
    cnt_iogroup=subp_count.add_argument_group(title='Optional arguments for input and output',description='')
    cnt_iogroup.add_argument('--sample-label',default='',help='Sample labels, separated by comma (,). Must be equal to the number of samples provided (in --fastq option). Default "sample1,sample2,...".')
    cnt_iogroup.add_argument('-n','--output-prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
    cnt_iogroup.add_argument('--unmapped-to-file',action='store_true',help='Save unmapped reads to file, with sgRNA lengths specified by --sgrna-len option.')
    cnt_iogroup.add_argument('--keep-tmp',action='store_true',help='Keep intermediate files.')
    cnt_iogroup.add_argument('--test-run',action='store_true',help='Test running. If this option is on, MAGeCK will only process the first 1M records for each file.')

    cnt_fqgroup=subp_count.add_argument_group(title='Optional arguments for processing fastq files',description='')
    cnt_fqgroup.add_argument('--fastq-2', nargs='+',
                                                            help='Paired sample fastq files (or fastq.gz files), the order of which should be consistent with that in fastq option.')
    cnt_fqgroup.add_argument('--count-pair', default='False',
                                                     help='Report all valid alignments per read or pair (default: False). ')
    cnt_fqgroup.add_argument('--trim-5',default='AUTO',help='Length of trimming the 5\' of the reads. Users can specify multiple trimming lengths, separated by comma (,); for example, "7,8". Use "AUTO" to allow MAGeCK to automatically determine the trimming length. Default AUTO.')
    cnt_fqgroup.add_argument('--sgrna-len',type=int,default=20,help='Length of the sgRNA. Default 20. ATTENTION: after v0.5.3, the program will automatically determine the sgRNA length from library file; so only use this if you turn on the --unmapped-to-file option.')
    cnt_fqgroup.add_argument('--count-n',action='store_true',help='Count sgRNAs with Ns. By default, sgRNAs containing N will be discarded.')
    cnt_fqgroup.add_argument('--reverse-complement',action='store_true',help='Reverse complement the sequences in library for read mapping.')
    
    cnt_qcgroup=subp_count.add_argument_group(title='Optional arguments for quality controls',description='')
    cnt_qcgroup.add_argument('--pdf-report',action='store_true',help='Generate pdf report of the fastq files.')
    cnt_qcgroup.add_argument('--day0-label',help='Turn on the negative selection QC and specify the label for control sample (usually day 0 or plasmid). For every other sample label, the negative selection QC will compare it with day0 sample, and estimate the degree of negative selections in essential genes.')
    cnt_qcgroup.add_argument('--gmt-file',help='The pathway file used for QC, in GMT format. By default it will use the GMT file provided by MAGeCK.')
    
def arg_mle(subparser):
    """
    Adding arguments for MLE
    """
    subm_mle=subparser.add_parser('mle',help='Perform MLE estimation of gene essentiality.')
    # required parameters
    reqgroup=subm_mle.add_argument_group(title='Required arguments',description='')
    reqgroup.add_argument('-k','--count-table',required=True,help='Provide a tab-separated count table. Each line in the table should include sgRNA name (1st column), target gene (2nd column) and read counts in each sample.')
    reqgroup.add_argument('-s', '--screen-sequence-map',help='Sequence to screen map file.')
    reqgroup.add_argument('-d','--design-matrix',help='Provide a design matrix, either a file name or a quoted string of the design matrix. For example, "1,1;1,0". The row of the design matrix must match the order of the samples in the count table (if --include-samples is not specified), or the order of the samples by the --include-samples option.')
    # optional arguments
    #opgroup=subm_mle.add_argument_group(title='Optional arguments',description='Optional arguments')
    ## input and output
    iogroup=subm_mle.add_argument_group(title='Optional arguments for input and output',description='')
    iogroup.add_argument('-n','--output-prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
    iogroup.add_argument('-i', '--include-samples', help='Specify the sample labels if the design matrix is not given by file in the --design-matrix option. Sample labels are separated by ",", and must match the labels in the count table.')
    iogroup.add_argument('-b', '--beta-labels', help='Specify the labels of the variables (i.e., beta), if the design matrix is not given by file in the --design-matrix option. Should be separated by ",", and the number of labels must equal to (# columns of design matrix), including baseline labels. Default value: "bata_0,beta_1,beta_2,...".')
    #iogroup.add_argument('--control-sgrna',help='A list of control sgRNAs. Permutation will also be done from a list of control sgRNAs (instead of all sgRNAs).')
    iogroup_ctrl=iogroup.add_mutually_exclusive_group()
    iogroup_ctrl.add_argument('--control-sgrna',help='A list of control sgRNAs for normalization and for generating the null distribution of RRA.')
    iogroup_ctrl.add_argument('--control-gene',help='A list of genes whose sgRNAs are used as control sgRNAs for normalization and for generating the null distribution of RRA.')
    # iogroup.add_argument('--cnv-norm',help='A matrix of copy number variation data across cell lines to normalize CNV-biased BetaScores.')
    ## Optional CNV correction arguments
    cnvcorgroup=subm_mle.add_argument_group(title='Optional arguments for CNV correction',description='')
    cnvcorgroup.add_argument('--cnv-norm',help='A matrix of copy number variation data across cell lines to normalize CNV-biased sgRNA scores prior to gene ranking.')
    cnvcorgroup.add_argument('--cell-line',help='The name of the cell line to be used for copy number variation normalization. Must match one of the column names in the file provided by --cnv-norm. This option will overwrite the default case where cell line names are inferred from the columns of the design matrix.')
    cnvcorgroup.add_argument('--cnv-est',help='Estimate CNV profiles from screening data. A BED file with gene positions are required as input. The CNVs of these genes are to be estimated and used for copy number bias correction.')

    ## General options
    mlegroup=subm_mle.add_argument_group(title='Optional arguments for MLE module',description='')
    mlegroup.add_argument('--debug',action='store_true',help='Debug mode to output detailed information of the running.')
    mlegroup.add_argument('--debug-gene',help='Debug mode to only run one gene with specified ID.')
    mlegroup.add_argument('--norm-method',choices=['none','median','total','control'],default='median',help='Method for normalization, including "none" (no normalization), "median" (median normalization, default), "total" (normalization by total read counts), "control" (normalization by control sgRNAs specified by the --control-sgrna option).')
    mlegroup.add_argument('--genes-varmodeling',type=int,default=0,help='The number of genes for mean-variance modeling. Default 0.')
    mlegroup.add_argument('--permutation-round',type=int,default=2,help='The rounds for permutation (interger). The permutation time is (# genes)*x for x rounds of permutation. Suggested value: 10 (may take longer time). Default 2.')
    mlegroup.add_argument('--no-permutation-by-group', action='store_true',help='By default, gene permutation is performed separately, by their number of sgRNAs. Turning this option will perform permutation on all genes together. This makes the program faster, but the p value estimation is accurate only if the number of sgRNAs per gene is approximately the same.')
    mlegroup.add_argument('--max-sgrnapergene-permutation',type=int,default=40,help='Do not calculate beta scores or p vales if the number of sgRNAs per gene is greater than this number. This will save a lot of time if some isolated regions are targeted by a large number of sgRNAs (usually hundreds). Must be an integer. Default 40.')
    mlegroup.add_argument('--remove-outliers', action='store_true', help='Try to remove outliers. Turning this option on will slow the algorithm.')
    mlegroup.add_argument('--threads', type=int,default=1, help='Using multiple threads to run the algorithm. Default using only 1 thread.')
    mlegroup.add_argument('--adjust-method',choices=['fdr','holm','pounds'],default='fdr',help='Method for sgrna-level p-value adjustment, including false discovery rate (fdr), holm\'s method (holm), or pounds\'s method (pounds).')
    ## EM options
    emgroup=subm_mle.add_argument_group(title='Optional arguments for the EM iteration',description='')
    emgroup.add_argument('--sgrna-efficiency',help='An optional file of sgRNA efficiency prediction. The efficiency prediction will be used as an initial guess of the probability an sgRNA is efficient. Must contain at least two columns, one containing sgRNA ID, the other containing sgRNA efficiency prediction.')
    emgroup.add_argument('--sgrna-eff-name-column',type=int,default=0,help='The sgRNA ID column in sgRNA efficiency prediction file (specified by the --sgrna-efficiency option). Default is 0 (the first column).')
    emgroup.add_argument('--sgrna-eff-score-column',type=int,default=1,help='The sgRNA efficiency prediction column in sgRNA efficiency prediction file (specified by the --sgrna-efficiency option). Default is 1 (the second column).')
    emgroup.add_argument('--update-efficiency',action='store_true',help='Iteratively update sgRNA efficiency during EM iteration.')

def create_design_matrix(args, control_name='pDNA'):
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
            if row['ModelID'] != control_name: # exclude pDNA (control)
                replicate_model_dict[row['SequenceID']]=row['ModelID']

    # create design matrix
    read_df = pd.read_csv(args.count_table)
    sequence_id = np.unique(read_df.columns[2:])
    sequence_id = np.intersect1d(sequence_id, screen_sequence_map['SequenceID'])
    
    ## get model id from sequence id
    model_id = []

    for i in range(len(sequence_id)):
        if sequence_id[i] not in model_id:
            if control_name not in sequence_id[i]:
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

    desmat.index.name = 'Samples'

    # Create a temporary file
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    desmat.to_csv(temp_file.name, sep='\t')

    # update design matrix argument
    args.design_matrix = temp_file.name
    
    return args

# Remove ID from genes
def header_cleanup(df, index=True):
    if index:
        df_header = list(df.index)
    else:
        df_header = list(df.columns)

    for gene_idx in range(0, len(df_header)):
        gene = df_header[gene_idx]
        gene = gene.split(' (')[0]
        df_header[gene_idx] = gene

    if index:
        df.index = df_header
    else:
        df.columns = df_header

    return(df)

def format_CNV_data(args):
    '''
    post-processing of argument parsing
    '''
    # read CNV data
    cnv_data = pd.read_csv(args.cnv_norm, index_col=0).T

    # remove ID from genes
    cnv_data = header_cleanup(cnv_data)
    cnv_data = cnv_data.fillna(0)
    cnv_data.index.name = 'SYMBOL'
    
    # Create a temporary file
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    cnv_data.to_csv(temp_file.name, sep='\t')
    
    # update CNV data argument
    args.cnv_norm = temp_file.name
    
    return args

def crisprseq_parseargs():
    """
    Parsing mageck arguments.
    """
    parser=argparse.ArgumentParser(description='mageck: performs sgRNA, gene and pathway analysis on CRISPR-Cas9 screening data.')
    # definition of sub commands
    subparser=parser.add_subparsers(help='commands to run mageck',dest='subcmd')
    
    # countonly: only collect counts
    arg_count(subparser)
    
    # MLE
    arg_mle(subparser)
    
    # args=parser.parse_args(['mle', 
    #                         '-k', '/group/iorio/Alessandro/CN_benchmark/leukemia.new.csv',
    #                         '-d', '/group/iorio/Alessandro/CN_benchmark/designmat.txt',
    #                         '-n', 'beta_leukemia',
    #                         '--cnv-norm', '/group/iorio/Alessandro/CN_benchmark/cnv_data.txt',
    #                         '--no-permutation-by-group'])
    
    # args=parser.parse_args(['mle', 
    #                         '-k', '/group/iorio/Alessandro/CN_benchmark/test_readcounts.csv',
    #                         '-s', '/group/iorio/Alessandro/CN_benchmark/data/ScreenSequenceMap.csv',
    #                         '-n', 'data/corrected/avana_supervised',
    #                         '--cnv-norm', '/group/iorio/Alessandro/CN_benchmark/data/OmicsCNGene.csv',
    #                         '--permutation-round', '10',
    #                         '--no-permutation-by-group'])
    
    # create design matrix if not provided
    if args.design_matrix == None and args.screen_sequence_map is not None:
        args = create_design_matrix(args)
    elif args.design_matrix == None and args.screen_sequence_map is None:
        print('Please specify a design matrix if a screen to sequence map file is not provided.')
        sys.exit(-1)

    # format CNV data if provided
    if args.cnv_norm is not None:
        args = format_CNV_data(args)
    
    design_matrix_path = args.design_matrix
    cnv_norm_path = args.cnv_norm

    if args.subcmd == 'mle':
        mageckmle_main(parsedargs=args); # ignoring the script path, and the sub command
    else:
        parser.print_help()
        sys.exit(0)
    
    # delete temporary files
    if design_matrix_path is not None and 'tmp' in design_matrix_path:
        os.remove(design_matrix_path)
    
    if cnv_norm_path is not None and 'tmp' in cnv_norm_path:
        os.remove(cnv_norm_path)


if __name__ == '__main__':
    crisprseq_parseargs()

