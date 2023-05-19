import argparse
import chronos
import pandas as pd


def main():
    # create ArgParser object
    parser = argparse.ArgumentParser(description='Run Chronos')

    # add arguments
    parser.add_argument('--lfc', dest = 'lfc', type = str, required = True,
                        help = 'LFC dataset to correct')

    parser.add_argument('--cn', dest = 'cn', type = str, required = True,
                        help = 'Copy number dataset to use for correction')
    
    parser.add_argument('-o', '--out', dest = 'out', type = str, required = True,
                        help = 'Path to save corrected file')
    
    args = parser.parse_args()

    # Load data
    gene_effects = pd.read_csv(args.lfc, index_col=0).T

    cn = pd.read_csv(args.cn, index_col=0)
    cn.columns = cn.columns.str.replace(r'\s\(.*\)', '', regex=True)

    common_genes = list(set(gene_effects.columns) & set(cn.columns))
    cn = cn[common_genes]

    ## For genes not in CN data, assume no effect
    for col in set(gene_effects.columns) - set(cn.columns):
        cn = pd.concat([cn, pd.DataFrame({col: [1]*cn.shape[0]}, index=cn.index)], axis=1)

    # Run Chronos
    corrected, _ = chronos.alternate_CN(gene_effects, cn)

    # Save results
    corrected.T.to_csv(args.out)


if __name__ == '__main__':
    main()