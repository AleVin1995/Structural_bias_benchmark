from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import WhiteKernel, ConstantKernel, RBF
import argparse
import numpy as np
import pandas as pd


# Merge guide map and copy number data
def intersect_guide_cn(guide, cn_cell_line):
    ## convert coords to int
    cn_cell_line.loc[:, 'Start'] = cn_cell_line.loc[:, 'Start'].astype(int)
    cn_cell_line.loc[:, 'End'] = cn_cell_line.loc[:, 'End'].astype(int)

    guide.loc[:, 'Start_sgRNA'] = guide.loc[:, 'Start_sgRNA'].astype(int)
    guide.loc[:, 'sample'] = cn_cell_line.loc[:, 'sample'].unique()[0]

    ## add columns to guide
    guide['Start'] = np.nan
    guide['End'] = np.nan
    guide['ratio'] = np.nan

    for i in range(cn_cell_line.shape[0]):
        cn_seg_chr = cn_cell_line.loc[i, 'Chr']
        cn_seg_start = cn_cell_line.loc[i, 'Start']
        cn_seg_end = cn_cell_line.loc[i, 'End']
        cn_seg_ratio = cn_cell_line.loc[i, 'ratio']

        ## get guide indeces within each copy number segment
        guide_bool = [(guide['Chr'] == cn_seg_chr) &
                      (guide['Start_sgRNA'] >= cn_seg_start) & 
                      (guide['Start_sgRNA'] <= cn_seg_end)]
        guide_idx = np.where(np.array(guide_bool)[0])[0]

        if len(guide_idx) > 0:
            ## check which guide columns are not nan
            guide_not_nan = guide.loc[guide_idx, ['Start', 'End', 'ratio']].notna().any(axis=1)
            guide_not_nan_idx = np.array(guide_not_nan.index)[np.where(guide_not_nan)[0]]

            ## remove indeces that are not nan
            guide_idx = np.setdiff1d(guide_idx, guide_not_nan_idx)

            if len(guide_idx) > 0:
                guide.loc[guide_idx, 'Start'] = cn_seg_start
                guide.loc[guide_idx, 'End'] = cn_seg_end
                guide.loc[guide_idx, 'ratio'] = cn_seg_ratio

            ## duplicate rows that have multiple guides
            if len(guide_not_nan_idx) > 0:
                start_idx = guide.shape[0]
                end_idx = start_idx + len(guide_not_nan_idx)
                guide = pd.concat([guide, guide.loc[guide_not_nan_idx, :]], axis=0).reset_index(drop=True)

                guide.loc[start_idx:end_idx, 'Start'] = cn_seg_start
                guide.loc[start_idx:end_idx, 'End'] = cn_seg_end
                guide.loc[start_idx:end_idx, 'ratio'] = cn_seg_ratio
        
    guide = guide.dropna(subset=['Start', 'End', 'ratio']).reset_index(drop=True)
    return guide


# Correct sgRNA-level LFC
def correct_lfc(df, n_sgrna=10):
    df = df.dropna(subset=['Chr', 'Start', 'End', 'ratio', 'LFC'])
    df = df.copy()
    
    # Fit Gaussian Process Regression on segment LFC
    gpr = CrispyGaussian(df, n_sgrna=n_sgrna)
    gpr = gpr.fit(x='ratio', y='LFC')

    df.loc[:, "gp_mean"] = gpr.predict(df['ratio'])

    # Correct fold-change by subtracting the estimated mean
    df.loc[:, "corrected"] = df.eval("LFC - gp_mean")

    return df


class CrispyGaussian(GaussianProcessRegressor):
    SEGMENT_COLUMNS = ["Chr", "Start", "End"]

    SEGMENT_AGG_FUN = dict(
        LFC=np.mean,
        ratio=np.mean,
        sgRNA="count",
    )

    def __init__(
        self,
        bed_seg,
        kernel=None,
        n_sgrna=10,
        alpha=1e-10,
        n_restarts_optimizer=3,
        optimizer="fmin_l_bfgs_b",
        normalize_y=False,
        copy_X_train=False,
        random_state=None,
    ):
        
        self.bed_seg = bed_seg.groupby(self.SEGMENT_COLUMNS).agg(self.SEGMENT_AGG_FUN)
        self.n_sgrna = n_sgrna

        super().__init__(
            alpha=alpha,
            kernel=self.get_default_kernel() if kernel is None else kernel,
            optimizer=optimizer,
            normalize_y=normalize_y,
            random_state=random_state,
            copy_X_train=copy_X_train,
            n_restarts_optimizer=n_restarts_optimizer,
        )

    @property
    def _constructor(self):
        return CrispyGaussian

    @staticmethod
    def get_default_kernel():
        """
        Deafult Gaussian Process kernel

        :return: sklearn.gaussian_process.kernels.Sum
        """
        return ConstantKernel() * RBF() + WhiteKernel()

    def fit(self, x, y):
        x = self.bed_seg.query(f"sgRNA >= {self.n_sgrna}")[x]
        y = self.bed_seg.query(f"sgRNA >= {self.n_sgrna}")[y]

        if x.ndim == 1:
            x = x.values.reshape(-1, 1)
        
        if y.ndim == 1:
            y = y.values.reshape(-1, 1)

        return super().fit(x, y)

    def predict(self, x=None, return_std=False, return_cov=False):
        if x is None:
            x = self.bed_seg[["ratio"]]
        
        if x.ndim == 1:
            x = x.values.reshape(-1, 1)

        return super().predict(x, return_std=return_std, return_cov=return_cov)


def main():
    # create ArgParser object
    parser = argparse.ArgumentParser(description='Run Chronos')

    # add arguments
    parser.add_argument('--lfc', dest = 'lfc', type = str, required = True,
                        help = 'raw (sgRNA-level) LFC dataset')

    parser.add_argument('--cn', dest = 'cn', type = str, required = True,
                        help = 'Copy number segment dataset to use for correction')
    
    parser.add_argument('--map', dest = 'map', type = str, required = True,
                        help = 'Mapping file to map profile IDs to sample IDs')
    
    parser.add_argument('--lib', dest = 'lib', type = str, required = True,
                        help = 'Library file to map sgRNA IDs to gene IDs')
    
    parser.add_argument('-o', '--out', dest = 'out', type = str, required = True,
                        help = 'Path to save corrected file')
    
    args = parser.parse_args()

    # Load data
    sgrna_effects = pd.read_csv(args.lfc, index_col=0)
    guide_map = pd.read_csv(args.lib)
    cn = pd.read_csv(args.cn)
    profile = pd.read_csv(args.map)

    # Format guide map
    guide_map = guide_map[guide_map["UsedByChronos"] == True]
    guide_map.loc[:, 'Gene'] = guide_map.loc[:, 'Gene'].str.replace(r'\s\(.*\)', '', regex=True)
    guide_map = guide_map[['sgRNA', 'Gene', 'GenomeAlignment']]
    guide_map[['Chr', 'Start', 'Strand']] = guide_map.loc[:, 'GenomeAlignment'].str.split('_', expand=True)
    guide_map = guide_map.drop(columns=['Strand', 'GenomeAlignment'])
    guide_map = guide_map.rename(columns={'Start': 'Start_sgRNA'})
    guide_map = guide_map.reset_index(drop=True)

    # Format copy number data
    cn = pd.merge(cn, profile, how='left', on='ProfileID')
    cn = cn[["ModelID", "Chromosome", "Start", "End", "SegmentMean"]]
    cn = cn.rename(columns={'Chromosome': 'Chr', 'SegmentMean': 'ratio', 'ModelID': 'sample'})
    cn["Chr"] = "chr" + cn["Chr"]

    # Apply correction cell-wise
    for i, cell_line in enumerate(sgrna_effects.columns):
        print(f'Processing cell line: {cell_line} {i+1}/{sgrna_effects.shape[1]}')

        fc = sgrna_effects.iloc[:, i]
        cn_cell_line = cn[cn["sample"] == cell_line]
        cn_cell_line = cn_cell_line.reset_index(drop=True)

        if cn_cell_line.shape[0] > 0:
            ## Intersect guide map and copy number data
            guide_cn = intersect_guide_cn(guide_map, cn_cell_line)

            ## Intersect fold change and copy number data
            fc_cn = pd.merge(fc, guide_cn, how='left', on='sgRNA')
            fc_cn = fc_cn.rename(columns={fc_cn.columns.to_list()[1]: 'LFC'})

            ## Run Crispy
            corrected_fc = correct_lfc(fc_cn)
            corrected_fc = corrected_fc[['sgRNA', 'Gene', 'corrected']]
            corrected_fc = corrected_fc.rename(columns={'corrected': cell_line})
            corrected_fc = corrected_fc.drop_duplicates()

            ## convert to gene-level LFC
            corrected_fc = corrected_fc.drop(columns=['sgRNA'])
            corrected_fc = corrected_fc.groupby('Gene').agg('median')

            if i == 0:
                corrected_fcs = corrected_fc
            else:
                corrected_fcs = pd.merge(corrected_fcs, corrected_fc, how='outer', on=['Gene'])

    # Save results
    corrected_fcs.to_csv(args.out)


if __name__ == '__main__':
    main()