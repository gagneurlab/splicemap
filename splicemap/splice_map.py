import pandas as pd
import re


class SpliceMap:

    def __init__(self, df, name):
        self.df = df
        self.name = name
        # self.method = self._infer_method(self.ref_tables)
        self.method = self._infer_method(self.df)

    @classmethod
    def read_csv(cls, path, **kwargs):
        with open(path) as f:
            line = f.readline()
            assert line.startswith('# name: '), \
                'Name field not defined in metadata'

            # TODO: split by first :
            name = line.split(':')[1].strip()

            import pdb
            pdb.set_trace()

            return cls(pd.read_csv(f, **kwargs), name)

    # def save_combined_ref_tables(self, save_path):
    #     self.df.reset_index().to_csv(save_path, index=False)

    # @staticmethod
    # def _infer_method(ref_tables):
    #     method = list()
    #     for df in ref_tables:
    #         if 'k' in df.columns and 'n' in df.columns:
    #             method.append('kn')
    #         elif 'alpha' in df.columns and 'beta' in df.columns:
    #             method.append('bb')
    #     return method

    @staticmethod
    def _infer_method(df):
        if 'k' in df.columns and 'n' in df.columns:
            method = 'kn'
        elif 'alpha' in df.columns and 'beta' in df.columns:
            method = 'bb'
        elif 'std' in df.columns:
            method = 'mean'
        return method

    def to_csv(path):
        with open(path, 'w') as f:
            f.write(f'# name: {self.name}')
        df.to_csv(path, mode='a')

    @property
    def junctions(self):
        return self.df.index

    @staticmethod
    def valid():
        raise NotImplementedError()

    @staticmethod
    def download(tissue_name):
        raise NotImplementedError()

    @staticmethod
    def fetch(tissue_name):
        raise NotImplementedError()


# df_gtf_junc = gr_gtf.features.introns(by='transcript').df
