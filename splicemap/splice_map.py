import pandas as pd
import gzip


class SpliceMap:

    def __init__(self, df, name):
        self.df = df
        self.name = name
        self.method = self._infer_method(self.df)

    @classmethod
    def read_csv(cls, path, **kwargs):
        with gzip.open(path, 'rt') as f:
            line = f.readline()
            assert line.startswith('# name: '), \
                'Name field not defined in metadata'

            # TODO: split by first :
            name = line.split(':')[1].strip()
            return cls(pd.read_csv(f, **kwargs), name)

    @staticmethod
    def _infer_method(df):
        if 'k' in df.columns and 'n' in df.columns:
            method = 'kn'
        elif 'alpha' in df.columns and 'beta' in df.columns:
            method = 'bb'
        elif 'std' in df.columns:
            method = 'mean'
        return method

    def to_csv(self, path):
        with gzip.open(path, 'wt') as f:
            f.write(f'# name: {self.name}\n')
        df = self.df.copy()
        if not isinstance(df.index, pd.RangeIndex):
            df = df.reset_index()
        df.to_csv(path, mode='a', index=False)

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
