import colorsys
from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm
import seaborn as sns
import pyranges as pr
import matplotlib.pyplot as plt
import matplotlib.colors as mc
from fastbetabino import fit_alpha_beta
from sklearn.mixture import GaussianMixture
from kipoiseq.extractors import FastaStringExtractor
from splicemap.dataclasses import Junction
from splicemap.utils import get_variants_around_junction, \
    remove_chr_from_chrom_annotation
from splicemap.splice_map import SpliceMap


gt_mapping = {0: 'AA',  1: 'Aa', 3: 'aa',  2: 'NN'}


dinucleotide_motif = {
    ('GT', 'AG'): '+', ('CT', 'AC'): '-',
    ('GC', 'AG'): '+', ('CT', 'GC'): '-',
    ('AT', 'AC'): '+', ('GT', 'AT'): '-',
    ('AT', 'AA'): '+', ('TT', 'AT'): '-',
    ('AT', 'AG'): '+', ('CT', 'AT'): '-'
}

dinucleotide_start_motif = {
    start: strand
    for (start, end), strand in dinucleotide_motif.items()
}

dinucleotide_end_motif = {
    end: strand
    for (start, end), strand in dinucleotide_motif.items()
}


def infer_junction_strand(junction, fasta, default='.'):
    junction = Junction.from_str(junction) if type(
        junction) == str else junction
    dinucleotide = junction.dinucleotide_region()
    motif = (fasta.extract(dinucleotide[0]), fasta.extract(dinucleotide[1]))

    if motif in dinucleotide_motif:
        return dinucleotide_motif[motif]
    elif motif[0] == 'GT' or motif[1] == 'AG':
        return '+'
    elif motif[0] == 'CT' or motif[1] == 'AC':
        return '-'
    elif motif[0] in dinucleotide_start_motif:
        return dinucleotide_start_motif[motif[0]]
    elif motif[1] in dinucleotide_end_motif:
        return dinucleotide_end_motif[motif[1]]
    else:
        return default


def df_to_interval_str(df):
    return df['Chromosome'].astype(str) + ':' + df['Start'].astype(str) + '-' \
        + df['End'].astype(str) + ':' + df['Strand'].astype(str)


class SpliceCountTable:
    required_columns = ('Chromosome', 'Start', 'End', 'Strand')
    # TODO: implement optional columns

    def __init__(self, df, name, gene_expression=None):
        '''
        Args:
          df: pd.DataFrame containing 'Chromosome', 'Start', 'End', 'Strand'
            and sample counts as columns.
          name: name of count table will be stored in metadata \
            (containing dataset, tissue, version .etc)
        '''
        self.df = df
        self.name = name
        assert self.validate(df), \
            f'First 4 columns need to be {SpliceCountTable.required_columns}'
        self.df = self.df.astype(self._dtype(self.samples))
        self._set_index()

        self._gene_expression = None or gene_expression
        if gene_expression is not None:
            self.validate_expression()

    def _set_index(self):
        df = self.df.reset_index(drop=True)
        df['junctions'] = df_to_interval_str(df)
        self.df = df.set_index('junctions')
        self._splice_site5 = None
        self._splice_site3 = None
        self._event5 = None
        self._event3 = None
        self._event5_counts = None
        self._event3_counts = None
        self._psi5 = None
        self._psi3 = None
        self._annotation = None
        self._annotation_exploded = None
        self._gene_expression_median = None

    @staticmethod
    def validate(df):
        return SpliceCountTable._validate_columns(df.columns)

    def validate_expression(self):
        if self._gene_expression.index.name != 'gene_id':
            if 'gene_id' in self._gene_expression.columns:
                self._gene_expression = self._gene_expression.set_index(
                    'gene_id')
            else:
                raise ValueError('`gene_id` is not index or columns')

        if not self._gene_expression.index.is_unique:
            raise ValueError(
                '`gene_id` in the expression file need to be unique!')

    @staticmethod
    def _validate_columns(columns):
        return tuple(columns[:4]) == ('Chromosome', 'Start', 'End', 'Strand')

    @staticmethod
    def _dtype(samples):
        dtype = {
            'Chromosome': 'category',
            'Start': 'int32',
            'End': 'int32',
            'Strand': 'category'
        }
        for i in samples:
            dtype[i] = 'int32'
        return dtype

    @classmethod
    def read_csv(cls, csv_path, name=None, gene_expression_path=None):
        '''
        Args:
          csv_path: cvs file containing 'Chromosome', 'Start', 'End', 'Strand'
            and sample counts as columns.
          name: name of count table will be stored in metadata \
            (containing dataset, tissue, version .etc)
        '''
        with open(csv_path) as f:
            columns = next(f).strip().split(',')
            assert SpliceCountTable._validate_columns(columns), \
                f'First 4 columns need to be {SpliceCountTable.required_columns}'
            samples = columns[4:]

        df = pd.read_csv(csv_path, dtype=cls._dtype(samples))

        if name is None:
            name = Path(csv_path).name.split('.')[0]

        if gene_expression_path is not None:
            df_exp = pd.read_csv(gene_expression_path).set_index('gene_id')
        else:
            df_exp = None

        return cls(df, name, gene_expression=df_exp)

    @property
    def junctions(self):
        return self.df.index

    @property
    def gene_expression_median(self):
        if self._gene_expression_median is None:
            if type(self._gene_expression) is dict:
                self._gene_expression_median = pd.Series(
                    self._gene_expression, name='gene_tpm')
            elif type(self._gene_expression) is pd.Series:
                self._gene_expression_median = self._gene_expression
                self._gene_expression_median.name = 'gene_tpm'
                self._gene_expression_median.index.name = None
            elif type(self._gene_expression) is pd.DataFrame:
                self._gene_expression_median = self._gene_expression.median(
                    axis=1)
                self._gene_expression_median.name = 'gene_tpm'
                self._gene_expression_median.index.name = None
        return self._gene_expression_median

    @property
    def junction_df(self):
        return self.df[self.df.columns[:4]]

    @property
    def samples(self):
        return self.df.columns[4:].tolist()

    @property
    def annotation(self):
        if self._annotation is not None:
            return self._annotation
        else:
            raise AttributeError('SpliceCountTable does not have annotation '
                                 'unless `infer_annotation` is called.')

    def update_samples(self, mapping):
        self.df = self.df.rename(columns=mapping)

    @property
    def counts(self):
        return self.df[self.samples]

    def infer_strand(self, fasta_file, default='+', progress=False):
        fasta = FastaStringExtractor(fasta_file)
        chroms = FastaStringExtractor(fasta_file).fasta.keys()
        chroms = set(chroms).intersection(self.df['Chromosome'].unique())

        if len(chroms) == 0:
            raise ValueError('Chromosome annotation in fasta file does not match'
                             ' with count table chromosome annotation.')
        else:
            self.df = self.df[self.df['Chromosome'].isin(chroms)]

        junctions = tqdm(self.junctions) if progress else self.junctions
        self.df['Strand'] = [
            infer_junction_strand(i, fasta, default=default)
            for i in junctions
        ]
        self._set_index()

    @staticmethod
    def _splice_site5_from(df):
        df_pos = df[df['Strand'] == '+']
        splice_site_pos = df_pos['Chromosome'].astype(str) + ':' \
            + df_pos['Start'].astype(str) + ':' + \
            df_pos['Strand'].astype(str)

        df_neg = df[df['Strand'] == '-']
        splice_site_neg = df_neg['Chromosome'].astype(str) + ':' \
            + df_neg['End'].astype(str) + ':' + \
            df_neg['Strand'].astype(str)

        splice_site = pd.concat([splice_site_pos, splice_site_neg])
        return splice_site.to_frame().rename(columns={0: 'splice_site'})

    @property
    def splice_site5(self):
        if self._splice_site5 is None:
            self._splice_site5 = self._splice_site5_from(self.df)
        return self._splice_site5

    @staticmethod
    def _splice_site3_from(df):
        df_pos = df[df['Strand'] == '+']
        splice_site_pos = df_pos['Chromosome'].astype(str) + ':' \
            + df_pos['End'].astype(str) + ':' + \
            df_pos['Strand'].astype(str)

        df_neg = df[df['Strand'] == '-']
        splice_site_neg = df_neg['Chromosome'].astype(str) + ':' \
            + df_neg['Start'].astype(str) + ':' + \
            df_neg['Strand'].astype(str)

        splice_site = pd.concat([splice_site_pos, splice_site_neg])
        return splice_site.to_frame().rename(columns={0: 'splice_site'})

    @property
    def splice_site3(self):
        if self._splice_site3 is None:
            self._splice_site3 = self._splice_site3_from(self.df)
        return self._splice_site3

    def _event(self, splice_site):
        df = splice_site.reset_index().groupby('splice_site').agg(list)
        df = pd.DataFrame({
            'junctions': df['junctions'],
            'events': df['junctions']
        }).explode('junctions').set_index('junctions')
        df['events'] = df['events'].str.join(';')
        return df.loc[self.junctions]

    @property
    def event5(self):
        if self._event5 is None:
            self._event5 = self._event(self.splice_site5)
        return self._event5

    @property
    def event3(self):
        if self._event3 is None:
            self._event3 = self._event(self.splice_site3)
        return self._event3

    def _event_counts(self, event):
        return self.counts.join(event).groupby('events').sum()

    @property
    def event5_counts(self):
        if self._event5_counts is None:
            self._event5_counts = self._event_counts(self.event5)
        return self._event5_counts

    @property
    def event3_counts(self):
        if self._event3_counts is None:
            self._event3_counts = self._event_counts(self.event3)
        return self._event3_counts

    def k(self, junctions):
        return self.df.loc[junctions, self.samples]

    def n5(self, junctions):
        df = self.event5_counts.loc[self.event5.loc[junctions]['events']]
        if type(junctions) == str:
            return df
        df['junctions'] = junctions
        return df.set_index('junctions')

    def n3(self, junctions):
        df = self.event3_counts.loc[self.event3.loc[junctions]['events']]
        if type(junctions) == str:
            return df
        df['junctions'] = junctions
        return df.set_index('junctions')

    def _kn(self, junction, n):
        df = pd.DataFrame({
            'k': self.k(junction).astype('int64'),
            'n': n
        })
        df.index.name = 'sample'
        return df

    def kn5(self, junction):
        return self._kn(junction, self.n5(junction))

    def kn3(self, junction):
        return self._kn(junction, self.n3(junction))

    def _plot_kn(self, df, log=False, highlight=None,
                 highlight_name='outlier'):
        if highlight:
            df[highlight_name] = df.index.isin(highlight)

            ax = sns.scatterplot(
                x='n', y='k', data=df,
                hue=highlight_name, style=highlight_name,
                size=highlight_name, sizes=(200, 20))
        else:
            ax = sns.scatterplot(x='n', y='k', data=df)

        if log:
            ax.set(xscale='log', yscale='log')
        return ax

    def plot_kn5(self, junction, log=False, highlight=None,
                 highlight_name='outlier'):
        return self._plot_kn(self.kn5(junction), log=log, highlight=highlight,
                             highlight_name=highlight_name)

    def plot_kn3(self, junction, log=False, highlight=None,
                 highlight_name='outlier'):
        return self._plot_kn(self.kn3(junction), log=log, highlight=highlight,
                             highlight_name=highlight_name)

    def _lighten_color(self, color, amount=0.5):
        try:
            c = mc.cnames[color]
        except:
            c = color
        c = colorsys.rgb_to_hls(*mc.to_rgb(c))
        return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

    def _plot_psi(self, df, x=None, hue=None, min_read=1, swarm=False,
                  plot_type='boxplot', ylim=True):
        df = df[df['n'] >= min_read]
        if df.shape[0] == 0:
            raise ValueError('psi is not defined for all the junctions'
                             'because of the low n')
        df['ref_psi'] = df['k'] / df['n']

        if ylim:
            plt.ylim((-0.1, 1.1))

        if plot_type == 'boxplot':
            plot = sns.boxplot
        elif plot_type == 'violinplot':
            plot = sns.violinplot
        else:
            raise ValueError('Plot type is not support. '
                             'Supported types are "violinplot" and "boxplot"')

        kwargs = {}
        if x:
            kwargs['x'] = x
        if hue:
            kwargs['hue'] = hue
            kwargs['dodge'] = True
            kwargs['hue_order'] = ['AA', 'Aa', 'aa', 'NN']
        ax = plot(y='ref_psi', data=df, **kwargs)

        if plot_type == 'boxplot':
            for i, artist in enumerate(ax.artists):
                # Set the linecolor on the artist to the facecolor, and set the facecolor to None
                col = self._lighten_color(artist.get_facecolor(), 1.2)
                artist.set_edgecolor(col)

                # Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
                # Loop over them here, and use the same colour as above
                for j in range(i*6, i*6+6):
                    line = ax.lines[j]
                    line.set_color(col)
                    line.set_mfc(col)
                    line.set_mec(col)
                    line.set_linewidth(0.5)

        if swarm:
            sns.swarmplot(y='ref_psi', data=df,
                          color=".25", alpha=0.5, **kwargs)
        return ax

    def _plot_psi_samples(self, df, samples=None, category_name='outliers',
                          min_read=1, swarm=False, plot_type='boxplot',
                          ylim=True):
        if samples:
            df[category_name] = df.index.isin(samples)
        else:
            category_name = None
        return self._plot_psi(df, x=category_name, min_read=min_read,
                              swarm=swarm, plot_type=plot_type, ylim=ylim)

    def plot_psi5(self, junction, samples=None, category_name='outliers',
                  min_read=1, swarm=False, plot_type='boxplot'):
        return self._plot_psi_samples(
            self.kn5(junction), samples=samples, category_name=category_name,
            min_read=min_read, swarm=swarm, plot_type=plot_type)

    def plot_psi3(self, junction, samples=None, category_name='outliers',
                  min_read=1, swarm=False, plot_type='boxplot'):
        return self._plot_psi_samples(
            self.kn3(junction), samples=samples, category_name=category_name,
            min_read=min_read, swarm=swarm, plot_type=plot_type)

    def _psi_variants(self, df, vcf, variants):
        dfs = list()
        for v in variants:
            _df = df.copy()
            v_samples = dict(zip(vcf.samples, v.source.gt_types))
            genotype = list()
            for sample in df.index:
                if sample in vcf.samples:
                    genotype.append(gt_mapping[v_samples[sample]])
                else:
                    genotype.append('NN')
            _df['genotype'] = genotype
            _df['variant'] = str(v)
            dfs.append(_df)
        return pd.concat(dfs, axis=0)

    def _plot_psi_variants(self, df, vcf, variants, min_read=1, swarm=False,
                           plot_type='boxplot', ylim=True):
        df = self._psi_variants(df, vcf, variants)
        return self._plot_psi(df, x='variant', hue='genotype',
                              min_read=min_read, swarm=swarm,
                              plot_type=plot_type, ylim=ylim)

    def plot_psi5_variants(self, junction, vcf, variants=None, min_read=1,
                           swarm=False, plot_type='boxplot',
                           overhang=(100, 100), variant_filter=None):
        if not variants:
            donor_variants, acceptor_variants = get_variants_around_junction(
                vcf, junction, overhang=overhang,
                variant_filter=variant_filter)
            variants = acceptor_variants
        else:
            variants = [vcf.get_variant(v) for v in variants]
        return self._plot_psi_variants(
            self.kn5(junction), vcf, variants=variants, min_read=min_read,
            swarm=swarm, plot_type=plot_type)

    def plot_psi3_variants(self, junction, vcf, variants=None, min_read=1,
                           swarm=False, plot_type='boxplot',
                           overhang=(100, 100), variant_filter=None):
        if not variants:
            donor_variants, acceptor_variants = get_variants_around_junction(
                vcf, junction, overhang=overhang,
                variant_filter=variant_filter)
            variants = donor_variants
        else:
            variants = [vcf.get_variant(v) for v in variants]
        return self._plot_psi_variants(
            self.kn3(junction), vcf, variants=variants, min_read=min_read,
            swarm=swarm, plot_type=plot_type)

    def filter(self, junctions):
        return SpliceCountTable(self.df.loc[junctions], name=self.name,
                                gene_expression=self._gene_expression)

    def _filter_event(self, junctions, events):
        keep_events = set(events.loc[junctions]['events'])
        event_filter = events['events'].isin(keep_events)
        return SpliceCountTable(self.df.loc[event_filter], name=self.name,
                                gene_expression=self._gene_expression)

    def filter_event5(self, junctions):
        return self._filter_event(junctions, self.event5)

    def filter_event3(self, junctions):
        return self._filter_event(junctions, self.event3)

    def _plot_median_read_hist_event(self, event_counts):
        x = np.median(event_counts, axis=1)
        x = x[x > 0]
        plt.xlabel('log(N) in median sample')
        plt.ylabel('Number of junctions')
        plt.hist(np.log(x))

    def plot_median_read_hist_event5(self):
        self._plot_median_read_hist_event(self.event5_counts)

    def plot_median_read_hist_event3(self):
        self._plot_median_read_hist_event(self.event3_counts)

    def quantile_filter(self, quantile=95, min_read=1):
        '''
        Filter junction which observed at least %5 of samples
        '''
        percentile_filter = np.percentile(
            self.counts, quantile, axis=1) >= min_read
        return SpliceCountTable(self.df[percentile_filter], name=self.name,
                                gene_expression=self._gene_expression)

    def _median_filter_event_counts(self, event_counts, cutoff=1):
        return event_counts[event_counts.median(axis=1) >= cutoff]

    def _median_filter(self, event_counts, event, cutoff=1):
        expressed_events = self._median_filter_event_counts(
            event_counts, cutoff).index
        ct = SpliceCountTable(
            self.df.loc[event['events'].isin(expressed_events)],
            name=self.name,
            gene_expression=self._gene_expression)
        return ct

    def event5_median_filter(self, cutoff=1):
        return self._median_filter(
            self.event5_counts, self.event5, cutoff)

    def event3_median_filter(self, cutoff=1):
        return self._median_filter(
            self.event3_counts, self.event3, cutoff)

    def _is_expressed_events(self, event_counts):
        mat = np.log(event_counts.median(axis=1)).values.reshape((-1, 1))
        gmm = GaussianMixture(n_components=2).fit(mat)
        labels = gmm.predict(mat)
        is_expressed = labels == labels[mat.argmax()]
        cutoff = mat[~is_expressed].max()
        return is_expressed, np.e**cutoff

    def _event_count_filter(self, event_counts, event):
        event_counts = self._median_filter_event_counts(event_counts)
        is_expressed, cutoff = self._is_expressed_events(event_counts)
        expressed_events = event_counts[is_expressed].index
        ct = SpliceCountTable(
            self.df.loc[event['events'].isin(expressed_events)],
            name=self.name,
            gene_expression=self._gene_expression)
        return ct, cutoff

    def event5_count_filter(self):
        return self._event_count_filter(
            self.event5_counts, self.event5)

    def event3_count_filter(self):
        return self._event_count_filter(
            self.event3_counts, self.event3)

    def _join_count_with_event_counts(self, counts, event_counts, event):
        return counts.join(event) \
                     .join(event_counts, on='events', rsuffix='_event')

    @property
    def _event_samples(self):
        return ['%s_event' % i for i in self.samples]

    def _ref_psi_with_beta_binomial(self, counts, event_counts,
                                    event, progress=False, niter=1000):
        count_rows = self._join_count_with_event_counts(
            counts, event_counts, event).iterrows()
        if progress:
            count_rows = tqdm(count_rows, total=counts.shape[0])

        for junc, row in count_rows:
            k = row[self.samples].tolist()
            n = row[self._event_samples].tolist()
            sum_k = np.sum(k)
            sum_n = np.sum(n)
            median_n = np.median(n)
            if ';' in row['events']:
                alpha, beta = fit_alpha_beta(n, k, niter=niter)
            else:
                alpha, beta = 100, 0.01
            psi = alpha / (alpha + beta)
            yield junc, psi, alpha, beta, sum_k, sum_n, median_n

    def _ref_psi_with_kn(self, counts, event_counts, event):
        count_rows = self._join_count_with_event_counts(
            counts, event_counts, event)
        k = count_rows[self.samples].sum(axis=1)
        n = count_rows[self._event_samples].sum(axis=1)
        median_n = count_rows[self._event_samples].median(axis=1)
        return pd.DataFrame({
            'junctions': count_rows.index,
            'ref_psi': k / n,
            'k': k,
            'n': n,
            'median_n': median_n
        }).set_index('junctions')

    def _ref_psi(self, event_counts, event, splice_site, method, annotation=True):
        if method == 'beta_binomial':
            df = pd.DataFrame([
                i for i in self._ref_psi_with_beta_binomial(
                    self.counts, event_counts, event)
            ], columns=['junctions', 'ref_psi', 'alpha', 'beta', 'k', 'n', 'median_n']) \
                .set_index('junctions')
        elif method == 'k/n':
            df = self._ref_psi_with_kn(
                self.counts, event_counts, event)
        else:
            raise ValueError('method name %s is valid' % method)

        df = self.junction_df.join(splice_site).join(event).join(df)

        if annotation:
            df = df.join(self.annotation)
            df = df[~df.index.get_level_values('gene_id').isna()]
            # df = df[~df['gene_name'].isna()]

        if self.gene_expression_median is not None:
            df = df.join(self.gene_expression_median, on='gene_id')

        return df

    def ref_psi5(self, method='k/n', annotation=True):
        return SpliceMap(self._ref_psi(self.event5_counts, self.event5, self.splice_site5,
                                       method=method, annotation=annotation),
                         name=self.name)

    def ref_psi3(self, method='k/n', annotation=True):
        return SpliceMap(self._ref_psi(self.event3_counts, self.event3, self.splice_site3,
                                       method=method, annotation=annotation),
                         name=self.name)

    def _psi(self, event_counts, event):
        count_rows = self._join_count_with_event_counts(
            self.counts, event_counts, event)
        k = count_rows[self.samples].values
        n = count_rows[self._event_samples].values
        return pd.DataFrame(k / n, index=self.junctions, columns=self.samples)

    @property
    def psi5(self):
        if self._psi5 is None:
            self._psi5 = self._psi(self.event5_counts, self.event5)
        return self._psi5

    @property
    def psi3(self):
        if self._psi3 is None:
            self._psi3 = self._psi(self.event3_counts, self.event3)
        return self._psi3

    def to_csv(self, path):
        self.df.to_csv(path, index=False)

    def junction_report(self, junction_id):
        # TODO: beta, alpha, pval, and delta_psi
        df = pd.DataFrame({
            'sample': self.samples,
            'count': self.counts.loc[junction_id, self.samples],
            'psi5': self.psi5.loc[junction_id, self.samples],
            'psi3': self.psi3.loc[junction_id, self.samples],
        }).set_index('sample')
        # if styled:
        #     df = df.style \
        #            .applymap(lambda x: 'color: red' if x <= 0.05 else '',
        #                      subset=['pval_psi5', 'pval_psi3']) \
        #         .bar(subset=['psi5', 'psi3'],
        #              color='#FFA07A', vmin=0, vmax=1) \
        #         .bar(subset=['expected_psi5', 'expected_psi3'],
        #              color='#ee1f5f', vmin=0, vmax=1)
        return df

    def _pr_genes_from_gtf(self, gr_gtf):
        pr_genes = gr_gtf.subset(lambda df: df['Feature'] == 'gene')

        if not any('chr' in i for i in self.df['Chromosome'].unique()):
            pr_genes = remove_chr_from_chrom_annotation(pr_genes)

        pr_genes = self._infer_gene_type(pr_genes)
        return pr_genes

    @staticmethod
    def _infer_gene_type(gr):
        # Aggregate junctions which can be mapped to multiple genes.
        if 'gene_type' in gr.columns:
            pass
        elif 'gene_biotype' in gr.columns:
            gr.gene_type = gr.gene_biotype
        else:
            raise ValueError('gene_type can not be inferred from gtf file')
        return gr

    def _gene_junction_overlap(self, pr_genes):
        pr_junctions = pr.PyRanges(self.junction_df.reset_index())
        # Overlap genes and junctions
        df_gene_junc = pr_genes.join(pr_junctions).df
        # Filter inter-genenic junctions
        df_gene_junc = df_gene_junc[
            (df_gene_junc['Start'] < df_gene_junc['Start_b'])
            & (df_gene_junc['End'] > df_gene_junc['End_b'])
        ]

        if 'Strand' not in df_gene_junc.columns:
            raise ValueError(
                'Strand is missing. You have to call "infer_strand()" first.')
        df_gene_junc = df_gene_junc[
            (df_gene_junc['Strand'].astype(str) ==
             df_gene_junc['Strand_b'].astype(str))
        ]

        df_gene_junc = df_gene_junc[[
            'junctions', 'gene_id', 'gene_name', 'gene_type'
        ]].drop_duplicates()

        return df_gene_junc.set_index('junctions')

    def _load_junction_from_gtf(self, gr_gtf):
        gr_gtf_junc = gr_gtf.features.introns(by='transcript')
        gr_gtf_junc = self._infer_gene_type(gr_gtf_junc)
        if not any('chr' in i for i in self.df['Chromosome'].unique()):
            gr_gtf_junc = remove_chr_from_chrom_annotation(gr_gtf_junc)

        df_gtf_junc = gr_gtf_junc.df
        df_gtf_junc['Chromosome'] = df_gtf_junc['Chromosome'].astype(str)
        df_gtf_junc['Strand'] = df_gtf_junc['Strand'].astype(str)

        df_gtf_junc['junctions'] = df_to_interval_str(df_gtf_junc)

        cols_index = ['junctions', 'Chromosome', 'Start', 'End', 'Strand',
                      'gene_id', 'gene_name', 'gene_type']
        cols_agg = ['transcript_id']

        df_gtf_junc = df_gtf_junc[[*cols_index, *cols_agg]] \
            .groupby(cols_index).agg(lambda x: ';'.join(set(x)))

        return df_gtf_junc.reset_index().set_index('junctions').drop_duplicates()

    def _infer_weak_novel(self, df_gene_junc, df_gtf_junc):
        ss5 = self.splice_site5
        df_gene_junc['splice_site5'] = ss5['splice_site']
        ss3 = self.splice_site3
        df_gene_junc['splice_site3'] = ss3['splice_site']

        ss5_gtf = self._splice_site5_from(df_gtf_junc)
        ss3_gtf = self._splice_site3_from(df_gtf_junc)

        novel_junction = set(df_gene_junc.index).difference(
            set(df_gtf_junc.index))
        weak_site_donor = set(ss5['splice_site']).difference(
            ss5_gtf['splice_site'])
        weak_site_acceptor = set(ss3['splice_site']).difference(
            ss3_gtf['splice_site'])

        df_gene_junc['novel_junction'] = df_gene_junc.index.isin(
            novel_junction)
        df_gene_junc['weak_site_donor'] = df_gene_junc['splice_site5'].isin(
            weak_site_donor)
        df_gene_junc['weak_site_acceptor'] = df_gene_junc['splice_site3'].isin(
            weak_site_acceptor)

        del df_gene_junc['splice_site5']
        del df_gene_junc['splice_site3']
        return df_gene_junc

    def _main_gene_id(self, gr_gtf):
        df_gtf = gr_gtf.df
        df_gtf['gene_id'] = df_gtf['gene_id'].apply(lambda x: x.split('.')[0])
        return pr.PyRanges(df_gtf)

    # TODO: add blacklist for regions that are enriched for splicing outliers
    def infer_annotation(self, gtf_file, protein_coding=False, main_gene_id=True):
        gr_gtf = pr.read_gtf(gtf_file)

        if protein_coding:
            gr_gtf = self._infer_gene_type(gr_gtf)
            gr_gtf = gr_gtf.subset(
                lambda df: df['gene_type'] == 'protein_coding')

        if main_gene_id:
            gr_gtf = self._main_gene_id(gr_gtf)

        gr_gene = self._pr_genes_from_gtf(gr_gtf)

        df_gene_junc = self._gene_junction_overlap(gr_gene)
        df_gtf_junc = self._load_junction_from_gtf(gr_gtf)

        df_gene_junc = self._infer_weak_novel(df_gene_junc, df_gtf_junc)

        # If not novel_junction, use genes from gtf
        df_gene_junc = df_gene_junc.set_index('gene_id', append=True)
        df_gtf_junc = df_gtf_junc.set_index('gene_id', append=True)

        gene_junc = df_gene_junc.index
        gtf_junc_gene = df_gtf_junc.index

        df_gene_junc = df_gene_junc[
            df_gene_junc['novel_junction']
            | gene_junc.isin(gtf_junc_gene)
        ]

        col = 'transcript_id'
        not_novel = ~df_gene_junc['novel_junction']
        df_gene_junc.loc[not_novel, col] = df_gtf_junc.loc[
            df_gene_junc[not_novel].index, col]

        self._annotation = df_gene_junc
        return self._annotation

    def join(self, ct, suffix='_other'):
        df = self.df.join(ct.df, rsuffix=suffix, how='outer')

        for i in self.required_columns:
            other_col = f'{i}{suffix}'
            df[i] = np.where(~df[i].isna(), df[i], df[other_col])
            del df[other_col]

        return SpliceCountTable(df.fillna(0), name=self.name,
                                gene_expression=self._gene_expression)
