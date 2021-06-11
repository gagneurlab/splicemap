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
from count_table.dataclasses import Junction
from count_table.utils import get_variants_around_junction, \
    read_genes_from_gtf, remove_chr_from_chrom_annotation
from count_table.splice_map import SpliceMap
from count_table.dataclasses import EventPSI5, EventPSI3
# from mmsplice.utils import logit
# from mmsplice_scripts.data.utils import clip

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


class CountTable:
    required_columns = ('Chromosome', 'Start', 'End', 'Strand')
    # TODO: implement optional columns

    def __init__(self, df, name):
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
            f'First 4 columns need to be {CountTable.required_columns}'
        self.df = self.df.astype(self._dtype(self.samples))
        self._set_index()

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

    @staticmethod
    def validate(df):
        return CountTable._validate_columns(df.columns)

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
    def read_csv(cls, csv_path, name=None):
        '''
        Args:
          csv_path: cvs file containing 'Chromosome', 'Start', 'End', 'Strand'
            and sample counts as columns.
          name: name of count table will be stored in metadata \
            (containing dataset, tissue, version .etc)
        '''
        with open(csv_path) as f:
            columns = next(f).strip().split(',')
            assert CountTable._validate_columns(columns), \
                f'First 4 columns need to be {CountTable.required_columns}'
            samples = columns[4:]

        df = pd.read_csv(csv_path, dtype=cls._dtype(samples))

        if name is None:
            name = Path(csv_path).name.split('.')[0]

        return cls(df, name)

    @property
    def junctions(self):
        return self.df.index

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
            raise AttributeError('CountTable does not has annotation '
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
            raise ValueError('Chromosome annotation in fasta file does match'
                             ' with count table chromosome annotation.')
        else:
            self.df = self.df[self.df['Chromosome'].isin(chroms)]

        junctions = tqdm(self.junctions) if progress else self.junctions
        self.df['Strand'] = [
            infer_junction_strand(i, fasta, default=default)
            for i in junctions
        ]
        self._set_index()

    @property
    def splice_site5(self):
        if self._splice_site5 is None:
            df_pos = self.df[self.df['Strand'] == '+']
            splice_site_pos = df_pos['Chromosome'].astype(str) + ':' \
                + df_pos['Start'].astype(str) + ':' + \
                df_pos['Strand'].astype(str)

            df_neg = self.df[self.df['Strand'] == '-']
            splice_site_neg = df_neg['Chromosome'].astype(str) + ':' \
                + df_neg['End'].astype(str) + ':' + \
                df_neg['Strand'].astype(str)

            splice_site = pd.concat([splice_site_pos, splice_site_neg])
            self._splice_site5 = pd.DataFrame({
                'junctions': self.junctions,
                'splice_site': splice_site.loc[self.junctions]
            }).set_index('junctions')
        return self._splice_site5

    @property
    def splice_site3(self):
        if self._splice_site3 is None:
            df_pos = self.df[self.df['Strand'] == '+']
            splice_site_pos = df_pos['Chromosome'].astype(str) + ':' \
                + df_pos['End'].astype(str) + ':' + \
                df_pos['Strand'].astype(str)

            df_neg = self.df[self.df['Strand'] == '-']
            splice_site_neg = df_neg['Chromosome'].astype(str) + ':' \
                + df_neg['Start'].astype(str) + ':' + \
                df_neg['Strand'].astype(str)

            splice_site = pd.concat([splice_site_pos, splice_site_neg])
            self._splice_site3 = pd.DataFrame({
                'junctions': self.junctions,
                'splice_site': splice_site.loc[self.junctions]
            }).set_index('junctions')
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
        return CountTable(self.df.loc[junctions], name=self.name)

    def _filter_event(self, junctions, events):
        keep_events = set(events.loc[junctions]['events'])
        event_filter = events['events'].isin(keep_events)
        return CountTable(self.df.loc[event_filter], name=self.name)

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
        return CountTable(self.df[percentile_filter], name=self.name)

    def _median_filter_event_counts(self, event_counts, cutoff=1):
        return event_counts[event_counts.median(axis=1) >= cutoff]

    def _median_filter(self, event_counts, event, cutoff=1):
        expressed_events = self._median_filter_event_counts(
            event_counts, cutoff).index
        ct = CountTable(self.df.loc[event['events'].isin(
            expressed_events)], name=self.name)
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
        ct = CountTable(self.df.loc[event['events'].isin(
            expressed_events)], name=self.name)
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

    def _ref_psi_with_mean_std(self, counts, event_counts, event):
        count_rows = self._join_count_with_event_counts(
            counts, event_counts, event)
        k = count_rows[self.samples].values
        n = count_rows[self._event_samples].values
        psi = k / n
        return pd.DataFrame({
            'junctions': count_rows.index,
            'ref_psi': np.nanmean(psi, axis=1),
            'std': np.nanstd(psi, axis=1)
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
        elif method == 'mean':
            df = self._ref_psi_with_mean_std(
                self.counts, event_counts, event)
        else:
            raise ValueError('method name %s is valid' % method)

        df = self.junction_df.join(splice_site).join(event).join(df)

        if annotation:
            df = df.join(self.annotation)
            df = df[~df['gene_name'].isna()]
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

        df_gene_junc = df_gene_junc[[
            'junctions', 'gene_id', 'gene_name', 'gene_type'
        ]].groupby(by='junctions').agg(list)

        return df_gene_junc

    def _load_junction_from_gtf(self, gr_gtf):
        gr_junc = gr_gtf.features.introns(by='transcript')
        gr_junc = self._infer_gene_type(gr_junc)

        df_junc = gr_junc.df
        df_junc['Strand'] = df_junc['Strand'].astype(str)

        df_junc['junctions'] = df_to_interval_str(df_junc)

        cols_index = ['junctions', 'Chromosome', 'Start', 'End', 'Strand']
        cols_agg = ['gene_id', 'gene_name', 'transcript_id', 'gene_type']

        df_junc = df_junc[[*cols_index, *cols_agg]] \
            .groupby(cols_index).agg(lambda x: list(set(x)))

        return df_junc.reset_index().set_index('junctions')

    def infer_annotation(self, gtf_file):
        gr_gtf = pr.read_gtf(gtf_file)
        gr_gene = self._pr_genes_from_gtf(gr_gtf)

        df_gene_junc = self._gene_junction_overlap(gr_gene)
        df_gtf_junc = self._load_junction_from_gtf(gr_gtf)

        # df_gene_junc['weak'] = ~df_gene_junc.index.isin(df_gtf_junc.index)

        # # If not weak, use genes from gtf
        # for col in ['gene_id', 'gene_name', 'transcript_id', 'gene_type']:
        #     df_gene_junc.at[~df_gene_junc['weak'], col] = df_gtf_junc.loc[
        #         df_gene_junc[~df_gene_junc['weak']].index, col]

        # for col in ['gene_id', 'gene_name', 'transcript_id', 'gene_type']:
        #     df_gene_junc[col] = df_gene_junc[col].str.join(';')

        df_gene_junc['splice_site_psi5'] = df_gene_junc.index.map(
            lambda x: EventPSI5.from_str(x).donor_str)
        df_gene_junc['splice_site_psi3'] = df_gene_junc.index.map(
            lambda x: EventPSI3.from_str(x).acceptor_str)
        splice_sites5_gtf = set(df_gtf_junc.index.map(
            lambda x: EventPSI5.from_str(x).donor_str))
        splice_sites3_gtf = set(df_gtf_junc.index.map(
            lambda x: EventPSI3.from_str(x).acceptor_str))

        novel_junction = set(df_gene_junc.index).difference(set(df_gtf_junc.index))
        weak_site_donor = set(
            df_gene_junc['splice_site_psi5']).difference(splice_sites5_gtf)
        weak_site_acceptor = set(
            df_gene_junc['splice_site_psi3']).difference(splice_sites3_gtf)

        df_gene_junc['novel_junction'] = ~df_gene_junc.index.isin(df_gtf_junc.index)
        df_gene_junc['weak_site_donor'] = df_gene_junc['splice_site_psi5'].isin(
            weak_site_donor)
        df_gene_junc['weak_site_acceptor'] = df_gene_junc['splice_site_psi3'].isin(
            weak_site_acceptor)

        # If not novel_junction, use genes from gtf
        for col in ['gene_id', 'gene_name', 'transcript_id', 'gene_type']:
            df_gene_junc.at[~df_gene_junc['novel_junction'], col] = df_gtf_junc.loc[
                df_gene_junc[~df_gene_junc['novel_junction']].index, col]

        for col in ['gene_id', 'gene_name', 'transcript_id', 'gene_type']:
            df_gene_junc[col] = df_gene_junc[col].str.join(';')

        df_gene_junc = df_gene_junc.drop(
            columns={'splice_site_psi5', 'splice_site_psi3'})

        self._annotation = df_gene_junc
        return self._annotation

    def join(self, ct, suffix='_other'):
        df = self.df.join(ct.df, rsuffix=suffix, how='outer')

        for i in self.required_columns:
            other_col = f'{i}{suffix}'
            df[i] = np.where(~df[i].isna(), df[i], df[other_col])
            del df[other_col]

        return CountTable(df.fillna(0), name=self.name)

    # def _delta_logit_psi(self, psi, ref_psi, clip_threshold=0.01):
    #     ref_psi = clip(ref_psi['ref_psi'].values.reshape((-1, 1)),
    #                    threshold=clip_threshold),
    #     ref_psi = clip(psi.values, clip_threshold=clip_threshold)
    #     return logit(psi.values) - logit(ref_psi)

    # def delta_logit_psi5(self, method='k/n'):
    #     return self._delta_logit_psi(
    #         self.psi5, self.ref_psi5(method=method))

    # def delta_logit_psi3(self, method='k/n'):
    #     return self._delta_logit_psi(
    #         self.psi3, self.ref_psi3(method=method))
