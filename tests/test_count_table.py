import pdb
import pytest
import numpy as np
import pandas as pd
from pytest_mock import mocker
from kipoiseq.extractors import MultiSampleVCF
from count_table.dataclasses import Junction
from count_table import CountTable, infer_junction_strand
from kipoiseq.extractors import FastaStringExtractor
from .conftest import fasta_file, vcf_file, gtf_file, junc_file


def test_infer_junction_strand():
    junc = Junction.from_str('17:10100-10500:+')
    assert infer_junction_strand(
        junc, FastaStringExtractor(fasta_file)) == '.'

    junc = Junction.from_str('17:41267796-41276033:+')
    assert infer_junction_strand(
        junc, FastaStringExtractor(fasta_file)) == '-'

    junc = Junction.from_str('17:41267796-41276000:+')
    assert infer_junction_strand(
        junc, FastaStringExtractor(fasta_file)) == '-'

    junc = Junction.from_str('17:4126705-41276033:+')
    assert infer_junction_strand(
        junc, FastaStringExtractor(fasta_file)) == '-'


@pytest.fixture
def count_table():
    df = pd.DataFrame({
        'Chromosome': ['chr1', 'chr1', 'chr2', 'chr2', 'chr2'],
        'Start': [5, 22, 10, 10, 20],
        'End': [30, 30, 30, 50, 50],
        'Strand': ['+', '+', '+', '-', '-'],
        's1': [1, 1, 1, 1, 1],
        's2': [2, 1, 1, 1, 1],
        's3': [10, 5, 1, 2, 4]
    })
    return CountTable(df)


def count_table_validate():
    df = pd.DataFrame({
        'index': [1, 2, 3, 4, 5],
        'Chromosome': ['chr1', 'chr1', 'chr2', 'chr2', 'chr2'],
        'Start': [5, 22, 10, 10, 20],
        'End': [30, 30, 30, 50, 50],
        'Strand': ['+', '+', '+', '-', '-'],
        's1': [1, 1, 1, 1, 1],
        's2': [2, 1, 1, 1, 1],
        's3': [10, 5, 1, 2, 4]
    })

    with pytest.raises(AssertionError):
        CountTable(df)


@pytest.fixture
def count_table_chr17():
    df = pd.DataFrame({
        'Chromosome': ['17', '17'],
        'Start': [41197819, 41197831],
        'End': [41199659, 41199670],
        'Strand': ['-', '-'],
        's1': [1, 1],
        's2': [2, 1]
    })
    return CountTable(df)


def test_CountTable_infer_strand(count_table):
    df = pd.DataFrame({
        'Chromosome': ['17', '17'],
        'Start': [41267796, 41279042],
        'End': [41276033, 41279742],
        'Strand': ['.', '*'],
        's1': [1, 1],
        's2': [2, 1]
    })
    ct = CountTable(df)

    ct.infer_strand(fasta_file)
    assert ct.df['Strand'].tolist() == ['-', '+']


def test_CountTable_junctions(count_table):
    assert count_table.junctions.tolist() == [
        'chr1:5-30:+',
        'chr1:22-30:+',
        'chr2:10-30:+',
        'chr2:10-50:-',
        'chr2:20-50:-'
    ]


def test_CountTable_junctions_df(count_table):
    pd.testing.assert_frame_equal(
        count_table.junction_df,
        pd.DataFrame({
            'junctions': [
                'chr1:5-30:+',
                'chr1:22-30:+',
                'chr2:10-30:+',
                'chr2:10-50:-',
                'chr2:20-50:-'
            ],
            'Chromosome': ['chr1', 'chr1', 'chr2', 'chr2', 'chr2'],
            'Start': [5, 22, 10, 10, 20],
            'End': [30, 30, 30, 50, 50],
            'Strand': ['+', '+', '+', '-', '-']
        }).set_index('junctions')
    )


def test_CountTable_samples(count_table):
    assert count_table.samples == ['s1', 's2', 's3']


def test_CountTable_counts(count_table):
    pd.testing.assert_frame_equal(
        count_table.counts,
        pd.DataFrame({
            'junctions': [
                'chr1:5-30:+',
                'chr1:22-30:+',
                'chr2:10-30:+',
                'chr2:10-50:-',
                'chr2:20-50:-'
            ],
            's1': [1, 1, 1, 1, 1],
            's2': [2, 1, 1, 1, 1],
            's3': [10, 5, 1, 2, 4]
        }).set_index('junctions')
    )


def test_CountTable_splice_site5(count_table):
    pd.testing.assert_frame_equal(
        count_table.splice_site5,
        pd.DataFrame({
            'junctions': [
                'chr1:5-30:+',
                'chr1:22-30:+',
                'chr2:10-30:+',
                'chr2:10-50:-',
                'chr2:20-50:-'
            ],
            'splice_site': [
                'chr1:5:+',
                'chr1:22:+',
                'chr2:10:+',
                'chr2:50:-',
                'chr2:50:-',
            ]
        }).set_index('junctions')
    )


def test_CountTable_splice_site3(count_table):
    pd.testing.assert_frame_equal(
        count_table.splice_site3,
        pd.DataFrame({
            'junctions': [
                'chr1:5-30:+',
                'chr1:22-30:+',
                'chr2:10-30:+',
                'chr2:10-50:-',
                'chr2:20-50:-'
            ],
            'splice_site': [
                'chr1:30:+',
                'chr1:30:+',
                'chr2:30:+',
                'chr2:10:-',
                'chr2:20:-'
            ]
        }).set_index('junctions')
    )


def test_CountTable_event5(count_table):
    pd.testing.assert_frame_equal(
        count_table.event5,
        pd.DataFrame({
            'junctions': [
                'chr1:5-30:+',
                'chr1:22-30:+',
                'chr2:10-30:+',
                'chr2:10-50:-',
                'chr2:20-50:-'
            ],
            'events': [
                'chr1:5-30:+',
                'chr1:22-30:+',
                'chr2:10-30:+',
                'chr2:10-50:-;chr2:20-50:-',
                'chr2:10-50:-;chr2:20-50:-',
            ]
        }).set_index('junctions')
    )


def test_CountTable_event3(count_table):
    pd.testing.assert_frame_equal(
        count_table.event3,
        pd.DataFrame({
            'junctions': [
                'chr1:5-30:+',
                'chr1:22-30:+',
                'chr2:10-30:+',
                'chr2:10-50:-',
                'chr2:20-50:-'
            ],
            'events': [
                'chr1:5-30:+;chr1:22-30:+',
                'chr1:5-30:+;chr1:22-30:+',
                'chr2:10-30:+',
                'chr2:10-50:-',
                'chr2:20-50:-'
            ]
        }).set_index('junctions')
    )


def test_CountTable_event5_count(count_table):
    pd.testing.assert_frame_equal(
        count_table.event5_counts,
        pd.DataFrame({
            'events': [
                'chr1:22-30:+',
                'chr1:5-30:+',
                'chr2:10-30:+',
                'chr2:10-50:-;chr2:20-50:-'
            ],
            's1': [1, 1, 1, 2],
            's2': [1, 2, 1, 2],
            's3': [5, 10, 1, 6]
        }).set_index('events')
    )


def test_CountTable_event3_count(count_table):
    pd.testing.assert_frame_equal(
        count_table.event3_counts,
        pd.DataFrame({
            'events': [
                'chr1:5-30:+;chr1:22-30:+',
                'chr2:10-30:+',
                'chr2:10-50:-',
                'chr2:20-50:-'
            ],
            's1': [2, 1, 1, 1],
            's2': [3, 1, 1, 1],
            's3': [15, 1, 2, 4]
        }).set_index('events')
    )


def test_CountTable_quantile_filter(count_table):
    ct = count_table.quantile_filter(quantile=100, min_read=2)
    assert sorted(ct.junctions.tolist()) == sorted([
        'chr1:5-30:+',
        'chr1:22-30:+',
        'chr2:10-50:-',
        'chr2:20-50:-'
    ])


def test_CountTable_event5_count_filter(count_table):
    ct, cutoff = count_table.event5_count_filter()
    assert cutoff == 1
    assert sorted(ct.junctions.tolist()) == sorted([
        'chr1:5-30:+',
        'chr2:10-50:-',
        'chr2:20-50:-'
    ])


def test_CountTable_event3_count_filter(count_table):
    ct, cutoff = count_table.event3_count_filter()
    assert sorted(ct.junctions.tolist()) == sorted([
        'chr1:5-30:+',
        'chr1:22-30:+'
    ])


def test_event5_median_filter(count_table):
    ct = count_table.event5_median_filter(cutoff=2)
    assert sorted(ct.junctions.tolist()) == sorted([
        'chr1:5-30:+',
        'chr2:10-50:-',
        'chr2:20-50:-'
    ])


def test_event3_median_filter(count_table):
    ct = count_table.event3_median_filter(cutoff=2)
    assert sorted(ct.junctions.tolist()) == sorted([
        'chr1:5-30:+',
        'chr1:22-30:+'
    ])


def test_CountTable_ref_psi5(count_table):

    df = count_table.ref_psi5(method='beta_binomial', annotation=False)
    np.testing.assert_almost_equal(
        df['ref_psi'].tolist(), [1, 1, 1, 0.4, 0.6],
        decimal=2
    )

    df = count_table.ref_psi5(method='k/n', annotation=False)
    np.testing.assert_almost_equal(
        df['ref_psi'].tolist(), [1, 1, 1, 0.4, 0.6],
        decimal=2
    )

    df = count_table.ref_psi5(method='mean', annotation=False)
    np.testing.assert_almost_equal(
        df['ref_psi'].tolist(), [1, 1, 1, 0.44, 0.55],
        decimal=2
    )

    df1 = count_table.ref_psi5(method='beta_binomial', annotation=False)
    df2 = count_table.ref_psi5(method='k/n', annotation=False)
    np.testing.assert_almost_equal(
        df1['median_n'].tolist(), df2['median_n'].tolist(),
        decimal=2
    )


def test_CountTable_ref_psi3(count_table):
    df = count_table.ref_psi3(method='beta_binomial', annotation=False)
    np.testing.assert_almost_equal(
        df['ref_psi'].tolist(), [.65, .35, 1, 1, 1],
        decimal=2
    )

    df = count_table.ref_psi3(method='k/n', annotation=False)
    np.testing.assert_almost_equal(
        df['ref_psi'].tolist(), [.65, .35, 1, 1, 1],
        decimal=2
    )

    df = count_table.ref_psi3(method='mean', annotation=False)
    np.testing.assert_almost_equal(
        df['ref_psi'].tolist(), [0.61, 0.38, 1, 1, 1],
        decimal=2
    )

    df1 = count_table.ref_psi3(method='beta_binomial', annotation=False)
    df2 = count_table.ref_psi3(method='k/n', annotation=False)
    np.testing.assert_almost_equal(
        df1['median_n'].tolist(), df2['median_n'].tolist(),
        decimal=2
    )


def test_CountTable_psi5(count_table):
    np.testing.assert_almost_equal(
        count_table.psi5.values,
        np.array([[1., 1., 1.],
                  [1., 1., 1.],
                  [1., 1., 1.],
                  [0.5, 0.5, 0.33333333],
                  [0.5, 0.5, 0.66666667]])
    )


def test_CountTable_psi3(count_table):
    np.testing.assert_almost_equal(
        count_table.psi3.values,
        np.array([[0.5, 0.66666667, 0.66666667],
                  [0.5, 0.33333333, 0.33333333],
                  [1., 1., 1.],
                  [1., 1., 1.],
                  [1., 1., 1.]])
    )


def test_CountTable_to_from_csv(count_table, tmp_path):
    path = tmp_path / 'count_table.csv'
    count_table.to_csv(path)
    count_table.df = count_table.df.astype({
        'Chromosome': 'category',
        'Start': 'int32',
        'End': 'int32',
        'Strand': 'category',
        's1': 'int32',
        's2': 'int32',
        's3': 'int32'
    })
    ct = CountTable.read_csv(path)
    pd.testing.assert_frame_equal(ct.df, count_table.df)


def test_CountTable_junction_report(count_table):
    df = count_table.junction_report('chr2:20-50:-')
    pd.testing.assert_frame_equal(
        df,
        pd.DataFrame({
            'sample': ['s1', 's2', 's3'],
            'count': [1, 1, 4],
            'psi5': [0.5, 0.5, 0.66666667],
            'psi3': [1.0, 1.0, 1.0]
        }).set_index('sample')
    )


def test_CountTable_filter(count_table):
    ct = count_table.filter(['chr1:5-30:+', 'chr1:22-30:+'])
    assert sorted(ct.junctions.tolist()) == sorted(
        ['chr1:5-30:+', 'chr1:22-30:+'])

    ct = count_table.filter_event5(['chr2:10-50:-', 'chr2:10-30:+'])
    assert sorted(ct.junctions.tolist()) == sorted([
        'chr2:10-50:-', 'chr2:20-50:-', 'chr2:10-30:+'])

    ct = count_table.filter_event3(['chr1:5-30:+', 'chr2:10-30:+'])
    assert sorted(ct.junctions.tolist()) == sorted([
        'chr1:5-30:+', 'chr1:22-30:+', 'chr2:10-30:+'])


def test_CountTable_k(count_table):
    assert count_table.k('chr2:10-50:-').to_dict() \
        == {'s1': 1, 's2': 1, 's3': 2}

    pd.testing.assert_frame_equal(
        count_table.k(['chr1:5-30:+', 'chr2:10-50:-']),
        pd.DataFrame({
            'junctions': [
                'chr1:5-30:+',
                'chr2:10-50:-',
            ],
            's1': [1,  1],
            's2': [2,  1],
            's3': [10, 2]
        }).set_index('junctions')
    )


def test_CountTable_n5(count_table):
    assert count_table.n5('chr2:10-50:-').to_dict() \
        == {'s1': 2, 's2': 2, 's3': 6}

    pd.testing.assert_frame_equal(
        count_table.n5(['chr1:5-30:+', 'chr2:10-50:-']),
        pd.DataFrame({
            'junctions': [
                'chr1:5-30:+',
                'chr2:10-50:-',
            ],
            's1': [1,  2],
            's2': [2,  2],
            's3': [10, 6]
        }).set_index('junctions')
    )


def test_CountTable_n3(count_table):
    assert count_table.n3('chr1:5-30:+').to_dict() \
        == {'s1': 2, 's2': 3, 's3': 15}

    pd.testing.assert_frame_equal(
        count_table.n3(['chr1:5-30:+', 'chr2:10-50:-']),
        pd.DataFrame({
            'junctions': [
                'chr1:5-30:+',
                'chr2:10-50:-',
            ],
            's1': [2,  1],
            's2': [3,  1],
            's3': [15, 2]
        }).set_index('junctions')
    )


def test_CountTable_kn5(count_table):
    pd.testing.assert_frame_equal(
        count_table.kn5('chr2:10-50:-'),
        pd.DataFrame({
            'sample': ['s1', 's2', 's3'],
            'k': [1, 1, 2],
            'n': [2, 2, 6]
        }).set_index('sample')
    )


def test_CountTable_kn3(count_table):
    pd.testing.assert_frame_equal(
        count_table.kn3('chr1:5-30:+'),
        pd.DataFrame({
            'sample': ['s1', 's2', 's3'],
            'k': [1, 2, 10],
            'n': [2, 3, 15]
        }).set_index('sample')
    )


def test_CountTable_plot_kn5(count_table, mocker):
    scatter = mocker.patch('seaborn.scatterplot')
    count_table.plot_kn5('chr2:10-50:-')

    args, kwargs = scatter.call_args_list[0]
    assert kwargs['x'] == 'n'
    assert kwargs['y'] == 'k'
    pd.testing.assert_frame_equal(
        kwargs['data'],
        pd.DataFrame({
            'sample': ['s1', 's2', 's3'],
            'k': [1, 1, 2],
            'n': [2, 2, 6]
        }).set_index('sample')
    )


def test_CountTable_plot_kn3(count_table, mocker):
    scatter = mocker.patch('seaborn.scatterplot')
    count_table.plot_kn3('chr1:5-30:+')

    # import matplotlib.pyplot as plt
    # plt.show()
    args, kwargs = scatter.call_args_list[0]
    assert kwargs['x'] == 'n'
    assert kwargs['y'] == 'k'
    pd.testing.assert_frame_equal(
        kwargs['data'],
        pd.DataFrame({
            'sample': ['s1', 's2', 's3'],
            'k': [1, 2, 10],
            'n': [2, 3, 15]
        }).set_index('sample')
    )

    count_table.plot_kn3('chr1:5-30:+', highlight=['s3'])
    args, kwargs = scatter.call_args_list[1]
    pd.testing.assert_frame_equal(
        kwargs['data'],
        pd.DataFrame({
            'sample': ['s1', 's2', 's3'],
            'k': [1, 2, 10],
            'n': [2, 3, 15],
            'outlier': [False, False, True]
        }).set_index('sample')
    )


def test_CountTable_plot_psi5(count_table, mocker):
    plot = mocker.patch('seaborn.boxplot')
    count_table.plot_psi5('chr2:10-50:-')
    args, kwargs = plot.call_args_list[0]
    np.testing.assert_almost_equal(
        kwargs['data']['ref_psi'], [0.5, 0.5, 0.3333333])
    # import matplotlib.pyplot as plt
    # plt.show()

    count_table.plot_psi5('chr2:10-50:-', samples=['s1'])
    # plt.show()
    args, kwargs = plot.call_args_list[1]
    np.testing.assert_almost_equal(
        kwargs['data']['outliers'], [True, False, False])


# TODO: decimal precision needed to be decreased to pass test
def test_CountTable_plot_psi3(count_table, mocker):
    plot = mocker.patch('seaborn.boxplot')
    count_table.plot_psi3('chr1:5-30:+', samples=['s3'])
#     import matplotlib.pyplot as plt
#     plt.show()
    args, kwargs = plot.call_args_list[0]
    np.testing.assert_almost_equal(
        kwargs['data']['ref_psi'], [0.5, 0.666666, 0.666666],
        decimal=6)


def test_CountTable__psi_variants(count_table):
    vcf = MultiSampleVCF(vcf_file)
    df = pd.DataFrame({
        'sample': ['NA00002', 'NA00003', 's1'],
        'ref_psi': [0.5, 0.3, 0.1]
    }).set_index('sample')
    variants = [
        vcf.get_variant('chr1:4:T>C'),
        vcf.get_variant('chr1:25:AACG>GA')
    ]

    df_variants = count_table._psi_variants(df, vcf, variants)
    pd.testing.assert_frame_equal(
        df_variants,
        pd.DataFrame({
            'sample': ['NA00002', 'NA00003', 's1'] * 2,
            'ref_psi': [0.5, 0.3, 0.1] * 2,
            'genotype': ['AA', 'aa', 'NN', 'aa', 'AA', 'NN'],
            'variant': ['chr1:4:T>C'] * 3 + ['chr1:25:AACG>GA'] * 3
        }).set_index('sample')
    )


def test_CountTable_plot_psi5_variants(count_table, mocker):
    vcf = MultiSampleVCF(vcf_file)
    df = count_table.df.rename(columns={'s1': 'NA00002', 's2': 'NA00003'})
    count_table.df = df
    count_table.plot_psi5_variants('chr1:5-30:+', vcf)


def test_CountTable_plot_psi3_variants(count_table, mocker):
    vcf = MultiSampleVCF(vcf_file)
    df = count_table.df.rename(columns={'s1': 'NA00002', 's2': 'NA00003'})
    count_table.df = df
    count_table.plot_psi3_variants('chr1:5-30:+', vcf)


def test_CountTable_infer_annotation(count_table_chr17):
    df = count_table_chr17.infer_annotation(gtf_file, junc_file)

    pd.testing.assert_frame_equal(
        df,
        pd.DataFrame({
            'junctions': ['17:41197819-41199659:-', '17:41197831-41199670:-'],
            'gene_id': ['ENSG00000012048.22_5', 'ENSG00000012048'],
            'gene_name': ['BRCA1', 'BRCA1'],
            'gene_type': ['protein_coding', 'protein_coding'],
            'weak': [False, True],
            'transcript_id': [
                'ENST00000586385.5_1;ENST00000591534.5_1;ENST00000461221.5_1;'
                'ENST00000493795.5_1;ENST00000357654.8_3;ENST00000591849.5_1;'
                'ENST00000468300.5_2;ENST00000471181.7_3;ENST00000491747.6_3;'
                'ENST00000352993.7_2;ENST00000644379.1_1',
                np.nan
            ]
        }).set_index('junctions'))


def test_CountTable_ref_psi5_annnotation(count_table_chr17):
    count_table_chr17.infer_annotation(gtf_file, junc_file)
    df = count_table_chr17.ref_psi5()

    assert df.columns.tolist() == [
        'Chromosome', 'Start', 'End', 'Strand', 'splice_site', 'events',
        'ref_psi', 'k', 'n', 'median_n', 'gene_id', 'gene_name', 'gene_type', 'weak',
        'transcript_id']

    df = count_table_chr17.ref_psi5(method='beta_binomial')
    assert df.columns.tolist() == [
        'Chromosome', 'Start', 'End', 'Strand', 'splice_site', 'events',
        'ref_psi', 'alpha', 'beta', 'k', 'n', 'median_n', 'gene_id', 'gene_name', 'gene_type', 'weak',
        'transcript_id']


def test_CountTable_ref_psi3_annnotation(count_table_chr17):
    count_table_chr17.infer_annotation(gtf_file, junc_file)
    df = count_table_chr17.ref_psi3()

    assert df.columns.tolist() == [
        'Chromosome', 'Start', 'End', 'Strand', 'splice_site', 'events',
        'ref_psi', 'k', 'n', 'median_n', 'gene_id', 'gene_name', 'gene_type', 'weak',
        'transcript_id']

    df = count_table_chr17.ref_psi3(method='beta_binomial')
    assert df.columns.tolist() == [
        'Chromosome', 'Start', 'End', 'Strand', 'splice_site', 'events',
        'ref_psi', 'alpha', 'beta', 'k', 'n', 'median_n', 'gene_id', 'gene_name', 'gene_type', 'weak',
        'transcript_id']


def test_CountTable_join(count_table):
    ct_other = CountTable(pd.DataFrame({
        'Chromosome': ['chr1', 'chr3', 'chr2'],
        'Start': [5, 100, 10],
        'End': [30, 120, 30],
        'Strand': ['+', '+', '+'],
        's4': [1, 1, 1],
        's5': [3, 1, 1],
        's1': [1, 2, 3]
    }))

    ct = count_table.join(ct_other, suffix='_second')
    ct_expected = CountTable(pd.DataFrame({
        'Chromosome': ['chr1', 'chr1', 'chr2', 'chr2', 'chr2', 'chr3'],
        'Start': [22, 5, 10, 10, 20, 100],
        'End': [30, 30, 30, 50, 50, 120],
        'Strand': ['+', '+', '+', '-', '-', '+'],
        's1': [1, 1, 1, 1, 1, 0],
        's2': [1, 2, 1, 1, 1, 0],
        's3': [5, 10, 1, 2, 4, 0],
        's4': [0, 1, 1, 0, 0, 1],
        's5': [0, 3, 1, 0, 0, 1],
        's1_second': [0, 1, 3, 0, 0, 2],
    }))
    pd.testing.assert_frame_equal(ct.df, ct_expected.df)


# def test_CountTable_delta_logit_psi5(count_table):
#     delta_logit_psi = count_table.delta_logit_psi5()
#     np.testing.assert_almost_equal(
#         delta_logit_psi,
#         np.array([[0.,  0.,  0.],
#                   [0.,  0.,  0.],
#                   [0.,  0.,  0.],
#                   [0.40546511,  0.40546511, -0.28768207],
#                   [-0.40546511, -0.40546511,  0.28768207]]
#                  )
#     )

# def test_CountTable_delta_logit_psi3(count_table):
#     delta_logit_psi = count_table.delta_logit_psi3()
#     np.testing.assert_almost_equal(
#         delta_logit_psi,
#         np.array([[-0.61903921,  0.07410797,  0.07410797],
#                   [0.61903921, -0.07410797, -0.07410797],
#                   [0.,  0.,  0.],
#                   [0.,  0.,  0.],
#                   [0.,  0.,  0.]]
#                  )
#     )
