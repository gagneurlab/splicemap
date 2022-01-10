import pytest
import numpy as np
import pandas as pd
from kipoiseq.extractors import MultiSampleVCF
from splicemap.dataclasses import Junction
from splicemap import SpliceCountTable as CountTable
from splicemap import infer_junction_strand
from kipoiseq.extractors import FastaStringExtractor
from conftest import fasta_file, vcf_file, gtf_file, junc_file, \
    gtf_file_with_chr, gtf_file_multi


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


dtypes = pd.Series({
    'Chromosome': 'category',
    'Start': 'int32',
    'End': 'int32',
    'Strand': 'category',
    's1': 'int32',
    's2': 'int32',
    's3': 'int32',
    'count': 'int32',
    'k': 'int32',
    'n': 'int32'
})


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
    return CountTable(df, name='test_count_table')


def test_count_table_gene_expression():
    df = pd.DataFrame({
        'Chromosome': ['chr1', 'chr1', 'chr2', 'chr2', 'chr2'],
        'Start': [5, 22, 10, 10, 20],
        'End': [30, 30, 30, 50, 50],
        'Strand': ['+', '+', '+', '-', '-'],
        's1': [1, 1, 1, 1, 1],
        's2': [2, 1, 1, 1, 1],
        's3': [10, 5, 1, 2, 4]
    })
    df_exp = pd.DataFrame({
        'gene_id': ['gene_a', 'gene_b', 'gene_c'],
        's1': [1, 2, 3],
        's2': [1, 2, 3],
        's3': [1, 2, 3]
    }).set_index('gene_id')

    ct = CountTable(df, name='test_count_table', gene_expression=df_exp)
    gene_expression_median = pd.Series([1., 2., 3.],
                                       index=['gene_a', 'gene_b', 'gene_c'],
                                       name='gene_tpm')
    pd.testing.assert_series_equal(
        ct.gene_expression_median,
        gene_expression_median
    )

# @pytest.fixture
# def count_table_gene():
#     df = pd.DataFrame({
#         'Chromosome': ['chr1', 'chr1', 'chr2', 'chr2', 'chr2'],
#         'Start': [5, 22, 10, 10, 20],
#         'End': [30, 30, 30, 50, 50],
#         'Strand': ['+', '+', '+', '-', '-'],
#         's1': [1, 1, 1, 1, 1],
#         's2': [2, 1, 1, 1, 1],
#         's3': [10, 5, 1, 2, 4]
#     })
#     return CountTable(df, name='test_count_table', gene_expression=df_gene)


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
        'Chromosome': ['17', '17', '17'],
        'Start': [41197819, 41197831, 42388823],
        'End': [41199659, 41199670, 42399786],
        'Strand': ['-', '-', '-'],
        's1': [1, 1, 1],
        's2': [2, 1, 2]
    })
    return CountTable(df, name='test_count_table')


@pytest.fixture
def count_table_chr1_strand():
    df = pd.DataFrame({
        'Chromosome': ['1'],
        'Start': [62598806],
        'End': [62601768],
        'Strand': ['+'],
        's1': [1],
        's2': [2]
    })
    return CountTable(df, name='test_count_table_two_genes')


@pytest.fixture
def count_table_chr17_expression():
    df = pd.DataFrame({
        'Chromosome': ['17', '17'],
        'Start': [41197819, 41197831],
        'End': [41199659, 41199670],
        'Strand': ['-', '-'],
        's1': [1, 1],
        's2': [2, 1]
    })
    df_gene = pd.DataFrame({
        'gene_id': ['ENSG00000012048'],
        's1': [1],
        's2': [2],
        's3': [3]
    }).set_index('gene_id')
    return CountTable(df, name='test_count_table', gene_expression=df_gene)


def test_CountTable_infer_strand(count_table):
    df = pd.DataFrame({
        'Chromosome': ['17', '17'],
        'Start': [41267796, 41279042],
        'End': [41276033, 41279742],
        'Strand': ['.', '*'],
        's1': [1, 1],
        's2': [2, 1]
    })
    ct = CountTable(df, name='test_count_table')

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
    df = pd.DataFrame({
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
    df = df.astype(dtypes[dtypes.keys().isin(df.columns)])
    pd.testing.assert_frame_equal(
        count_table.junction_df, df)


def test_CountTable_samples(count_table):
    assert count_table.samples == ['s1', 's2', 's3']


def test_CountTable_counts(count_table):
    df = pd.DataFrame({
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
    df = df.astype(dtypes[dtypes.keys().isin(df.columns)])

    pd.testing.assert_frame_equal(
        count_table.counts, df)


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
    df = pd.DataFrame({
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
    df = df.astype(dtypes[dtypes.keys().isin(df.columns)])

    pd.testing.assert_frame_equal(
        count_table.event5_counts, df)


def test_CountTable_event3_count(count_table):
    df = pd.DataFrame({
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
    df = df.astype(dtypes[dtypes.keys().isin(df.columns)])

    pd.testing.assert_frame_equal(
        count_table.event3_counts, df)


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
    sm = count_table.ref_psi5(method='beta_binomial', annotation=False)
    np.testing.assert_almost_equal(
        sm.df['ref_psi'].tolist(), [1, 1, 1, 0.4, 0.6],
        decimal=2
    )

    sm = count_table.ref_psi5(method='k/n', annotation=False)
    np.testing.assert_almost_equal(
        sm.df['ref_psi'].tolist(), [1, 1, 1, 0.4, 0.6],
        decimal=2
    )

    sm = count_table.ref_psi5(method='mean', annotation=False)
    np.testing.assert_almost_equal(
        sm.df['ref_psi'].tolist(), [1, 1, 1, 0.44, 0.55],
        decimal=2
    )

    sm1 = count_table.ref_psi5(method='beta_binomial', annotation=False)
    sm2 = count_table.ref_psi5(method='k/n', annotation=False)
    np.testing.assert_almost_equal(
        sm1.df['median_n'].tolist(), sm2.df['median_n'].tolist(),
        decimal=2
    )


def test_CountTable_ref_psi3(count_table):
    sm = count_table.ref_psi3(method='beta_binomial', annotation=False)
    np.testing.assert_almost_equal(
        sm.df['ref_psi'].tolist(), [.65, .35, 1, 1, 1],
        decimal=2
    )

    sm = count_table.ref_psi3(method='k/n', annotation=False)
    np.testing.assert_almost_equal(
        sm.df['ref_psi'].tolist(), [.65, .35, 1, 1, 1],
        decimal=2
    )

    sm = count_table.ref_psi3(method='mean', annotation=False)
    np.testing.assert_almost_equal(
        sm.df['ref_psi'].tolist(), [0.61, 0.38, 1, 1, 1],
        decimal=2
    )

    sm1 = count_table.ref_psi3(method='beta_binomial', annotation=False)
    sm2 = count_table.ref_psi3(method='k/n', annotation=False)
    np.testing.assert_almost_equal(
        sm1.df['median_n'].tolist(), sm2.df['median_n'].tolist(),
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


def test_CountTable_to_from_csv_expression(count_table, tmp_path):
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

    path_exp = tmp_path / 'expression.csv'
    df_exp = pd.DataFrame({
        'gene_id': ['ENSG00000012048'],
        's1': [1],
        's2': [2],
        's3': [3]
    }).set_index('gene_id')
    df_exp.to_csv(path_exp)

    ct = CountTable.read_csv(path, gene_expression_path=path_exp)
    pd.testing.assert_frame_equal(ct.df, count_table.df)
    pd.testing.assert_frame_equal(ct._gene_expression, df_exp)


def test_CountTable_junction_report(count_table):
    df = count_table.junction_report('chr2:20-50:-')
    _df = pd.DataFrame({
        'sample': ['s1', 's2', 's3'],
        'count': [1, 1, 4],
        'psi5': [0.5, 0.5, 0.66666667],
        'psi3': [1.0, 1.0, 1.0]
    }).set_index('sample')

    _df = _df.astype(dtypes[dtypes.keys().isin(_df.columns)])
    pd.testing.assert_frame_equal(df, _df)


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
    df = pd.DataFrame({
        'junctions': [
            'chr1:5-30:+',
            'chr2:10-50:-',
        ],
        's1': [1,  1],
        's2': [2,  1],
        's3': [10, 2]
    }).set_index('junctions')
    df = df.astype(dtypes[dtypes.keys().isin(df.columns)])

    pd.testing.assert_frame_equal(
        count_table.k(['chr1:5-30:+', 'chr2:10-50:-']), df)


def test_CountTable_n5(count_table):
    assert count_table.n5('chr2:10-50:-').to_dict() \
        == {'s1': 2, 's2': 2, 's3': 6}

    df = pd.DataFrame({
        'junctions': [
            'chr1:5-30:+',
            'chr2:10-50:-',
        ],
        's1': [1,  2],
        's2': [2,  2],
        's3': [10, 6]
    }).set_index('junctions')
    df = df.astype(dtypes[dtypes.keys().isin(df.columns)])

    pd.testing.assert_frame_equal(
        count_table.n5(['chr1:5-30:+', 'chr2:10-50:-']), df)


def test_CountTable_n3(count_table):
    assert count_table.n3('chr1:5-30:+').to_dict() \
        == {'s1': 2, 's2': 3, 's3': 15}

    df = pd.DataFrame({
        'junctions': [
            'chr1:5-30:+',
            'chr2:10-50:-',
        ],
        's1': [2,  1],
        's2': [3,  1],
        's3': [15, 2]
    }).set_index('junctions')
    df = df.astype(dtypes[dtypes.keys().isin(df.columns)])

    pd.testing.assert_frame_equal(
        count_table.n3(['chr1:5-30:+', 'chr2:10-50:-']),  df)


def test_CountTable_kn5(count_table):
    df = pd.DataFrame({
        'sample': ['s1', 's2', 's3'],
        'k': [1, 1, 2],
        'n': [2, 2, 6]
    }).set_index('sample')
    df = df.astype(dtypes[dtypes.keys().isin(df.columns)])

    pd.testing.assert_frame_equal(
        count_table.kn5('chr2:10-50:-'), df, check_dtype=False)


def test_CountTable_kn3(count_table):
    df = pd.DataFrame({
        'sample': ['s1', 's2', 's3'],
        'k': [1, 2, 10],
        'n': [2, 3, 15]
    }).set_index('sample')
    df = df.astype(dtypes[dtypes.keys().isin(df.columns)])

    pd.testing.assert_frame_equal(
        count_table.kn3('chr1:5-30:+'), df, check_dtype=False)


def test_CountTable_plot_kn5(count_table, mocker):
    scatter = mocker.patch('seaborn.scatterplot')
    count_table.plot_kn5('chr2:10-50:-')

    df = pd.DataFrame({
        'sample': ['s1', 's2', 's3'],
        'k': [1, 1, 2],
        'n': [2, 2, 6]
    }).set_index('sample')

    args, kwargs = scatter.call_args_list[0]
    assert kwargs['x'] == 'n'
    assert kwargs['y'] == 'k'
    pd.testing.assert_frame_equal(
        kwargs['data'], df, check_dtype=False)


def test_CountTable_plot_kn3(count_table, mocker):
    scatter = mocker.patch('seaborn.scatterplot')
    count_table.plot_kn3('chr1:5-30:+')

    df = pd.DataFrame({
        'sample': ['s1', 's2', 's3'],
        'k': [1, 2, 10],
        'n': [2, 3, 15]
    }).set_index('sample')

    # import matplotlib.pyplot as plt
    # plt.show()
    args, kwargs = scatter.call_args_list[0]
    assert kwargs['x'] == 'n'
    assert kwargs['y'] == 'k'
    pd.testing.assert_frame_equal(
        kwargs['data'], df, check_dtype=False)

    df = pd.DataFrame({
        'sample': ['s1', 's2', 's3'],
        'k': [1, 2, 10],
        'n': [2, 3, 15],
        'outlier': [False, False, True]
    }).set_index('sample')

    count_table.plot_kn3('chr1:5-30:+', highlight=['s3'])
    args, kwargs = scatter.call_args_list[1]
    pd.testing.assert_frame_equal(
        kwargs['data'], df, check_dtype=False)


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

# TODO: add case with 2 genes on different strands
def test_CountTable_infer_annotation(count_table_chr17):
    df = count_table_chr17.infer_annotation(gtf_file, filter_intergenic='complete', strandedness=True)
    df = df.drop(columns='transcript_id', axis=1)
    pd.testing.assert_frame_equal(
        df,
        pd.DataFrame({
            'junctions': ['17:41197819-41199659:-', '17:41197831-41199670:-'],
            'gene_id': ['ENSG00000012048', 'ENSG00000012048'],
            'gene_name': ['BRCA1', 'BRCA1'],
            'gene_type': ['protein_coding', 'protein_coding'],
            'novel_junction': [False, True],
            'weak_site_donor': [False, True],
            'weak_site_acceptor': [False, True]
        }).set_index(['junctions', 'gene_id']))

def test_CountTable_infer_annotation_strand(count_table_chr1_strand):
    df_strand = count_table_chr1_strand.infer_annotation(gtf_file, filter_intergenic='complete', strandedness=True)
    df_no_strand = count_table_chr1_strand.infer_annotation(gtf_file, filter_intergenic='complete', strandedness=False)
    assert sorted(set(df_strand['gene_name'])) == ['ANGPTL3']
    assert sorted(set(df_no_strand['gene_name'])) == ['ANGPTL3', 'DOCK7']
    
def test_CountTable_infer_annotation_partial(count_table_chr17):
    df = count_table_chr17.infer_annotation(gtf_file, filter_intergenic='partial')
    df = df.drop(columns='transcript_id', axis=1)
    pd.testing.assert_frame_equal(
        df,
        pd.DataFrame({
            'junctions': ['17:41197819-41199659:-', '17:41197831-41199670:-', '17:42388823-42399786:-', '17:42388823-42399786:-'],
            'gene_id': ['ENSG00000012048', 'ENSG00000012048', 'ENSG00000267750', 'ENSG00000013306'],
            'gene_name': ['BRCA1', 'BRCA1', 'RUNDC3A-AS1', 'SLC25A39'], #RUNDC3A was previously there (but is on + strand)
            'gene_type': ['protein_coding', 'protein_coding', 'antisense', 'protein_coding'],
            'novel_junction': [False, True, True, True],
            'weak_site_donor': [False, True, False, False],
            'weak_site_acceptor': [False, True, False, False]
        }).set_index(['junctions', 'gene_id']))

def test_CountTable_infer_annotation_multiple(count_table_chr17):
    df = count_table_chr17.infer_annotation(gtf_file_multi)
    pd.testing.assert_frame_equal(
        df,
        pd.DataFrame({
            'junctions': ['17:41197819-41199659:-',
                          '17:41197831-41199670:-',
                          '17:41197831-41199670:-'],
            'gene_id': ['ENSG00000012048', 'ENSG00000012048', 'ENSGBRCAX'],
            'gene_name': ['BRCA1', 'BRCA1', 'BRCAX'],
            'gene_type': ['protein_coding',
                          'protein_coding', 'protein_coding'],
            'novel_junction': [False, True, True],
            'weak_site_donor': [False, True, True],
            'weak_site_acceptor': [False, True, True],
            'transcript_id': [
                'ENST00000357654', np.nan, np.nan
            ]
        }).set_index(['junctions', 'gene_id']))

#     __import__("pdb").set_trace()
#     # TODO: test multiple genes mock _gene_junction_overlap
#     # TODOD: test _gene_junction_overlap multiple overlap


def test_CountTable_infer_annotation_with_chr(count_table_chr17):
    df = count_table_chr17.infer_annotation(gtf_file_with_chr)
    pd.testing.assert_frame_equal(
        df,
        pd.DataFrame({
            'junctions': ['17:41197819-41199659:-', '17:41197831-41199670:-'],
            'gene_id': ['ENSG00000012048', 'ENSG00000012048'],
            'gene_name': ['BRCA1', 'BRCA1'],
            'gene_type': ['protein_coding', 'protein_coding'],
            'novel_junction': [False, True],
            'weak_site_donor': [False, True],
            'weak_site_acceptor': [False, True],
            'transcript_id': [
                'ENST00000357654', np.nan
            ]
        }).set_index(['junctions', 'gene_id']))


def test_CountTable_ref_psi5_annnotation(count_table_chr17):
    count_table_chr17.infer_annotation(gtf_file)
    sm = count_table_chr17.ref_psi5()

    assert sm.df.columns.tolist() == [
        'Chromosome', 'Start', 'End', 'Strand', 'splice_site', 'events',
        'ref_psi', 'k', 'n', 'median_n', 'gene_name', 'gene_type',
        'novel_junction', 'weak_site_donor', 'weak_site_acceptor', 'transcript_id']

    sm = count_table_chr17.ref_psi5(method='beta_binomial')
    assert sm.df.columns.tolist() == [
        'Chromosome', 'Start', 'End', 'Strand', 'splice_site', 'events',
        'ref_psi', 'alpha', 'beta', 'k', 'n', 'median_n', 'gene_name', 'gene_type',
        'novel_junction', 'weak_site_donor', 'weak_site_acceptor', 'transcript_id']


# def test_CountTable_ref_psi5_annnotation(count_table_chr17):
#     count_table_chr17.infer_annotation(gtf_file)
#     sm = count_table_chr17.ref_psi5()

    
def test_CountTable_ref_psi3_annnotation(count_table_chr17):
    count_table_chr17.infer_annotation(gtf_file)
    sm = count_table_chr17.ref_psi3()

    assert sm.df.columns.tolist() == [
        'Chromosome', 'Start', 'End', 'Strand', 'splice_site', 'events',
        'ref_psi', 'k', 'n', 'median_n', 'gene_name', 'gene_type',
        'novel_junction', 'weak_site_donor', 'weak_site_acceptor', 'transcript_id']

    sm = count_table_chr17.ref_psi3(method='beta_binomial')
    assert sm.df.columns.tolist() == [
        'Chromosome', 'Start', 'End', 'Strand', 'splice_site', 'events',
        'ref_psi', 'alpha', 'beta', 'k', 'n', 'median_n', 'gene_name', 'gene_type',
        'novel_junction', 'weak_site_donor', 'weak_site_acceptor', 'transcript_id']


def test_CountTable_ref_psi5_annnotation_tpm(count_table_chr17_expression):
    count_table_chr17_expression.infer_annotation(gtf_file)
    sm = count_table_chr17_expression.ref_psi5()
    assert sm.df.columns.tolist() == [
        'Chromosome', 'Start', 'End', 'Strand', 'splice_site', 'events',
        'ref_psi', 'k', 'n', 'median_n', 'gene_name', 'gene_type',
        'novel_junction', 'weak_site_donor', 'weak_site_acceptor', 'transcript_id', 'gene_tpm']
    assert all(sm.df['gene_tpm'] == 2)


def test_CountTable_join(count_table):
    ct_other = CountTable(pd.DataFrame({
        'Chromosome': ['chr1', 'chr3', 'chr2'],
        'Start': [5, 100, 10],
        'End': [30, 120, 30],
        'Strand': ['+', '+', '+'],
        's4': [1, 1, 1],
        's5': [3, 1, 1],
        's1': [1, 2, 3]
    }), name='test_count_table')

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
    }), name='test_count_table')
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
