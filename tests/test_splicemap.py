import pytest
from splicemap.splice_map import SpliceMap
from conftest import ref_table5_kn_testis, ref_table3_kn_testis


@pytest.fixture
def splicemap5_kn():
    return SpliceMap.read_csv(ref_table5_kn_testis)


@pytest.fixture
def splicemap3_kn():
    return SpliceMap.read_csv(ref_table3_kn_testis)


def test_Splicemap__init__(splicemap5_kn, splicemap3_kn):
    assert splicemap5_kn.name == 'gtex-grch37-testis-psi5'
    assert splicemap3_kn.name == 'gtex-grch37-testis-psi3'

    assert splicemap5_kn.method == 'kn'
    assert splicemap5_kn.df.shape[0] == 24

    assert splicemap3_kn.method == 'kn'
    assert splicemap3_kn.df.shape[0] == 45


def test_Splicemap_valid(splicemap5_kn, splicemap3_kn):
    pass
