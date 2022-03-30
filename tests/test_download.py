from click.testing import CliRunner
from splicemap.splice_map import SpliceMap
from splicemap.main import splicemap_download


def tests_download_splicemap(tmp_path):
    runner = CliRunner()
    
    result = runner.invoke(splicemap_download, f'--version _test --splicemap_dir {str(tmp_path)}')
    assert result.exit_code == 0

    sm = SpliceMap.read_csv(str(tmp_path))
    assert sm.df.shape[0] > 0
    
