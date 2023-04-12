from click.testing import CliRunner
from splicemap.splice_map import SpliceMap
from splicemap.main import splicemap_download, gtex_v7_tissues, gtex_v8_tissues
import os

def tests_download_splicemap_hg38(tmp_path):
    runner = CliRunner()
    
    result = runner.invoke(splicemap_download, f'--version gtex_v8 --splicemap_dir {str(tmp_path)}')
    assert result.exit_code == 0
    
    for tissue in gtex_v8_tissues:
        for psi in ['psi5', 'psi3']:
            sm = SpliceMap.read_csv(os.path.join(
                    tmp_path, 
                    f'{tissue}_splicemap_{psi}_method=kn_event_filter=median_cutoff.csv.gz'))
            assert sm.df.shape[0] > 0
    
    
def tests_download_splicemap_single_tissue_hg38(tmp_path):
    runner = CliRunner()

    result = runner.invoke(
        splicemap_download, f'--version gtex_v8 --splicemap_dir {str(tmp_path)} --tissues Whole_Blood --tissues Lung')
    assert result.exit_code == 0

    for tissue in ['Whole_Blood', 'Lung']:
        for psi in ['psi5', 'psi3']:
            sm = SpliceMap.read_csv(os.path.join(
                    tmp_path, 
                    f'{tissue}_splicemap_{psi}_method=kn_event_filter=median_cutoff.csv.gz'))
            assert sm.df.shape[0] > 0
            
            
def tests_download_splicemap_hg19(tmp_path):
    runner = CliRunner()
    
    result = runner.invoke(splicemap_download, f'--version gtex_v7 --splicemap_dir {str(tmp_path)}')
    assert result.exit_code == 0
    
    for tissue in gtex_v7_tissues:
        for psi in ['psi5', 'psi3']:
            sm = SpliceMap.read_csv(os.path.join(
                    tmp_path, 
                    f'{tissue}_splicemap_{psi}_method=kn_event_filter=median_cutoff.csv.gz'))
            assert sm.df.shape[0] > 0
  
    
def tests_download_splicemap_single_tissue_hg19(tmp_path):
    runner = CliRunner()

    result = runner.invoke(
        splicemap_download, f'--version gtex_v7 --splicemap_dir {str(tmp_path)} --tissues Whole_Blood --tissues Lung')
    assert result.exit_code == 0

    for tissue in ['Whole_Blood', 'Lung']:
        for psi in ['psi5', 'psi3']:
            sm = SpliceMap.read_csv(os.path.join(
                    tmp_path, 
                    f'{tissue}_splicemap_{psi}_method=kn_event_filter=median_cutoff.csv.gz'))
            assert sm.df.shape[0] > 0
