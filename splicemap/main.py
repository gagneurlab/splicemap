import wget
import click
import tarfile
import requests
from tqdm import tqdm
from pathlib import Path
from urllib.parse import unquote

gtex_v8_tissues = [
        'Adipose_Subcutaneous',
        'Adipose_Visceral_Omentum',
        'Adrenal_Gland',
        'Artery_Aorta',
        'Artery_Coronary',
        'Artery_Tibial',
        'Brain_Amygdala',
        'Brain_Anterior_cingulate_cortex_BA24',
        'Brain_Caudate_basal_ganglia',
        'Brain_Cerebellar_Hemisphere',
        'Brain_Cerebellum',
        'Brain_Cortex',
        'Brain_Frontal_Cortex_BA9',
        'Brain_Hippocampus',
        'Brain_Hypothalamus',
        'Brain_Nucleus_accumbens_basal_ganglia',
        'Brain_Putamen_basal_ganglia',
        'Brain_Spinal_cord_cervical_c_1',
        'Brain_Substantia_nigra',
        'Breast_Mammary_Tissue',
        'Cells_Cultured_fibroblasts',
        'Cells_EBV_transformed_lymphocytes',
        'Colon_Sigmoid',
        'Colon_Transverse',
        'Esophagus_Gastroesophageal_Junction',
        'Esophagus_Mucosa',
        'Esophagus_Muscularis',
        'Heart_Atrial_Appendage',
        'Heart_Left_Ventricle',
        'Kidney_Cortex',
        'Liver',
        'Lung',
        'Minor_Salivary_Gland',
        'Muscle_Skeletal',
        'Nerve_Tibial',
        'Ovary',
        'Pancreas',
        'Pituitary',
        'Prostate',
        'Skin_Not_Sun_Exposed_Suprapubic',
        'Skin_Sun_Exposed_Lower_leg',
        'Small_Intestine_Terminal_Ileum',
        'Spleen',
        'Stomach',
        'Testis',
        'Thyroid',
        'Uterus',
        'Vagina',
        'Whole_Blood',
]

zenodo_base_url = 'https://zenodo.org/record/6387938/files/'

def complete_url(tissue, psi):
    return zenodo_base_url + tissue + '_splicemap_' + psi + '_method%3Dkn_event_filter%3Dmedian_cutoff.csv.gz?download=1'

splicemap_url = {
    'gtex_v8': [
        *[complete_url(tissue, 'psi5') for tissue in gtex_v8_tissues],
        *[complete_url(tissue, 'psi3') for tissue in gtex_v8_tissues]
    ]
}

def check_tissue_in_url(url, tissues):
    in_url = False
    for tissue in tissues:
        if tissue in url:
            in_url = True
    return in_url

def _download(url, splicemap_dir):
    r = requests.get(url)
    fname = Path(splicemap_dir)
    fname.mkdir(exist_ok=True)
    fname = fname / \
        unquote(Path(url).name).replace('?download=1', '')
    print(fname)
    with open(fname, 'wb') as fd:
        fd.write(r.content)

@click.command()
@click.option('--version', help='SpliceMap version (currently SpliceMaps from gtex_v8 (hg38) supported)')
@click.option('--splicemap_dir', help='Path to download SpliceMaps')
@click.option('--tissues', multiple=True, help='List of tissue names to download')
def splicemap_download(version, splicemap_dir, tissues=None):

    if version not in splicemap_url:
        raise(f'Version {version} is not supported.')

    print('Downloading SpliceMaps...')
    
    for url in tqdm(splicemap_url[version]):
        if tissues:
            if check_tissue_in_url(url, tissues):
                _download(url, splicemap_dir)
        else:
            _download(url, splicemap_dir)
