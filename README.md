# SpliceMap

A package to process RNA-seq split-read counts for splicing and generate tissue-specific splicing annotations (SpliceMaps) that quantify tissue-specific splice site usage and isoform competition.
SpliceMaps can be generated from RNA-seq split read counts (e.g. as provided by the software FRASER). The required input format is a table of split read-counts, with columns: 'Chromosome', 'Start', 'End', 'Strand', 'samples':

|   Chromosome |   Start |    End | Strand   | individual 1 | individual 2 | individual 3 | 
|-------------:|--------:|-------:|:---------|-------------:|-------------:|-------------:|
|            1 |   17729 |  17733 | +        |            9 |            0 |            0 |
|            1 |   30667 |  30975 | +        |            8 |            1 |            7 |
|            1 |  135802 | 137620 | +        |            0 |            2 |            2 |
|            1 |  320653 | 320880 | +        |            1 |            1 |            4 |
|            1 |  320653 | 324287 | +        |            0 |            2 |            8 |
|            1 |  320938 | 321031 | +        |            2 |            4 |            5 |
|            1 |  320938 | 322037 | +        |            8 |            5 |            4 |
|            1 |  322228 | 324287 | +        |           53 |           27 |           40 |
|            1 |  324345 | 324438 | +        |           99 |           54 |          101 |
|            1 |  324686 | 324718 | +        |            0 |            3 |            8 |

SpliceMaps contain the following information:

| column name | Description |
| --------  | ----------- |
| junctions | Intron with: chrom:donor site:acceptor site:strand |
| gene_id | Ensembl gene id |
| gene_name | Gene symbol |
| gene_type | Gene type |
| gene_tpm | Median expression of gene |
| transcript_id | List of transcripts that intron is part of |
| Chromosome | Chromosome |
| Start | Donor site of intron |
| End | Acceptor site of intron |
| Strand | Strand of intron |
| events | List of competing introns with shared splice site |
| splice_site | Fixed splice site in splicing event |
| ref_psi | Reference level of PSI |
| k | Sum of split-read counts supporting the intron (based on cohort that SpliceMap was computed on) |
| n | Sum of split-read counts supporting the splice-site (based on cohort that SpliceMap was computed on) |
| median_n | Median of n across individuals |
| novel_junction | Intron is not annotated in Gencode (binary) |
| weak_site_donor | Donor site is not annotated in Gencode (binary) |
| weak_site_acceptor | Acceptor site is not annotated in Gencode (binary) |


## Installation
-----------------
Clone git repository of splicemap:
```
git clone git@github.com:gagneurlab/splicemap.git

```
cd into repo directory:
```
cd splicemap
```

Install conda environment:
```
# Recommended if you have mamba installed
mamba env create -f environment.yaml
# otherwise
conda env create -f environment.yaml
```
Activate conda environment:
```
conda activate splicemap
```


### Example usage
-------------------

Check [notebooks/example.ipynb](https://github.com/gagneurlab/splicemap/blob/master/notebooks/example.ipynb)


### Download precomputed SpliceMaps from Zenodo
```bash
splicemap_download --version {version} --splicemap_dir {output_dir}
```
Supported versions: 'gtex_v8' (hg38) and 'gtex_v7' (hg19).