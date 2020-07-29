# CountTable
[![CircleCI]](https://gitlab.cmm.in.tum.de/gagneurlab/count_table/pipelines)
<!-- [![pypi](https://img.shields.io/pypi/v/mmsplice.svg)](https://pypi.python.org/pypi/mmsplice) -->

A package to process split-read counts for splicing. 
CountTable requires a specific format of RNA-seq split-read count data, namely: 'Chromosome', 'Start', 'End', 'Strand', 'samples':
.
|   Chromosome |   Start |    End | Strand   |   GTEX-111FC |   GTEX-1128S |   GTEX-117XS |   GTEX-1192X |   GTEX-11DXW |   GTEX-11DXY |
|-------------:|--------:|-------:|:---------|-------------:|-------------:|-------------:|-------------:|-------------:|-------------:|
|            1 |   17729 |  17733 | +        |            9 |            0 |            0 |            0 |            0 |            7 |
|            1 |   30667 |  30975 | +        |            8 |            1 |            7 |            3 |            5 |            5 |
|            1 |  135802 | 137620 | +        |            0 |            2 |            2 |            0 |            3 |            1 |
|            1 |  320653 | 320880 | +        |            1 |            1 |            4 |            1 |            2 |            4 |
|            1 |  320653 | 324287 | +        |            0 |            2 |            8 |            1 |            0 |            2 |
|            1 |  320938 | 321031 | +        |            2 |            4 |            5 |            3 |            5 |            1 |
|            1 |  320938 | 322037 | +        |            8 |            5 |            4 |            6 |            8 |            6 |
|            1 |  322228 | 324287 | +        |           53 |           27 |           40 |           17 |           35 |           33 |
|            1 |  324345 | 324438 | +        |           99 |           54 |          101 |           38 |           82 |           63 |
|            1 |  324686 | 324718 | +        |            0 |            3 |            8 |            2 |            1 |            2 |


## Installation
-----------------

Create virtual environment:
```bash
pip install virtualenv
virtualenv venv
source venv/bin/activate
```

External dependencies:
```bash
pip install cyvcf2 cython pytest pytest-runner
```

Install all packages from setup.py:
```bash
pip install -e .
```
-----------------

### Example code
-------------------

Check [notebooks/example.ipynb](https://gitlab.cmm.in.tum.de/gagneurlab/count_table/-/tree/master/notebooks/example.ipynb)