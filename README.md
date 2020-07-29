# CountTable
[![CircleCI]](https://gitlab.cmm.in.tum.de/gagneurlab/count_table/pipelines)
<!-- [![pypi](https://img.shields.io/pypi/v/mmsplice.svg)](https://pypi.python.org/pypi/mmsplice) -->

A package to process split read counts for splicing. 
CountTable requires a specific format of RNA-seq count data, namely: 'Chromosome', 'Start', 'End', 'Strand', samples.


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
