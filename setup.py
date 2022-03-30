#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'setuptools',
    'numpy',
    'pandas',
    'tqdm',
    'seaborn',
    'sklearn',
    'matplotlib',
    'kipoiseq>=0.3.0',
    'fastbetabino3',
    'pyranges>=0.0.71',
    'wget',
]


setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', 'pytest-mock', 'cyvcf2']

setup(
    author="M. Hasan Celik & Nils Wagner",
    author_email="wagnern@in.tum.de",
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'
    ],
    description="SpliceMap generation from split reads from RNA-seq",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    long_description_content_type='text/markdown',
    include_package_data=True,
    keywords='splicemap',
    name='splicemap',
    packages=find_packages(include=['splicemap']),
    setup_requires=setup_requirements,
    entry_points='''
        [console_scripts]
        splicemap_download=splicemap.main:splicemap_download
    ''',
    #     dependency_links=['http://github.com/lfiaschi/fastbetabino/tarball/master#egg=fastbetabino'],
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/gagneurlab/splicemap',
    version='0.0.1',
    zip_safe=False
)
