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
]


setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', 'pytest-benchmark']

setup(
    author="Hasan Celik",
    author_email="wagnern@in.tum.de",
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    description="Count table for sequencing data",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    long_description_content_type='text/markdown',
    include_package_data=True,
    keywords='count_table',
    name='count_table',
    packages=find_packages(include=['count_table']),
    setup_requires=setup_requirements,
    dependency_links=['http://github.com/lfiaschi/fastbetabino/tarball/master#egg=fastbetabino'],
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/gagneurlab/count_table',
    version='1.0.0',
    zip_safe=False
)
