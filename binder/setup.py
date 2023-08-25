# -*- coding: utf-8 -*-

from setuptools import setup
  
setup(
    name='pyDR Tutorial',
    version='0.1',
    description='pyDR Jupyter Tutorial',
    author='Albert Smith-Penzel',
    author_email='albert.smith-penzel@medizin.uni-leipzig.de',
    packages=['Tpl_basicNMR'],
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'MDAnalysis',
        'nglview'
    ],
)

import os
os.popen('git clone https://github.com/alsinmr/pyDR.git')