#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 15:31:49 2023

@author: albertsmith
"""



import sys
import os

if 'google.colab' in sys.modules:
    try:
        import MDAnalysis
    except:
        os.popen('pip3 install MDAnalysis')
    
    try:
        import pyDR
    except:
        os.popen('git clone https://github.com/alsinmr/pyDR.git')
            