#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 15:31:49 2023

@author: albertsmith
"""



import sys
import os
if 'google.colab' in sys.modules:
    _=help('modules')
    
    if not('MDAnalysis' in sys.modules):
        os.popen('pip3 install MDAnalysis')
    if not('pyDR' in sys.modules):
        os.popen('git clone https://github.com/alsinmr/pyDR.git')
    