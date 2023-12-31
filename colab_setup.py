#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 10:00:03 2023

@author: albertsmith
"""

import sys
import os

if 'google.colab' in sys.modules:
    try:
        import MDAnalysis
    except:
        os.popen('pip3 install MDAnalysis')
    
    os.chdir('/content')
    
    try:
        import pyDR
    except:
        os.popen('git clone https://github.com/alsinmr/pyDR.git')
    
    os.chdir('/content/pyDR_tutorial')
    
    
    # NGLviewer setup
    # os.popen('pip install -q nglview')
    # from google.colab import output
    # output.enable_custom_widget_manager()