#!/usr/bin/env python
# coding: utf-8

# # <font color="maroon">Chapter 3: A Basic NMR Analysis</font>

# <a href="https://githubtocolab.com/alsinmr/pyDR_tutorial/blob/main/JupyterBook/Ch3_basic_NMR.ipynb" target="_blank"><img src="https://colab.research.google.com/assets/colab-badge.svg"></a>

# In this example, we load an NMR data set and demonstrate how to go about fitting it. For this, we need a text file with the measured relaxation rates and experimental data in it. This is provided for this example, but one can also upload one's own data. Make sure to follow the prescribed file format. Entries are separated by tabs within a given line and by carriage returns over multiple lines. An example data file is printed out as an example.

# In[1]:


#Setup: Make sure packages are installed in google colab. 
#Also, make sure that we start in the main folder of the tutorial
import sys
import os
cwd=os.getcwd()
if 'google.colab' in sys.modules:
    get_ipython().system('git clone https://github.com/alsinmr/pyDR_tutorial.git')
    from pyDR_tutorial import colab_setup
else:
    os.chdir(os.path.join(cwd,'..')) #Go to the JupyterBook directory
    sys.path.append(os.path.join(cwd,'../..')) #Path to pyDR (make sure this is correct if working locally)


# In[2]:


#Imports
import pyDR


# ## Loading NMR data
# 
# First, we load a data object and discuss some of its contents, including plotting the experimental data.

# In[3]:


#Read-out data text file
with open('data/HETs_15N.txt','r') as f:
    for line in f:
        print(line.strip())


# In[4]:


data=pyDR.IO.readNMR('data/HETs_15N.txt')


# First, we discuss the data object briefly to understand the contents of the data object. There are a few key components:
# * Relaxation rate constants: **data.R**
# * Standard deviation of the relaxation rate constants: **data.Rstd**
# * Order parameters, given as $S^2$ (these are optional): **data.S2**
# * Standard deviation of $S^2$: **data.S2std**
# * Labels for the data: **data.label**
# * Sensitivity of the data: **data.sens**
# * Detector object for processing the data: **data.detect**
# * Experimental parameters: **data.info**
#     * This is a reference to data.sens.info
# * Details on the data / processing so far: **data.details**
# * Selection object, which connects data to a pdb of the molecular structure: **data.select**
#     * Not required: we'll connect the data to a pdb later
# 
# A number of other functions and data objects are also attached to data, but these are for data processing, or analysis of the processing after an initial fit.
# 
# Below, you can investigate these different components. We show data.info, which provides the various relevant experimental parameters for determining the sensitivities of the 8 experiments. Note that this data set also includes $S^2$, which is treated separately and does not show up in data.info.

# In[5]:


data.info


# The parameters found by default in data.info (can change if user provides their own sensitivity functions). Not all parameters are required for all experiments.
# * v0: External magnetic field, given as the $^1$H frequency in MHz
# * v1: Spin-locking strength for $R_{1\rho}$ experiments (applied to the relaxing nucleus), given in kHz
# * vr: Magic angle spinning frequency, given in kHz
# * offset: Offset applied to the spin-lock, given in kHz
# * stdev: Median standard deviation for data set
# * med_val: Median relaxation rate constant for data set
# * Nuc: Nucleus being relaxed (atomic number, symbol)
# * Nuc1: Nucleus (or list of nuclei) relaxing Nuc via dipole coupling
# * dXY: Dipole coupling (or list of couplings), given in Hz Same length as Nuc1
#     * Given as the full anisotropy of the dipole coupling, which is twice the coupling constant
# * CSA: Chemical shift anisotropy of nucleus being relaxed (z-component in ppm)
# * eta: Asymmetry of the CSA, which is unitless
# * CSoff: Unimplemented
# * QC: Quadrupole coupling in Hz
# * etaQ: Asymmetry of the quadrupole coupling
# * theta: Angle between the CSA and dipole couplings (for cross-correlated cross-relaxation)

# A graphical summary of the data object is obtained via data.plot (can use various plotting options, such as plot type, etc.). The plt_obj provides a variety of functions for manipulating the plot.

# In[6]:


plt_obj=data.plot(style='bar')
plt_obj.fig.set_size_inches([8,10])


# ## Processing NMR experimental data
# 
# Now, we want to process this data set. To do so, we first must optimize a set of detectors to analyze it with. The detector object is already attached to the data object, data.detect. So, we just need to run the optimization algorithm. We also will include $S^2$ in the data analysis, which we must also specify. We will use 4 detectors plus one more for $S^2$ data, to yield 5 total. Once the detector analysis is run, we may simply run the "fit" function to obtain the detector analysis, and finally plot the results. Note that this is essentially the same analysis that appeared
# 
# A. A. Smith, M. Ernst, S. Riniker, B. H. Meier. [Localized and collective motions in HET-s(218-289) Fibrils from Combined NMR Relaxation and MD Simulation.](https://onlinelibrary.wiley.com/doi/full/10.1002/anie.201901929) Angew. Chem. Int. Ed. 2019, 58, 9383-9388.
# 
# although a few improvements have been made since the initial analysis.

# In[7]:


data.detect.r_auto(4).inclS2()  #Optimize the detectors
fit=data.fit()  #Fit the data

plt_obj=fit.plot(style='bar')
plt_obj.fig.set_size_inches([8,10])
plt_obj.show_tc()
_=plt_obj.ax[-1].set_xlabel('Residue')


# Finally, we can check the fit quality by comparing the back-calculated relaxation rate constants to the original experimental relaxation rate constants.

# In[8]:


fig=fit.plot_fit()[0].axes.figure
fig.set_size_inches([12,10])


# In[ ]:




