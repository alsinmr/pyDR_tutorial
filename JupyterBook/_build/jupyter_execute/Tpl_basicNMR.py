#!/usr/bin/env python
# coding: utf-8

# # <font color="maroon">Template: Basic NMR</font>

# <a href="https://githubtocolab.com/alsinmr/pyDR_tutorial/blob/main/ColabNotebooks/Tpl_basicNMR.ipynb" target="_blank"><img src="https://colab.research.google.com/assets/colab-badge.svg"></a>

# To use this template, you need to prepare a text file with NMR relaxation data in it (see [HETs_15N.txt](https://github.com/alsinmr/pyDR_tutorial/raw/main/data/HETs_15N.txt) or [ubi_soln.txt](https://github.com/alsinmr/pyDR_tutorial/raw/main/data/ubi_soln.txt) for example). If you intend to run on a local pyDR installation, then you just need to point pyDR to the file. If you want to run in Google Colab, the file needs to be available online somehow. The suggested options are in your Google drive, where we will mount Google drive in the notebook, or via a shareable weblink (for example, also available in Dropbox, Google Drive, etc.).

# In[1]:


# SETUP pyDR
import os
os.chdir('..')
import sys
sys.path.append('../') # Path to pyDR location


# In[2]:


#Imports
import pyDR


# ## Load NMR Data
# The best way to get load your data depends if you're running locally (just point to the file) or if you're running online (download from somewhere, mount Google Drive)
# 
#  Note that you can use commands such as 'ls', 'cd', and 'pwd' if you're confused where your files are.

# ### Loading data v1: Running locally
# 
# Just point to file

# In[3]:


data=pyDR.IO.readNMR('data/HETs_15N.txt')  #Change path to your file


# ### Loading data v2: Download from online source
# Just point to a url online. Note: except for Google Drive, the url needs to be the actual file and not a html viewer showing the file. For example, this [GitHub webpage](https://github.com/alsinmr/pyDR_tutorial/blob/main/data/HETs_15N.txt) shows our HET-s data. However, it is a viewer for the data, and not directly useable as a text file. However, the button in the upper right of the data viewer that says "raw" can be used to obtain the actual [text file](https://github.com/alsinmr/pyDR_tutorial/raw/main/data/HETs_15N.txt).
# 
# For Google Drive, pyDR corrects the link internally, so you can just provide the share link. Note that without mounting Google drive, the file must be viewable to anyone with the link.

# In[4]:


# Example: raw file from Github
data=pyDR.IO.readNMR('https://github.com/alsinmr/pyDR_tutorial/raw/main/data/HETs_15N.txt')


# In[5]:


#Example: Google Drive share link
data=pyDR.IO.readNMR('https://drive.google.com/file/d/1w_0cR6ykjL7xvdxU2W90fRXvZ8XfLFc3/view?usp=share_link')


# ### Loading data v3: Colab with mounted drive
# Mounting Google Drive will first cause a prompt in this window, followed by a popup window for logging in. An advantage here is you can save results/figures to Google Drive.

# In[6]:


if 'google.colab' in sys.modules:
    from google.colab import drive
    drive.mount('/content/drive')
    data=pyDR.IO.readNMR(os.path.join('/content/drive/MyDrive','path_to_data'))


# ## Put data into a project
# 
# Projects are convenient ways to manage a lot of data, and provide convenient tools for overlaying data in 2D plots, as well as visualizing data in 3D in ChimeraX (ChimeraX doesn't work on Colab).

# In[7]:


proj=pyDR.Project(directory=None)    #Include a directory to save the project
proj.append_data(data)


# ## Attach structure to the data
# Download a pdb and attach it to the data object. For protein backbone dynamics, we recommend using the labels to provide the residue number, making this step easier.
# 
# Available bonds ('Nuc')
# * N,15N,N15       : Backbone N and the directly bonded hydrogen 
# * C,CO,13CO,CO13  : Backbone carbonyl carbon and the carbonyl oxygen
# * CA,13CA,CA13    : Backbone CA and the directly bonded hydrogen (only HA1 for glycine)
# * CACB            : Backbone CA and CB (not usually relaxation relevant)
# * IVL/IVLA/CH3    : Methyl groups in Isoleucine/Valine/Leucine, or ILV+Alanine, or simply all methyl groups. Each methyl group returns 3 pairs, corresponding to each hydrogen
# * IVL1/IVLA1/CH31 : Same as above, except only one pair
# * IVLl/IVLAl/CH3l : Same as above, but with only the 'left' leucine and valine methyl group
# * IVLr/IVLAr/CH3r : Same as above, but selects the 'right' methyl group
# * FY_d,FY_e,FY_z  : Phenylalanine and Tyrosine H–C pairs at either the delta, epsilon, or  zeta positions.
# * FY_d1,FY_e1,FY_z1:Same as above, but only one pair returned for each amino acid
# 
# We can also filter based on residues, segments, and a filter string ([MDAnalysis](https://docs.mdanalysis.org/stable/documentation_pages/selections.html) format).

# In[8]:


data.select=pyDR.MolSelect(topo='2KJ3')
data.select.select_bond(Nuc='N',resids=data.label,segids='B')


# ## Plot the data

# In[9]:


plt_obj=data.plot(style='bar')
plt_obj.fig.set_size_inches([8,10])


# ## Process NMR data

# In[10]:


data.detect.r_auto(4)    #Set number of detectors here
data.detect.inclS2() #Uncomment if including order parameters

fit=data.fit()  #Fit the data


# ## Plot the results

# In[11]:


proj.close_fig('all')
plt_obj=fit.plot(style='bar')
plt_obj.fig.set_size_inches([8,10])
plt_obj.show_tc()
_=plt_obj.ax[-1].set_xlabel('Residue')


# ## Plot the fit quality

# In[12]:


fig=fit.plot_fit()[0].axes.figure
fig.set_size_inches([12,10])


# ## Visualize with NGL viewer

# In[13]:


fit.nglview(0)

