#!/usr/bin/env python
# coding: utf-8

# # <font color="maroon">Chapter 4: Comparing MD and NMR Results</font>

# <a href="https://githubtocolab.com/alsinmr/pyDR_tutorial/blob/main/ColabNotebooks/Ch4_MDvsNMR.ipynb" target="_blank"><img src="https://colab.research.google.com/assets/colab-badge.svg"></a>

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


# In[3]:


#Download MD files
_=pyDR.IO.download('https://drive.google.com/file/d/1xgp5_BVeCh6Weu4tnsl1ToRZR49eM5hW/view?usp=share_link','HETs_md.xtc')
_=pyDR.IO.download('https://drive.google.com/file/d/1qkgKhQ7QMS8PAWS1AxyoaZUYej2qm_xT/view?usp=share_link','HETs_md.pdb')


# ## Load and Process NMR
# We'll do this in just a few lines. See the previous examples for processing NMR

# In[4]:


proj=pyDR.Project()
r_nmr=pyDR.IO.readNMR('https://drive.google.com/file/d/1w_0cR6ykjL7xvdxU2W90fRXvZ8XfLFc3/view?usp=share_link')
r_nmr.select=pyDR.MolSelect(topo='HETs_md.pdb').select_bond('N',segids='B',resids=r_nmr.label)
r_nmr.detect.r_auto(4).inclS2()
proj.append_data(r_nmr)
p_nmr=r_nmr.fit()  #Fit automatically goes into project


# ## Load MD and Process MD

# ###  Calculate Correlation functions

# In[5]:


sel=pyDR.MolSelect(topo='HETs_md.pdb',traj_files='/Users/albertsmith/Documents/GitHub.nosync/pyDR/examples/HETs15N/backboneB.xtc',project=proj).select_bond('N',segids='B')
r_md=pyDR.md2data(sel) #Automatically goes into project since project defined above


# ### MD Preprocessing (fit with no_opt detectors)

# In[25]:


r_md.detect.r_no_opt(12)
n_md=d_md.fit()


# ### Load above data from file (skip the slow download and computational steps)
# To avoid the calculations/downloads in the previous two cells, just run the cell below to get partially-processed MD data

# In[62]:


proj=pyDR.Project('HETs_MDvNMR/')  #This is already saved on google drive
n_md=proj['no_opt']['MD'][0]
r_nmr,p_nmr=proj['NMR']


# ### Process and plot MD data
# *Not for comparison to NMR*

# In[63]:


n_md.detect.r_auto(7)
p_md=n_md.fit()
o_md=p_md.opt2dist(rhoz_cleanup=True)


# In[64]:


proj.close_fig('all')
o_md.plot(style='bar')
proj.plot_obj.fig.set_size_inches([8,12])


# In[33]:


o_md.nglview(1,scaling=25)


# ## Compare NMR results to MD results

# ### Optimize MD detectors to match NMR detectors

# In[65]:


target=p_nmr.sens.rhoz
target[0][99:]=0 #Zero out rho0 for better MD performance
from matplotlib import pyplot as plt

_=plt.plot(target.T)


# In[74]:


n_md.detect.r_target(target,n=12)

ax=p_nmr.sens.plot_rhoz()[0].axes
_=n_md.detect.plot_rhoz(ax=ax,color='black',linestyle=':',index=range(5))


# In[75]:


p_md=n_md.fit()
o_md=p_md.opt2dist()


# In[88]:


get_ipython().run_line_magic('matplotlib', 'notebook')
proj.close_fig('all')
p_nmr.plot(style='bar',rho_index=range(3),errorbars=True)
o_md.plot(style='p')
proj.plot_obj.fig.set_size_inches([8,12])
for a,yl in zip(proj.plot_obj.ax,[.4,.08,.08]):a.set_ylim([0,yl])

