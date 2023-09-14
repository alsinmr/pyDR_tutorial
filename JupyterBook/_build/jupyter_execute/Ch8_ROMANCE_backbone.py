#!/usr/bin/env python
# coding: utf-8

# # <font color="maroon">Chapter 8: ROMANCE Analysis of Backbone Motion</font>

# <a href="https://githubtocolab.com/alsinmr/pyDR_tutorial/blob/main/ColabNotebooks/Ch8_ROMANCE_backbone.ipynb" target="_blank"><img src="https://colab.research.google.com/assets/colab-badge.svg"></a>

# Some text about [ROMANCE](https://doi.org/10.1016/j.jmro.2022.100045) analysis.

# ## Setup and data downloads
# Since we've learned now how pyDR is organized and allows us to manage larger data sets, we'll now use the full project functionality.

# In[ ]:


# SETUP pyDR
import os
os.chdir('../..')


# In[17]:


#Imports
import numpy as np
import pyDR


# In[7]:


# Project Creation and File loading
proj=pyDR.Project()

sel=pyDR.MolSelect(topo='pyDR/examples/HETs15N/backboneB.pdb',
                   traj_files='pyDR/examples/HETs15N/backboneB.xtc',
                   project=proj)  #Selection object

# Specify the bond select to analyze for MD
sel.select_bond('N')

sel.traj.step=50  #Take every tenth point for MD calculation (set to 1 for more accurate calculation)
pyDR.Defaults['ProgressBar']=False #Turns of the Progress bar (screws up webpage compilation)


# ## Run the frame analysis

# ### Setup the frame analysis

# In[8]:


from numpy import nan
bsheets=[(226,228),(230,234),(236,240),(242,245),
         (262,264),(266,270),(272,276),(278,281)]

# Define second reference frame (beta sheets)
sel0=[sel.uni.select_atoms(f'resid {bsheet[0]}-{bsheet[1]}') for bsheet in bsheets]

fo=pyDR.Frames.FrameObj(sel)
fo.tensor_frame(sel1=1,sel2=2)  #Define the frame of the interaction (dipole tensor)
fo.new_frame(Type='peptide_plane',resids=sel.label) #Define first reference frame
fo.new_frame(Type='superimpose',sel=sel0) #Define second reference frame


# ### Import the frames, process

# In[9]:


fo.frames2data()  #Process frames, send to data
proj['Frames'].detect.r_auto(6)  #Prepare detectors
proj['Frames'].fit().opt2dist(rhoz_cleanup=True)  #Fit and optimize fit


# ### Plot the results

# First, we validate the frame analysis quality by comparing the product of the individual correlation functions to the total correlation function

# In[12]:


proj.close_fig('all')
proj['opt_fit']['Frames']['Direct'].plot(style='bar')
proj['opt_fit']['Frames']['Product'].plot().fig.set_size_inches([8,12])


# Plot the results for the Principal Axis System to the peptide plane (H–N librational motion)

# In[13]:


proj.close_fig('all')
proj['opt_fit']['PAS>peptide_plane'].plot(style='bar').fig.set_size_inches([8,12])
for a in proj.plot_obj.ax[:-1]:a.set_ylim([0,.05])


# Plot the results for peptide plane to $\beta$-sheet superimposition (peptide plane dynamics)

# In[14]:


proj.close_fig('all')
proj['opt_fit']['peptide_plane>superimpose'].plot(style='bar').fig.set_size_inches([8,12])
for a in proj.plot_obj.ax[:-1]:a.set_ylim([0,.2])


# Plot the results for $\beta$-sheet superimposition to lab frame ($\beta$-sheet dynamics

# In[21]:


proj.close_fig('all')
# Just show the beta sheets via index
resids=proj['opt_fit']['Direct'].label
index=np.zeros(len(resids),dtype=bool)
for bsheet in bsheets:
    index[np.argwhere(bsheet[0]==resids)[0,0]:np.argwhere(bsheet[1]==resids)[0,0]+1]=1
proj['opt_fit']['superimpose>LF'].plot(style='bar',index=index).fig.set_size_inches([8,12])
for a in proj.plot_obj.ax[:-1]:a.set_ylim([0,.05])


# In[ ]:




