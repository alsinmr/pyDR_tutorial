#!/usr/bin/env python
# coding: utf-8

# In[1]:


import nglview as nv
from nglview.color import ColormakerRegistry as cm

func_str = '''
             this.atomColor = function (atom) {
             if (atom.serial < 1000) {
               return 0x0000FF // blue
             } else if (atom.serial > 2000) {
               return 0xFF0000 // red
             } else {
               return 0x00FF00 // green
             }
             }
         '''
cm.add_scheme_func('CustomColor',func_str)

func_str = '''
             this.atom.radius = function (atom) {
             if (atom.serial < 1000) {
               return .1 // blue
             } else if (atom.serial > 2000) {
               return .5 // red
             } else {
               return 1 // green
             }
             }
         '''
cm.add_scheme_func('CustomRadius',func_str)


# <a href="https://githubtocolab.com/alsinmr/pyDR_tutorial/blob/main/ColabNotebooks/Untitled.ipynb" target="_blank"><img src="https://colab.research.google.com/assets/colab-badge.svg"></a>

# In[2]:


import sys
sys.path.append('..')
import pyDR
data=pyDR.IO.readNMR('data/HETs_15N.txt')


# In[3]:


data.select=pyDR.MolSelect(topo='2kj3').select_bond('N',segids='B',resids=data.label)


# In[4]:


v1=nv.show_mdanalysis(data.select.uni,representations=[{"type":"ball+stick","params":{"color":"CustomColor","radius":1}}])
v1


# In[5]:


v1.shape.add_sphere(pos,data.select.uni.atoms[0].position,5)


# In[ ]:


v1


# In[ ]:


# v1.picked['atom1']['vdw']=6
v1.picked


# In[ ]:


data.select.uni.select_atoms('resid 292 and segid B and name HE1').positions


# In[ ]:


pos=[v1.picked['atom1'][key] for key in ['x','y','z']]
print(pos)


# In[ ]:



func_str = '''
             this.atom.radius = function (atom) {
             if (atom.serial < 1000) {
               return .1 
             } else if (atom.serial > 2000) {
               return 1 
             } else {
               return 5 
             }
             }
         '''
cm.add_scheme_func('radius',func_str)


# In[ ]:


v1._representations


# In[ ]:


v1._js()


# In[ ]:


v1._js('print(atom)')


# In[ ]:


string="""
var schemeID = NGL.atom.radius.addScheme(function (params){
    this.atom.radius = function (atom) {
             if (atom.serial < 1000) {
               return .1 
             } else if (atom.serial > 2000) {
               return 1 
             } else {
               return 5 
             }
             }
})

this._updateId(schemeId, 'radius')
"""


# In[ ]:


v1._js(string)


# In[ ]:


cm._js(string)


# In[ ]:


v1.get_state()['_ngl_msg_archive'].keys()


# In[ ]:


v1.


# In[ ]:


v1.parameters


# In[ ]:


v1.camera


# In[ ]:


v1._representations


# In[ ]:


v1._

