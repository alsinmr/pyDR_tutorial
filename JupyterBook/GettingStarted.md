# <font color="maroon">Getting Started</font>

There are currently two possiblities for running pyDR: with a local Python installation, or online via [Google Colab](https://colab.research.google.com/). The local installation provides the most functionality, including 3D visualization with [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) and in-browser visualization with [NGL viewer](https://nglviewer.org/). NGL viewer should in principle work in Google Colab, but at the moment fails for unknown reasons. A local installation of pyDR may be run in [Jupyter notebooks](https://jupyter.org/), as a script, or in a Python/iPython terminal (pyDR has been written to provide some extra functionality in iPython). Here we provide a few notes for a local installation or for running in Google Colab.


## Local Installation
pyDR does not currently exist for installation via [pip](https://pypi.org/) or [conda](https://docs.conda.io/en/latest/), and instead is downloaded from GitHub, with its containing folder added to the system path (or placed with other Python modules). The required modules must be installed separately by the user.

[**pyDR on GitHub**](https://github.com/alsinmr/pyDR/)

Versions listed tested for pyDR, although most recent module versions should work.

### Requirements
* Python 3 (3.7.3)
* numpy (1.19.2)
* matplotlib (3.4.2) 
* scipy (1.5.2)
* [MDanalysis](https://mdanalysis.org/) (1.0.0)

### Recommended Installations
* [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) (1.5): 3D visualization for detectors (this is a separate program– not installed via pip or conda)
* [pyfftw](https://pypi.org/project/pyFFTW/) (0.12.0): Faster Fourier transforms for correlation function calculation
* [Jupyter notebooks](https://jupyter.org/): Neat code organization based in a web browser
* [NGL viewer](https://nglviewer.org/) (3.0.3): 3D visualization in Jupyter notebooks

If you're starting from scratch, you may consider installing [Anaconda](https://anaconda.org), which will install Python 3, numpy, matplotlib, scipy, and Jupyter Notebooks for you. The remaining packages (MDAnalysis, etc.) may then be installed with conda (or less ideally, with pip). ChimeraX is not installed with Python, but has its own installer (you will need to tell pyDR where ChimeraX is installed, use pyDR.chimeraX.chimeraX_funs.set_chimera_path).

Once set up, any of the webpages with code in the tutorial can be downloaded with the button in the upper right (pick .ipynb) and run locally. 

## Google Colab
Colab has the advantage that there is no local installation requirement, with the base requirements for Python already setup. The notebooks that we provide include a few lines at the beginning that install the additional requirements for MDAnalysis and NGL viewer, so you don't have to know how to set it up– just execute the cells (shift+enter executes a cell).

For all notebooks in the tutorial, there is a button near the top that redirects you to Google Colab, where you can then edit and execute contents of the notebook: ![](https://colab.research.google.com/assets/colab-badge.svg)


Google Colab, however, is a little buggy at the moment. We hope for continued improvement, but have no real control over this. There are two major bugs that we observe.

1. We start by cloning the pyDR tutorial in GitHub, followed by importing 'colab_setup' from the pyDR tutorial. These two commands download and install all the required files for importing pyDR. Although the cell finishes running, it is often the case that the installation is not actually finished. Then, the import of pyDR fails. Usually, one simply needs to wait a minute or so until pyDR will import.

2. While visualization in NGLview is in principle supported in Google Colab, it is often the case that running data.nglview will either fail to return the NGLview widget, or will return a white molecule instead of a colored molecule in the widget. 