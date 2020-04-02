#=========================================================================================
# jetHTplots.py --------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# Author(s): Brendan Regnery -------------------------------------------------------------
#-----------------------------------------------------------------------------------------

# modules
import ROOT as root
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg') #prevents opening displays (fast), must use before pyplot
import matplotlib.pyplot as plt
import mplhep as hep
import uproot
#import uproot_methods
#import awkward

# enter batch mode in root (so python can access displays)
root.gROOT.SetBatch(True)

#=========================================================================================
# Load Data ------------------------------------------------------------------------------
#=========================================================================================

upTree = uproot.open("jetHT2017.root")["run/jetTree"]
jetDF = upTree.pandas.df(flatten=False)
print jetDF["HT"].count()

#=========================================================================================
# Make Plots -----------------------------------------------------------------------------
#=========================================================================================

# use the CMS plot style
plt.style.use(hep.style.ROOT)

# plot HT
plt.hist(jetDF["HT"], bins=30)
plt.ylabel('Events')
plt.xlabel(r'$H_{T}$ [GeV]')
hep.cms.cmslabel(data=True, paper=False, year='2017')
plt.savefig('plots/hist_HT.png')
plt.clf()



