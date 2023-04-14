import uproot

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

events1 = uproot.open("~/tafc/Tafc_Project/HeatMap1.root:HeatMap1;1")
Data = events1.arrays(['photon_z', 'photon_phi'], library="pd")

Data=Data.head(30000000)
my_cmap = plt.get_cmap('Spectral')
sns.histplot(Data, x=Data["photon_phi"], y=Data["photon_z"], binwidth=(0.5,0.6), cmap = my_cmap, cbar=True)
plt.ylim(0.02, 0.08)

plt.show()
