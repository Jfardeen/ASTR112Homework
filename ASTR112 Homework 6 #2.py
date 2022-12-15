#!/usr/bin/env python
# coding: utf-8

# In[2]:


from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import numpy as np
import pandas as pd
from matplotlib.pyplot import figure
import scipy.stats
from shapely.geometry import Point
from shapely.ops import cascaded_union
from scipy.interpolate import interp2d
from scipy.stats import gaussian_kde
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u
from astropy.constants import c, G
import scipy.optimize as opt
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib.colors import LogNorm
from matplotlib.pyplot import *
import turtle
Gaia.ROW_LIMIT = -1
plt.rcParams['font.size'] = 16
get_ipython().run_line_magic('matplotlib', 'inline')


# # GAIA QUERY
# 
# SELECT gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.a_g_val
# FROM gaiadr2.gaia_source 
# WHERE (gaiadr2.gaia_source.parallax>=10 AND gaiadr2.gaia_source.parallax_over_error>=50 AND gaiadr2.gaia_source.phot_g_mean_flux_over_error>=50 AND gaiadr2.gaia_source.phot_bp_mean_flux_over_error>=5 AND gaiadr2.gaia_source.phot_rp_mean_flux_over_error>=5)

# In[3]:


RAW = pd.read_csv(r'C:\Users\jbfar\ASTR112\Homework6\data\stars.csv')


# In[4]:


#Removing NaN
RAW = RAW.dropna(subset=['phot_g_mean_mag'])
RAW = RAW.dropna(subset=['bp_rp'])
print(RAW)


# In[87]:


#Making 1D Histograms
n2, bins2, patches2 = plt.hist(RAW['phot_g_mean_mag'], bins=252378)
print(n2)
z2 = np.array(n2)


# In[86]:


n1, bins1, patches1 = plt.hist(RAW['bp_rp'], bins=252378)
print(n1)
z1 = np.array(n1)


# In[36]:


#Normal Scatterplot
X = RAW['bp_rp']
Y = RAW['phot_g_mean_mag']
plt.figure(figsize=[12,7])
plt.scatter(X, Y, alpha=0.1)
plt.axis([min(X),max(X),max(Y),min(Y)])


# In[34]:


#2D Histogram
plt.figure(figsize=[20,15])
counts, xbins, ybins, image = plt.hist2d(X, Y, bins=2000)
plt.axis([min(X),max(X),max(Y),min(Y)])
plt.show()


# In[ ]:





# In[33]:


#An ABSOLUTELY BEUATIFUL CONTOUR PLOT
fig, ax = plt.subplots(figsize=[20,15])
contours = ax.contour(counts.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()])
plt.axis([min(X),max(X),max(Y),min(Y)])

plt.xlabel('BP-RP')
plt.ylabel('Absolute G Magnitude')
plt.title('CMD of Gaia DR2 Stars')
plt.text(1.5, 6, 'Red Giants', color='red')
plt.text(1.5, 5, 'Red Clump', color='green')
plt.text(0.4, 13, 'White Dwarf', color='blue')
plt.text(3, 10, 'Main Sequence', fontsize=20)

#Boxes highlighting specific regions
RedGiant = ptc.Rectangle((1, 2.5),0.5,5, fill = False, color = "red", linewidth = 2)
ax.add_patch(RedGiant)
RedClump = ptc.Rectangle((1.05, 4), 0.4, 2, color='green', fill=False, linewidth = 2)
ax.add_patch(RedClump)
WhiteDwarf = ptc.Rectangle((-0.5,13.5), 2, 5.5, color='blue', fill=False, linewidth = 2)
ax.add_patch(WhiteDwarf)
plt.show()

#Actually the best plot manking has ever seen >:)


# In[ ]:





# In[ ]:





# In[ ]:




