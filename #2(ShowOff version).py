#!/usr/bin/env python
# coding: utf-8

# In[2]:


from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import numpy as np
import pandas as pd
from matplotlib.pyplot import figure
from shapely.geometry import Point
from shapely.ops import cascaded_union
from scipy.interpolate import interp2d
from scipy.stats import gaussian_kde
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u
from astropy.constants import c, G
Gaia.ROW_LIMIT = -1
plt.rcParams['font.size'] = 16
get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


RAW = pd.read_csv (r'C:\Users\jbfar\ASTR112\Homework3\data\stars.csv')
Parallax = pd.DataFrame(RAW, columns= ['parallax'])
Parallax = pd.DataFrame.dropna(Parallax)
print(Parallax)


# In[4]:


#a)
ParallaxPositive = Parallax[~(Parallax['parallax'] <= 4.5)]
TrueParallax = ParallaxPositive[~(ParallaxPositive['parallax'] >= 7)]
#Shrink the Parallax Dataframe to contain only targetted range

plt.figure(figsize=[20, 5])
bincount, binedges, axes = plt.hist(TrueParallax, bins=50, density=False, edgecolor='k', facecolor='r')
[axes[i].set_facecolor('b') for i in range (11,19)]
#Histogram

peak_bin = binedges[np.argmax(bincount)]
plt.axvline(peak_bin+.775, ls="--", color="r",lw=2,
    label="Peak at ~{:0.3f} mas".format(peak_bin+0.775))
plt.axvspan(peak_bin+.5, peak_bin+.9, color="r", alpha=0.2,
    label="Peak range")
#^Shaded region ahnd dotted line in histogram

plt.xlabel('Parallax (mas)')
plt.ylabel('Number of Occurrences')
plt.title('Histogram plot of the Parallaxes')
plt.legend()
#For Looks


# In[5]:


#here I reduce the dataframe to contain only values of my lowest bin and higher and highest bin and lower
NewFrame = Parallax[~(Parallax['parallax'] <= 5.049923562286532)]
NewerFrame = NewFrame[~(NewFrame['parallax'] >= 5.399865361322014)]
print(NewerFrame)
#I looked up the actual distance and plotted it


# In[19]:


#b)
# I need to redo the dataframe reducing code but with my original (RAW) dataframe
# I need the pmra (Proper Motion Right ascension) and pmdec (Proper Motion Declination) columns
RAW1 = RAW.dropna(subset=['parallax', 'pmra', 'pmdec'])
RAW2 = RAW1[~(RAW1['parallax'] <= 5.049923562286532)]
RAW3 = RAW2[~(RAW2['parallax'] >= 5.399865361322014)]
pmra = pd.DataFrame(RAW3, columns= ['pmra'])
pmdec = pd.DataFrame(RAW3, columns= ['pmdec'])
pmra_new = pmra.dropna(subset=['pmra'])
pmdec_new = pmdec.dropna(subset=['pmdec'])
# print(pmra_new)
# print(pmdec_new)
# print(RAW3)


# In[7]:


plt.figure(figsize=[10, 10])
plt.hist2d(pmdec_new['pmdec'], pmra_new['pmra'], bins=2500)
plt.show()


# In[8]:


plt.figure(figsize=[15, 15])
plt.hist2d(pmdec_new['pmdec'], pmra_new['pmra'], bins=2500)
plt.xlim(-17,-7)
plt.ylim(-40,-30)
plt.xlabel('Proper motion in Declination Direction (mas.yr**-1)')
plt.ylabel('Proper motion in Right Ascension Direction (mas.yr**-1)')
plt.title('2D Histogram of the Star Cluster')


# In[9]:


#Range for pmdec is (-14, -11)
#Range for pmra is (-37, -34)
#A circle of radius 1.5 centered at (-12.5,-35.5)


# In[18]:


#c
RAWRef1 = RAW1[~(RAW1['pmdec'] <= -15)]
RAWRef2 = RAWRef1[~(RAWRef1['pmdec'] >= -10)]
RawRef3 = RAWRef2[~(RAWRef2['pmra'] <= -38)]
RAWRef = RawRef3[~(RawRef3['pmra'] >= -33)]
#Creating Dataframe where pmra and pmdec are within certain values
RAWNoClust = RAW1[~RAW1.index.isin(RAWRef.index)]
# print(RAWNoClust)
# print(RAWRef)


# In[16]:


plt.figure(figsize=[10, 10])
xclus = RAWRef['pmdec']
yclus = RAWRef['pmra']

xnoclus = RAWNoClust['pmdec']
ynoclus = RAWNoClust['pmra']

x = RAW1['pmdec']
y = RAW1['pmra']

R = 1.5; mask = np.sqrt((x+12.5) * (x+12.5) + (y+35.5) * (y+35.5)) < R
# plt.plot(xclus, yclus, 'o', markersize=1, color='k')
# plt.plot(xnoclus, ynoclus, 'o', markersize=0.1, color='r')
plt.plot(x[mask], y[mask], 'o', color='k', markersize=1)
plt.plot(x[~mask], y[~mask], 'o', color='r', markersize=0.1)
plt.xlim(-20,0)
plt.ylim(-40,0)
plt.xlabel('Proper motion in Declination Direction (mas.yr**-1)')
plt.ylabel('Proper motion in Right Ascension Direction (mas.yr**-1)')
plt.title('Scatterplot where the black region is the star cluster')
plt.show()


# In[17]:


# print(RAWRef)
RAWRef = RAWRef.dropna(subset=['phot_g_mean_mag', 'bp_rp'])

