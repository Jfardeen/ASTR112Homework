#!/usr/bin/env python
# coding: utf-8

# In[108]:


from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import numpy as np
import pandas as pd
from matplotlib.pyplot import figure
from shapely.geometry import Point
from shapely.ops import cascaded_union

get_ipython().run_line_magic('matplotlib', 'inline')


# In[162]:


RAW = pd.read_csv (r'C:\Users\jbfar\ASTR112\Homework1\data\stars.csv')
ParallaxSN_Table = pd.DataFrame(RAW, columns= ['parallax','parallax_error'] )

#Luminosity Temperature
LumTemp = pd.DataFrame(RAW, columns= ['teff_val', 'lum_val'])
LumTemp.dropna(subset = ["lum_val"], inplace=True)

#Color Magnitude
ColMag = pd.DataFrame(RAW, columns= ['parallax', 'phot_g_mean_mag', 'bp_rp', 'lum_val'])
ColMag.dropna(subset = ["lum_val"], inplace=True)
print(ParallaxSN_Table)


# In[3]:


#Division
ParallaxSN = ParallaxSN_Table['parallax']/ParallaxSN_Table['parallax_error']
print(ParallaxSN)


# In[140]:


#Histogram

plt.hist(ParallaxSN, bins=[10, 210, 410, 610, 810, 1010, 1210, 1410, 1610, 1810, 2010], density=False, edgecolor='k')
plt.xlabel('Parallax S/N')
plt.ylabel('Number of Stars')
plt.title('Histogram of Parallax S/N')
plt.xticks([10, 210, 410, 610, 810, 1010, 1210, 1410, 1610, 1810, 2010])
plt.figure(figsize=[20,20])


# In[170]:


x1 = LumTemp['teff_val']
y1 = LumTemp['lum_val']
colorspread = np.random.randint(0,2,2242)
colormap = np.array(['r', 'b'])
plt.scatter(x1, y1, alpha=0.3, s=1, c=colormap[colorspread])
plt.axis([max(x1)+100, min(x1)-100, min(y1), max(y1)])
plt.xlabel('Logarithm of Temperature (Kelvin)')
plt.ylabel('Logarithm of Luminosity (solLum)')
plt.yscale('log')
plt.xscale('log')
plt.title('Hertzsprung-Russel Diagram')
#lower white space (general goal)
#utilize range (max-min)
#certain modules allow me to have padding
#look into plotting on log scale - plot on linear scale
#look into plotting on a linear scale because of weirdness with units (minor)
plt.show


# In[168]:


x2 = ColMag['bp_rp']
y2 = ColMag['phot_g_mean_mag']
plt.scatter(x2,y2, s=1, alpha=0.3, c=colormap[colorspread])
plt.axis([min(x2)-0.1, max(x2)+0.1, max(y2)+0.25, min(y2)-0.25])
plt.xlabel('BP-RP Color (mag)')
plt.ylabel('Absolute G Band (mag)')
plt.title('Color-Magnitude Diagram')
plt.show


# In[172]:


Par = RAW['parallax']
lumi = RAW['lum_val']


# In[176]:


# M = m - 2.5log[(d/10)^2]
M = RAW['bp_rp'] - 2.5*np.log(np.power(np.power(RAW['parallax'], -1)*0.1, 2))
print(M)
#Absolute magnitude


# In[ ]:




