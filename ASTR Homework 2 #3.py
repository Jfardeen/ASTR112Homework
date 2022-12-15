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

get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


RAW = pd.read_csv (r'C:\Users\jbfar\ASTR112\Homework2\data\stars.csv')
Parallax = pd.DataFrame(RAW, columns= ['parallax'])
max = Parallax.max()
min = Parallax.min()
print(max)


# In[4]:


#Histogram

plt.hist(Parallax, bins=100, density=False, edgecolor='k', label='test')
plt.figtext(0.17, 0.55, 'peak')
plt.xlim(0.5, 22.1)
plt.xlabel('Parallax (mas)')
plt.ylabel('Number of Stars')
plt.title('Histogram of the Parallax')
plt.figure(figsize=[20,20])
#peak from 0-1.2 mas


# In[5]:


#c)
#The parallax range is  0.5-1.36 mas - or 0.0005-0.00136 as
#To convert to distance in parsecs - use Distance = 1/(parallax)
#This makes the range 735.29pc-2000pc


# In[6]:


RAW2 = pd.read_csv (r'C:\Users\jbfar\ASTR112\Homework2\data\stars2.csv')
CMD = pd.DataFrame(RAW2, columns=['phot_g_mean_mag', 'bp_rp'])
CMD.dropna(subset = ["bp_rp"], inplace=True)
print(CMD)


# In[11]:


plt.figure(figsize=[10,10])
x2 = CMD['bp_rp']
y2 = CMD['phot_g_mean_mag']
colormap = np.array(['r', 'b'])
colorspread = np.random.randint(0,2,3694)
plt.scatter(x2,y2, s=5, alpha=0.3, c=colormap[colorspread])
plt.xlabel('BP-RP Color (mag)')
plt.ylabel('Absolute G Band (mag)')
plt.title('Color-Magnitude Diagram')


# In[8]:


#e)
#Potentially I could use their relative distance from one another as a constraint.
#The stars in the cluster will be clustered together and their relative distance from one another will be quite small
#compared to that of a star not in the cluster.


# In[ ]:




