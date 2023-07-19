#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.stats import mad_std
from time import sleep
import numpy.random as ran
import scipy.optimize as opt
from scipy.optimize import curve_fit
import scipy.stats as stat
import scipy as scipy
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
import matplotlib.ticker as ticker

get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


dataFrame_01 = pd.read_csv (r'C:\Users\jbfar\ASTR135\Lab2\data\TRACE001.CSV', skiprows=11)
dataFrame_02 = pd.read_csv (r'C:\Users\jbfar\ASTR135\Lab2\data\TRACE007.CSV', skiprows=11)


# In[4]:


#Convert Frequency to Proper Velocity
dataFrame_01['Prop. Velocity(km/s)'] = (3*10**5)*(1-(dataFrame_01['Frequency(Hz)'])/(1000000*1420.4))
dataFrame_02['Prop. Velocity(km/s)'] = (3*10**5)*(1-(dataFrame_02['Frequency(Hz)'])/(1000000*1420.4))


# In[5]:


#Convert Amplitude to (Watts)
dataFrame_01['Amplitude(W)'] = 10**((dataFrame_01['Amplitude(dBm)']-30)/10)
dataFrame_02['Amplitude(W)'] = 10**((dataFrame_02['Amplitude(dBm)']-30)/10)


# In[6]:


#Convert to Watts/Hertz
dataFrame_01['Amplitude(W/Hz)'] = (1/30000)*dataFrame_01['Amplitude(W)']
dataFrame_02['Amplitude(W/Hz)'] = (1/30000)*dataFrame_02['Amplitude(W)']


# In[7]:


#Account for Antenna
dataFrame_01['Amplitude(W/Hz)at'] = (1.995/5011.872)*dataFrame_01['Amplitude(W/Hz)']
dataFrame_02['Amplitude(W/Hz)at'] = (19.953/5011.872)*dataFrame_02['Amplitude(W/Hz)']


# In[8]:


#Polarize
dataFrame_01['Amplitude(W/Hz)polar'] = 2*dataFrame_01['Amplitude(W/Hz)at']
dataFrame_02['Amplitude(W/Hz)polar'] = 2*dataFrame_02['Amplitude(W/Hz)at']


# In[9]:


#Trim Dataframe
dataFrame_01 = dataFrame_01[~(dataFrame_01['Prop. Velocity(km/s)'] >= 100)]
dataFrame_02 = dataFrame_02[~(dataFrame_02['Prop. Velocity(km/s)'] >= 100)]


# In[10]:


#Define Gaussians
def gaussian(x, *p):
    return (((p[0]/(p[2]*np.sqrt(2*np.pi)))*np.exp(-(1/2)*(((x-p[1])/p[2])**2)))+p[3])


# In[11]:


def _4gaussian(x, *p):
    return (((((p[0]/(p[2]*np.sqrt(2*np.pi)))*np.exp(-(1/2)*(((x-p[1])/p[2])**2)))+p[3])) +             ((((p[4]/(p[6]*np.sqrt(2*np.pi)))*np.exp(-(1/2)*(((x-p[5])/p[6])**2)))+p[7])) +             ((((p[8]/(p[10]*np.sqrt(2*np.pi)))*np.exp(-(1/2)*(((x-p[9])/p[10])**2)))+p[11])) +             ((((p[12]/(p[14]*np.sqrt(2*np.pi)))*np.exp(-(1/2)*(((x-p[13])/p[14])**2)))+p[15])))


# In[12]:


#Standard Deviation
dataFrame_01std = dataFrame_01[~(dataFrame_01['Prop. Velocity(km/s)'] >= -100)]
dataFrame_01std = dataFrame_01std[~(dataFrame_01std['Prop. Velocity(km/s)'] <= -150)]
dumb, std = stat.norm.fit(dataFrame_01std['Prop. Velocity(km/s)'])
print(std)


# In[13]:


#More Dataframe Trimming
dataFrame_01params1 = dataFrame_01[~(dataFrame_01['Prop. Velocity(km/s)'] <= 22)]
dataFrame_01params1 = dataFrame_01params1[~(dataFrame_01params1['Prop. Velocity(km/s)'] >= 75)]

dataFrame_01params2 = dataFrame_01[~(dataFrame_01['Prop. Velocity(km/s)'] <= -10)]
dataFrame_01params2 = dataFrame_01params2[~(dataFrame_01params2['Prop. Velocity(km/s)'] >= 22)]

dataFrame_01params3 = dataFrame_01[~(dataFrame_01['Prop. Velocity(km/s)'] <= -40)]
dataFrame_01params3 = dataFrame_01params3[~(dataFrame_01params3['Prop. Velocity(km/s)'] >= -10)]

dataFrame_01params4 = dataFrame_01[~(dataFrame_01['Prop. Velocity(km/s)'] <= -100)]
dataFrame_01params4 = dataFrame_01params4[~(dataFrame_01params4['Prop. Velocity(km/s)'] >= -39)]


# In[26]:


#Scatterplot
print(dataFrame_01['Amplitude(W/Hz)polar'])
plt.figure(figsize=[20,5])
x = dataFrame_01['Prop. Velocity(km/s)']
y = dataFrame_01['Amplitude(W/Hz)polar']
plt.scatter(x,y)
plt.ylim(1e-21, 2e-21)
plt.xticks(np.arange(-100, 100, 5.0))
plt.grid()

#Close Parameter Guesses
guess_params1 = np.array([13e-19, 35, 8, 21e-20])
plt.plot(x, gaussian(x, *guess_params1))

guess_params2 = np.array([18e-19, 5, 9, 21e-20])
plt.plot(x, gaussian(x, *guess_params2))

guess_params3 = np.array([15.5e-19, -18, 8, 21e-20])
plt.plot(x, gaussian(x, *guess_params3))

guess_params4 = np.array([17e-19, -74, 14.5, 21e-20])
plt.plot(x, gaussian(x, *guess_params4))

#xsmooth1
xsmooth1 = np.linspace(32, 75, 339)
xsmooth1 = dataFrame_01params1['Prop. Velocity(km/s)']
ysmooth1 = dataFrame_01params1['Amplitude(W/Hz)polar']
popt1, pcov1 = opt.curve_fit(gaussian, xsmooth1, ysmooth1, p0=guess_params1)
print(popt1)
fsmooth1 = gaussian(x,*popt1)

#xsmooth2
xsmooth2 = dataFrame_01params2['Prop. Velocity(km/s)']
ysmooth2 = dataFrame_01params2['Amplitude(W/Hz)polar']
popt2, pcov2 = opt.curve_fit(gaussian, xsmooth2, ysmooth2, p0=guess_params2)
print(popt2)
fsmooth2 = gaussian(x,*popt2)

#xsmooth3
xsmooth3 = dataFrame_01params3['Prop. Velocity(km/s)']
ysmooth3 = dataFrame_01params3['Amplitude(W/Hz)polar']
popt3, pcov3 = opt.curve_fit(gaussian, xsmooth3, ysmooth3, p0=guess_params3)
print(popt3)
fsmooth3 = gaussian(x,*popt3)

#xsmooth4
xsmooth4 = dataFrame_01params4['Prop. Velocity(km/s)']
ysmooth4 = dataFrame_01params4['Amplitude(W/Hz)polar']
popt4, pcov4 = opt.curve_fit(gaussian, xsmooth4, ysmooth4, p0=guess_params4)
print(popt4)
fsmooth4 = gaussian(x,*popt4)
# plt.plot(x,fsmooth4,color='violet')

superpram = np.append(guess_params1, guess_params2)
ultrapram = np.append(guess_params3, guess_params4)
megapram = np.append(superpram, ultrapram)
print(megapram)

dataFrame_01trim = dataFrame_01[~(dataFrame_01['Prop. Velocity(km/s)'] <= -100)]
dataFrame_01trim = dataFrame_01trim[~(dataFrame_01trim['Prop. Velocity(km/s)'] >= 75)]
ultrax = dataFrame_01trim['Prop. Velocity(km/s)']
ultray = dataFrame_01trim['Amplitude(W/Hz)polar']
popt, pcov = opt.curve_fit(_4gaussian, ultrax, ultray, p0=megapram, maxfev = 20000)
print(popt)
fsmooth = _4gaussian(x, *popt)
plt.plot(x, _4gaussian(x, *popt), color='crimson')


# In[15]:


def _4gaussiannew(x, *p):
    return ((((p[0]/(p[2]*np.sqrt(2*np.pi)))*np.exp(-(1/2)*(((x-p[1])/p[2])**2)))+p[12]) + 
            (((p[3]/(p[5]*np.sqrt(2*np.pi)))*np.exp(-(1/2)*(((x-p[4])/p[5])**2)))+p[12]) + 
            (((p[6]/(p[8]*np.sqrt(2*np.pi)))*np.exp(-(1/2)*(((x-p[7])/p[8])**2)))+p[12]) + 
            (((p[9]/(p[11]*np.sqrt(2*np.pi)))*np.exp(-(1/2)*(((x-p[10])/p[11])**2)))+p[12]))
def _2gaussian(x, *p):
    return ((((p[0]/(p[2]*np.sqrt(2*np.pi)))*np.exp(-(1/2)*(((x-p[1])/p[2])**2)))+p[6]) + 
            (((p[3]/(p[5]*np.sqrt(2*np.pi)))*np.exp(-(1/2)*(((x-p[4])/p[5])**2)))+p[6]))
def _3gaussian(x, *p):
    return ((((p[0]/(p[2]*np.sqrt(2*np.pi)))*np.exp(-(1/2)*(((x-p[1])/p[2])**2)))+p[9]) + 
            (((p[3]/(p[5]*np.sqrt(2*np.pi)))*np.exp(-(1/2)*(((x-p[4])/p[5])**2)))+p[9]) + 
            (((p[6]/(p[8]*np.sqrt(2*np.pi)))*np.exp(-(1/2)*(((x-p[7])/p[8])**2)))+p[9]))
def linear(x, *p):
    return p[0]*x+p[1]


# In[16]:


print(dataFrame_02)


# In[17]:


plt.figure(figsize=[20,5])
x7 = dataFrame_02['Prop. Velocity(km/s)']
y7 = dataFrame_02['Amplitude(W/Hz)polar']
plt.scatter(x7,y7)
guess_params7 = np.array([8.5e-20, -7.4, 10, 1.39e-20])
#Original big peak gauss

dataFrame_02prams = np.array([2e-20, -85, 11, 2e-20, -35, 12, 6.95e-21])
pleasework = np.array([2e-20, -85, 11, 2e-20, -35, 12, 7.7e-20, -7.5, 9, 4.63e-21])
# small peak guesstimate

xsmooth7 = dataFrame_02['Prop. Velocity(km/s)']

plot7 = gaussian(xsmooth7, *guess_params7)

please = _3gaussian(xsmooth7, *pleasework)
plt.plot(xsmooth7, please, color='purple', label='Test Fit 3g')

#Curve Fitting Shenanigans
popt7, pcov7 = opt.curve_fit(gaussian, x7, y7, p0=guess_params7, maxfev = 20000)
print(popt7)
fsmooth7 = gaussian(xsmooth7, *popt7)

_2popt7, pcov7 = opt.curve_fit(_2gaussian, x7, y7, p0=dataFrame_02prams, maxfev = 20000)
print(_2popt7)
f2smooth7 = _2gaussian(xsmooth7, *_2popt7)

_3popt7, pcov7 = opt.curve_fit(_3gaussian, x7, y7, p0=pleasework, maxfev = 20000)
print(_3popt7)
f3smooth7 = _3gaussian(xsmooth7, *_3popt7)
plt.plot(xsmooth7, _3gaussian(xsmooth7, *_3popt7), color='crimson', label='Real Fit')
plt.legend()


# In[20]:


dataFrame_02['AmpSubBack'] = dataFrame_02['Amplitude(W/Hz)polar'] - linfsmooth7
plt.figure(figsize=[20,5])
plt.scatter(x7, dataFrame_02['AmpSubBack'])
dataFrame_02['Brightness Distribution'] = (dataFrame_02['Amplitude(W/Hz)polar']/(2*1.38e-23*0.5))
dataFrame_02['Brightness Distribution2'] = (dataFrame_02['AmpSubBack']/(2*1.38e-23*0.5))
plt.xlabel('Radial Velocity (km/s)')
plt.ylabel('Power (W/Hz)')
plt.title('Power vs. Radial Velocity')


# In[ ]:





# In[664]:


dataFrame_02['Brightness Distribution'] = (dataFrame_02['Amplitude(W/Hz)polar']/(2*1.38e-23*0.5))
print(dataFrame_02)


# In[665]:


plt.figure(figsize=[20,5])
plt.plot(dataFrame_02['Prop. Velocity(km/s)'], dataFrame_02['Brightness Distribution2'])
plt.xlabel('Radial Velocity (km/s)')
plt.ylabel('Brightness Temperature (K)')
plt.title('Brightness Temperature vs. Radial Velocity')
plt.xticks(np.arange(-35, 20, 5.0))
plt.grid()


# In[666]:


plt.figure(figsize=[20,5])
plt.plot(dataFrame_02['Prop. Velocity(km/s)'], dataFrame_02['Brightness Distribution2'])
plt.xlabel('Raddial Velocity (km/s)')
plt.ylabel('Brightness Temperature (K)')


# In[667]:


dataFrame_02['Column Density'] = 1.8224e18


# In[ ]:





# In[668]:


dataFrame_02 = dataFrame_02[~(dataFrame_02['Prop. Velocity(km/s)'] >= 100)]


# In[670]:


# print(dataFrame_02)
dataFrame_02.iat[338, 10]
dataFrame_02trim = dataFrame_02[~(dataFrame_02['Prop. Velocity(km/s)'] >= 10)]
dataFrame_02doubletrim = dataFrame_02trim[~(dataFrame_02trim['Prop. Velocity(km/s)'] <= -30)]
print(dataFrame_02doubletrim)
print(dataFrame_02doubletrim.iat[8, 2])


# In[674]:


absinthe = (3*10**5)*(1-((xsmooth7)/(1420.4e6)))
# print(absinthe)
n=43
sum1=0
for i in range(0,43):
    sum1=sum1 + dataFrame_02doubletrim.iat[i, 10]*(dataFrame_02doubletrim.iat[8, 3]-dataFrame_02doubletrim.iat[9, 3])
Column = sum1*1.8224e18
print(Column)


# In[ ]:





# In[ ]:




