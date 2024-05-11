# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:39:43 2024

@author: nlancaster
"""

agc100_file1 = 'P:/Projects/NML_TPC_2023_ThermoFisher_Phosphoproteomics_q-OT-Astral/20240306_HeLaMethodDev_Astral/20240306_NML_EGFHeLa_15min_2Th_1e4_AGC_01.raw'
agc500_file1 = 'P:/Projects/NML_TPC_2023_ThermoFisher_Phosphoproteomics_q-OT-Astral/20240306_HeLaMethodDev_Astral/20240306_NML_EGFHeLa_15min_2Th_5e4_AGC_01.raw'

from pymsfilereader import MSFileReader
agc100_raw = MSFileReader(agc100_file1)
agc500_raw = MSFileReader(agc500_file1)

def IonCount_MS2InjTime(raw_file):
    total_scans = raw_file.GetNumSpectra()
    charge_count = []
    inject_time = []
    for x in range(1,total_scans+1):
        if raw_file.GetMSOrderForScanNum(x)==2:
            charge_count.append(raw_file.GetScanHeaderInfoForScanNum(x)['TIC'])
            inject_time.append(raw_file.GetTrailerExtraForScanNum(x)['Ion Injection Time (ms)'])
    return charge_count, inject_time

agc100_TIC, agc100_InjTime = IonCount_MS2InjTime(agc100_raw)
agc500_TIC, agc500_InjTime = IonCount_MS2InjTime(agc500_raw)
#%%Sup Fig 2F
import matplotlib.pyplot as plt
import numpy as np
fig,ax = plt.subplots(figsize = (7,3))
ax.spines[['top','right']].set_visible(False)

labels = ['1e4','5e4']
colors =['#C1BFBF','#5084C4']
bins = list(np.logspace(0,np.log10(851914000),num = 100))
ax.hist(agc100_TIC,bins = bins, histtype = 'step',label = '1e4',color = colors[0])
ax.hist(agc500_TIC,bins = bins, histtype = 'step',label = '5e4',color = colors[1])
ax.legend()
ax.set_xscale('log')
ax.set_xlabel('MS2 Total Ion Current')
ax.set_ylabel('Frequency')

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
fig.tight_layout()
#%%Sup Fig 2G
fig,ax = plt.subplots(figsize = (7,3))
ax.spines[['top','right']].set_visible(False)

labels = ['1e4','5e4']
colors =['#C1BFBF','#5084C4']
bins = list(np.linspace(0,4.1,num= 100))
ax.hist(agc500_InjTime,bins = bins, histtype = 'step',label = '5e4',color = colors[1])
ax.hist(agc100_InjTime,bins = bins, histtype = 'step',label = '1e4',color = colors[0])



ax.legend(reverse= True)
ax.set_xlabel('MS2 Injection Time (ms)')
ax.set_ylabel('Frequency')
ax.set_xlim(0,)
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
fig.tight_layout()
