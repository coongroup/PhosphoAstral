# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 17:32:52 2023

@author: nlancaster
"""

#Supplemental Figure 1
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

phospho_report = pd.read_csv('P:/Projects/NML_TPC_2023_ThermoFisher_Phosphoproteomics_q-OT-Astral/SpectronautSearchOutputs/20230523_0p75_Pulsar_cutoff_searches/20230412_HEK239T_April_phosphobatch_250ng_0p75PulasrCutoff_30minDIA_DecemberLibS_-PTMSiteReport.tsv',sep='\t')
phospho_report = phospho_report[phospho_report['PTM.ModificationTitle']=='Phospho (STY)']
phospho_report = phospho_report[phospho_report['PTM.SiteProbability']>=0.75]
file_names = phospho_report['R.FileName'].unique()
HEK_results = {}

def phospho_list(data, file):
    from copy import deepcopy
    file_data = data[data['R.FileName']==file]
    
    
    df = deepcopy(file_data)
    df = df.sort_values('PTM.CollapseKey')
    df = df.drop_duplicates('PTM.Group')
    df = df.drop_duplicates('PTM.CollapseKey')
    total_unique = len(df['PTM.Group'])
    
    return total_unique
    



for x in file_names:
    HEK_results[x] = phospho_list(phospho_report,x)


phospho_7min_2Th = [HEK_results['20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep01'],HEK_results['20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep02'],HEK_results['20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep03']]
phospho_15min_2Th = [HEK_results['20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep01'],HEK_results['20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep02'],HEK_results['20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep03']] 
phospho_30min_2Th = [HEK_results['20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep01'],HEK_results['20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep02'],HEK_results['20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep03']] 


labels = [7,15,30]
raw_data = [phospho_7min_2Th, phospho_15min_2Th,phospho_30min_2Th]
def normalize_data(data_list):
    import numpy as np
    avgs = []
    for x in data_list:
        avgs.append(np.mean(x))
    max_val = max(avgs)
    output = []
    for x in data_list:
        output.append([y/max_val*100 for y in x ])
    return output
rel_data = raw_data
data = []
def avg_min_max(arg):
    import numpy as np
    output = [np.mean(arg),np.mean(arg)-min(arg),max(arg)-np.mean(arg)]
    return output
for x in rel_data:
    data.append(avg_min_max(x))    



import numpy as np


width = 0.5
colors = ['silver','cornflowerblue','#5084C4']
linewidth = 0.6
capsize = 5
spacing = width * 1
height_arr = np.array([data[0][0],data[1][0],data[2][0]])
err_array = np.array([[data[0][1],data[1][1],data[2][1]],[data[0][2],data[1][2],data[2][2]]])


chimerys_raw_data = [[31799,31284,31582],[41739,41771,42306],[50996,51676,51761]]

data = []
for x in chimerys_raw_data:
    data.append(avg_min_max(x))    





import numpy as np


width = 0.5
colors = ['silver','cornflowerblue','#5084C4']
linewidth = 0.6
capsize = 5
spacing = width * 1
chimerys_height_arr = np.array([data[0][0],data[1][0],data[2][0]])
chimerys_err_arr = np.array([[data[0][1],data[1][1],data[2][1]],[data[0][2],data[1][2],data[2][2]]])




width_multiplier = 0.7*0.7
height_multiplier = 1*0.7
fig = plt.figure()
ax = plt.subplot() 
x = np.arange(len(labels))
ax.spines[['top','right']].set_visible(False)
ax.set_yticks([0,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,55000])
ax.set_ylim(0,)
x = np.arange(3)
width = 0.35
ax.bar(x - 0.2,height_arr,width,color = '#5084C4',yerr = err_array,capsize = 5,label = 'Spectronaut 17')
ax.bar(x + 0.2,chimerys_height_arr,width,color = 'silver',yerr = chimerys_err_arr,capsize = 5, label = 'PD 3.1 (CHIMERYS)')
ax.set_xticks([0,1,2])
ax.set_xticklabels(['7','15','30'])
# ax.errorbar(labels,height_arr,yerr = err_array,capsize = 4,linestyle = '')
ax.legend()
ax.set_ylabel('Localized Phosphosites',fontweight = 'bold')
ax.set_xlabel('Gradient Length (min)',fontweight = 'bold')
file_save_path= 'C:\\Users\\nlancaster\OneDrive - UW-Madison\Astral_Phosphoproteomics_Manuscript\PythonFigures\\'
# fig.savefig(file_save_path + '20231205_CHIMERYS_Spectronaut_comparison.svg')
# fig.savefig(file_save_path + '20231205_CHIMERYS_Spectronaut_comparison.png',dpi=1000)
