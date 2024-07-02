# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 17:32:52 2023

@author: nlancaster
"""

#Supplemental Figure 2A
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
ax.set_yticks([0,10000,20000,30000,40000,50000])
ax.set_ylim(0,55000)
x = np.arange(3)
width = 0.35
ax.bar(x - 0.2,height_arr,width,color = '#5084C4',yerr = err_array,capsize = 5,label = 'Spectronaut 17')
ax.bar(x + 0.2,chimerys_height_arr,width,color = 'silver',yerr = chimerys_err_arr,capsize = 5, label = 'PD 3.1 (CHIMERYS)')
ax.set_xticks([0,1,2])
ax.set_xticklabels(['7','15','30'])
# ax.errorbar(labels,height_arr,yerr = err_array,capsize = 4,linestyle = '')


#spectronaut results dots
raw_data = rel_data
spacing = [[-1*width/3-0.2,0-0.2,1*width/3-0.2],[-1*width/3+0.8,0.8,1*width/3+0.8],[-1*width/3+1.8,1.8,1*width/3+1.8]]
for i in range(len(x)):
    # distribute scatter randomly across whole width of bar
    ax.scatter(spacing[i], raw_data[i], edgecolor='black', facecolors = 'none', s = 15 )




#chimerys results dots
raw_data = chimerys_raw_data
spacing = [[-1*width/3+0.2,0+0.2,1*width/3+0.2],[-1*width/3+1.2,1+0.2,1*width/3+1.2],[-1*width/3+2.2,2+0.2,1*width/3+2.2]]
for i in range(len(x)):
    # distribute scatter randomly across whole width of bar
    ax.scatter(spacing[i], raw_data[i], edgecolor='black', facecolors = 'none', s = 15 )





ax.legend()
ax.set_ylabel('Localized Phosphosites',fontsize = 12)
ax.set_xlabel('Gradient Length (min)',fontsize = 12)
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
fig.tight_layout()
file_save_path= 'C:\\Users\\nlancaster\OneDrive - UW-Madison\Astral_Phosphoproteomics_Manuscript\PythonFigures\\'
file_save_path= 'C:\\Users\\nlancaster\OneDrive - UW-Madison\Astral_Phosphoproteomics_Manuscript\\NatCommsSubmission\Review_Round2\DraftResponseDocument\\UpdatedManuscriptPythonFigures'
fig.savefig(file_save_path + '\\20240702_SFig2A_CHIMERYS_Spectronaut.svg')
#%% Supplementary Figure 2D + E Data Prep

#Method Dev Results
import pandas as pd
hela_method_dev = pd.read_csv('P:/Projects/NML_TPC_2023_ThermoFisher_Phosphoproteomics_q-OT-Astral/20240306_HeLaMethodDev_Astral/20240313_20240313_HeLa_MethodDev_dDIA/20240313_HeLa_MethodDev_dDIA_Report_PTMSiteReport (Normal).tsv',sep='\t')
files = hela_method_dev['R.FileName'].unique()
def phospho_list(data, file= False,loc_prob = 0.75):
    from copy import deepcopy
    if file:
        file_data = data[data['R.FileName']==file]
    else: file_data = data
    test_df = deepcopy(file_data)
    test_df = test_df[test_df['PTM.ModificationTitle']=='Phospho (STY)']
    test_df = test_df[test_df['PTM.SiteProbability']>=loc_prob]
    test_df = test_df.sort_values('PTM.CollapseKey')
    test_df = test_df.drop_duplicates(subset = 'PTM.Group')
    test_df = test_df.drop_duplicates(subset= 'PTM.CollapseKey')
    test_df = test_df.dropna(subset = ['PTM.Group'])
    return len(test_df['PTM.Group'])

results = []
for x in files:
    results.append(phospho_list(hela_method_dev,x))

results_df = pd.DataFrame({'File':files,'Count':results})
#%%Supplementary Figure 2D
# 380-980 m/z vs 480-1080 m/z 

import numpy as np


def avg_min_max(df,files):
    import numpy as np
    temp_df = df[df['File'].isin(files)]
    return [np.average(temp_df['Count']),max(temp_df['Count'])- np.average(temp_df['Count']),np.average(temp_df['Count'])-min(temp_df['Count'])]


def raw_data_for_plots(df,files):
    temp_df = df[df['File'].isin(files)]
    return temp_df['Count']

mz480_files = ['20240306_NML_EGFHeLa_15min_2Th_480to980mz_01','20240306_NML_EGFHeLa_15min_2Th_480to980mz_02']
mz380_files = ['20240306_NML_EGFHeLa_15min_2Th_5e4_AGC_01','20240306_NML_EGFHeLa_15min_2Th_5e4_AGC_02']
mz380 = avg_min_max(results_df,mz380_files)
mz480 = avg_min_max(results_df,mz480_files)
plot_mz = np.stack([mz380,mz480])

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize = (2.5,3) )
ax.spines[['top','right']].set_visible(False)
width = 0.8
linewidth = 1
colors =['#C1BFBF','#5084C4']
x = np.arange(2)
ax.bar(x,plot_mz[:,0], width = width, linewidth = linewidth, edgecolor = 'black',yerr = (plot_mz[:,1],plot_mz[:,2]),capsize = 10,color = colors)
ax.set_xticks(x,['380-980','480-1080'])

ax.set_ylabel('Localized Phosphosites')
ax.set_xlabel('DIA m/z range')

#add dots on plot for each injection
raw_data = [raw_data_for_plots(results_df, mz380_files),raw_data_for_plots(results_df, mz480_files)]

spacing = [[-1*width/3,1*width/3],[[-1*width/3+1,1+1*width/3]]]
for i in range(len(x)):
    # distribute scatter randomly across whole width of bar
    ax.scatter(spacing[i], raw_data[i], edgecolor='black', facecolors = 'none', s = 15 )



plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
fig.tight_layout()
file_save_path= 'C:\\Users\\nlancaster\OneDrive - UW-Madison\Astral_Phosphoproteomics_Manuscript\\NatCommsSubmission\Review_Round2\DraftResponseDocument\\UpdatedManuscriptPythonFigures'
fig.savefig(file_save_path + '\\20240702_SFig2D_mz_range.svg')

#%% Supplemental Figure 2E
#AGC Astral comparison
import numpy as np


def avg_min_max(df,files):
    import numpy as np
    temp_df = df[df['File'].isin(files)]
    return [np.average(temp_df['Count']),max(temp_df['Count'])- np.average(temp_df['Count']),np.average(temp_df['Count'])-min(temp_df['Count'])]

agc100_files = ['20240306_NML_EGFHeLa_15min_2Th_1e4_AGC_01','20240306_NML_EGFHeLa_15min_2Th_1e4_AGC_02']
agc500_files = ['20240306_NML_EGFHeLa_15min_2Th_5e4_AGC_01','20240306_NML_EGFHeLa_15min_2Th_5e4_AGC_02']
agc100  = avg_min_max(results_df,agc100_files)
agc500 = avg_min_max(results_df,agc500_files)
plot_agc = np.stack([agc100,agc500])

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize = (2.5,3) )
ax.spines[['top','right']].set_visible(False)
width = 0.8
linewidth = 1
colors =['#C1BFBF','#5084C4']
x = np.arange(2)
ax.bar(x,plot_agc[:,0], width = width, linewidth = linewidth, edgecolor = 'black',yerr = (plot_agc[:,1],plot_agc[:,2]),capsize = 10,color = colors)
ax.set_xticks(x,['1e4','5e4'])

ax.set_ylabel('Localized Phosphosites')
ax.set_xlabel('AGC Target (# of charges)')

#add dots on plot for each injection
raw_data = [raw_data_for_plots(results_df, agc100_files),raw_data_for_plots(results_df, agc500_files)]

spacing = [[-1*width/3,1*width/3],[[-1*width/3+1,1+1*width/3]]]
for i in range(len(x)):
    # distribute scatter randomly across whole width of bar
    ax.scatter(spacing[i], raw_data[i], edgecolor='black', facecolors = 'none', s = 15 )




plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
fig.tight_layout()
file_save_path= 'C:\\Users\\nlancaster\OneDrive - UW-Madison\Astral_Phosphoproteomics_Manuscript\\NatCommsSubmission\Review_Round2\DraftResponseDocument\\UpdatedManuscriptPythonFigures'
fig.savefig(file_save_path + '\\20240702_SFig2E_AGC.svg')
