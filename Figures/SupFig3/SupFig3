# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 17:34:18 2023

@author: nlancaster
"""

#Supplemental Figure 3
#LC Peak Width Evalutions
import pandas as pd
file_path = 'P:/Projects/NML_TPC_2023_ThermoFisher_Phosphoproteomics_q-OT-Astral/SpectronautSearchOutputs/20230412_HEK239T_April_phosphobatch_250ng_0p75PulasrCutoff_30minDIA_DecemberLibSearch_FullEG.tsv'
base_data = pd.read_csv(file_path,sep='\t')
#%%
file_names = {'7min_2Th':[    
    '20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep01',
    '20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep02',
    '20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep03'],
    '7min_4Th':[
    '20230412_HEK239T_April_phosphobatch_250ng_7min_4Th_500AGC_rep01',
    '20230412_HEK239T_April_phosphobatch_250ng_7min_4Th_500AGC_rep02',
    '20230412_HEK239T_April_phosphobatch_250ng_7min_4Th_500AGC_rep03'],
    '15min_2Th':[
    '20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep01',
    '20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep02',
    '20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep03'],
    '15min_4Th':[
    '20230412_HEK239T_April_phosphobatch_250ng_15min_4Th_500AGC_rep01',
    '20230412_HEK239T_April_phosphobatch_250ng_15min_4Th_500AGC_rep02',
    '20230412_HEK239T_April_phosphobatch_250ng_15min_4Th_500AGC_rep03'],
    '30min_2Th':[
    '20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep01',
    '20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep02',
    '20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep03'],
    '30min_4Th':[
    '20230412_HEK239T_April_phosphobatch_250ng_30min_4Th_500AGC_rep01',
    '20230412_HEK239T_April_phosphobatch_250ng_30min_4Th_500AGC_rep02',
    '20230412_HEK239T_April_phosphobatch_250ng_30min_4Th_500AGC_rep03'
    ]}


fwhm_baseline = 1.698643601
fwhm_dict ={}
for x in file_names.keys():
    temp_subset = base_data[base_data['R.FileName'].isin(file_names[x])]
    temp_subset = temp_subset.dropna( subset = ['EG.PTMPositions [Phospho (STY)]'])
    fwhm_dict[x] = [x*60 for x in temp_subset['EG.PeakWidth']]
peak_width_dict = {}
for x in file_names.keys():
    temp_subset = base_data[base_data['R.FileName'].isin(file_names[x])]
    peak_width_dict[x] = [x*60 for x in temp_subset['EG.PeakWidth']]
#Peak widths
ms1_points_dict = {}

for x in file_names.keys():
    temp_subset = base_data[base_data['R.FileName'].isin(file_names[x])]
    ms1_points_dict[x] = [x for x in temp_subset['EG.DatapointsPerPeak (MS1)']]



ms2_points_dict = {}

for x in file_names.keys():
    temp_subset = base_data[base_data['R.FileName'].isin(file_names[x])]
    ms2_points_dict[x] = [x for x in temp_subset['EG.DatapointsPerPeak']]


#%%Chromatographic Peak Widths

import matplotlib.pyplot as plt
labels = [7,15,30]
data = [[x*1.7 for x in fwhm_dict['7min_2Th']],[x*1.7 for x in fwhm_dict['15min_2Th']],[x*1.7 for x in fwhm_dict['30min_2Th']]]
fig, (ax1,ax2) = plt.subplots(2, 1,sharex = True,figsize = (3.3333,2.5))

import seaborn as sb
import numpy as np

medians = [str(round(np.median(x),1)) for x in data]
edge_color = 'black'
body_color= '#6086B6'
violin1 = ax1.violinplot(data,bw_method = 0.3,showextrema = False,showmedians = True)
for violin in violin1['bodies']:
    violin.set_facecolor(body_color)
    violin.set_edgecolor(edge_color)
    violin.set_alpha(1)

violin2 = ax2.violinplot(data,bw_method = 0.3,showextrema = False,showmedians = True)
for violin in violin2['bodies']:
    violin.set_facecolor(body_color)
    violin.set_edgecolor(edge_color)
    violin.set_alpha(1)
    
plt.setp(ax1.collections, edgecolor="k")
plt.setp(ax2.collections, edgecolor="k")
ax1.set_ylim(50,280)

ax2.set_ylim(0,50)



ax1.tick_params(labeltop=False)  # don't put tick labels at the top
ax1.tick_params(labelbottom = False)
ax1.spines[['top','right','bottom']].set_visible(False)

ax2.spines[['top','right']].set_visible(False)

ax2.tick_params(labeltop=False)
ax1.tick_params(bottom = False)
ax2.tick_params(top=False)  # don't put tick labels at the top


    


ax2.set_ylabel('Chromatographic\n Baseline Width (s)',fontweight = 'bold',fontsize = 12, ha = 'center')
ax2.set_xticks([1,2,3])
ax2.set_xticklabels(['7','15','30'])
ax2.set_xlabel('Gradient Length (min)',fontweight = 'bold',fontsize = 12)


ax1.annotate(medians[0]+'s',(1,250),fontsize = 10, fontweight = 'bold', ha = 'center')
ax1.annotate(medians[1]+'s',(2,250),fontsize = 10, fontweight = 'bold', ha = 'center')
ax1.annotate(medians[2]+'s',(3,250),fontsize = 10, fontweight = 'bold', ha = 'center')

d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
plt.subplots_adjust(hspace = 0.05)

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['figure.constrained_layout.use'] = True

file_save_path= 'C:\\Users\\nlancaster\OneDrive - UW-Madison\Astral_Phosphoproteomics_Manuscript\PythonFigures\\'
# fig.savefig(file_save_path + '20231101_fhwm_violin_poster.svg',bbox_inches = 'tight')


#%% Points Across the Peak violin
import matplotlib.pyplot as plt
labels = [7,15,30]
data = [[x*1.7 for x in ms1_points_dict['7min_2Th']],[x*1.7 for x in ms1_points_dict['15min_2Th']],[x*1.7 for x in ms1_points_dict['30min_2Th']]]
width_multiplier = 0.7*0.7
height_multiplier = 1*0.7
fig, (ax1,ax2) = plt.subplots(2, 1,sharex = True,figsize = (3.3333,2.5))

import seaborn as sb
import numpy as np

medians = [str(round(np.median(x),0))[:-2] for x in data]
edge_color = 'black'
body_color= '#6086B6'
violin1 = ax1.violinplot(data,bw_method = 0.3,showextrema = False,showmedians = True)
for violin in violin1['bodies']:
    violin.set_facecolor(body_color)
    violin.set_edgecolor(edge_color)
    violin.set_alpha(1)

violin2 = ax2.violinplot(data,bw_method = 0.3,showextrema = False,showmedians = True)
for violin in violin2['bodies']:
    violin.set_facecolor(body_color)
    violin.set_edgecolor(edge_color)
    violin.set_alpha(1)
    
plt.setp(ax1.collections, edgecolor="k")
plt.setp(ax2.collections, edgecolor="k")
ax1.set_ylim(50,450)

ax2.set_ylim(0,50)



ax1.tick_params(labeltop=False)  # don't put tick labels at the top
ax1.tick_params(labelbottom = False)
ax1.spines[['top','right','bottom']].set_visible(False)

ax2.spines[['top','right']].set_visible(False)

ax2.tick_params(labeltop=False)
ax1.tick_params(bottom = False)
ax2.tick_params(top=False)  # don't put tick labels at the top


    


ax2.set_ylabel('MS1 Points\n Across Baseline',fontweight = 'bold',fontsize = 12, ha = 'center')
ax2.set_xticks([1,2,3])
ax2.set_xticklabels(['7','15','30'])
ax2.set_xlabel('Gradient Length (min)',fontweight = 'bold',fontsize = 12)






ax1.annotate(medians[0],(1,400),fontsize = 10, fontweight = 'bold', ha = 'center')
ax1.annotate(medians[1],(2,400),fontsize = 10, fontweight = 'bold', ha = 'center')
ax1.annotate(medians[2],(3,400),fontsize = 10, fontweight = 'bold', ha = 'center')

d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
plt.subplots_adjust(hspace = 0.05)

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['figure.constrained_layout.use'] = True
file_save_path= 'C:\\Users\\nlancaster\OneDrive - UW-Madison\Astral_Phosphoproteomics_Manuscript\PythonFigures\\'
# fig.savefig(file_save_path + '20231101_ms1_points.svg',dpi=2000,bbox_inches = 'tight')


#%%MS2 points violin plot
import matplotlib.pyplot as plt
labels = [7,15,30]
data = [[x*1.7 for x in ms2_points_dict['7min_2Th']],[x*1.7 for x in ms2_points_dict['15min_2Th']],[x*1.7 for x in ms2_points_dict['30min_2Th']]]

fig, (ax1,ax2) = plt.subplots(2, 1,sharex = True,figsize = (3.3333,2.5))

import seaborn as sb
import numpy as np

medians = [str(round(np.median(x),0))[:-2] for x in data]
edge_color = 'black'
body_color= '#6086B6'
violin1 = ax1.violinplot(data,bw_method = 0.3,showextrema = False,showmedians = True)
for violin in violin1['bodies']:
    violin.set_facecolor(body_color)
    violin.set_edgecolor(edge_color)
    violin.set_alpha(1)

violin2 = ax2.violinplot(data,bw_method = 0.3,showextrema = False,showmedians = True)
for violin in violin2['bodies']:
    violin.set_facecolor(body_color)
    violin.set_edgecolor(edge_color)
    violin.set_alpha(1)
    
plt.setp(ax1.collections, edgecolor="k")
plt.setp(ax2.collections, edgecolor="k")
ax1.set_ylim(20,190)

ax2.set_ylim(0,20)



ax1.tick_params(labeltop=False)  # don't put tick labels at the top
ax1.tick_params(labelbottom = False)
ax1.spines[['top','right','bottom']].set_visible(False)

ax2.spines[['top','right']].set_visible(False)

ax2.tick_params(labeltop=False)
ax1.tick_params(bottom = False)
ax2.tick_params(top=False)  # don't put tick labels at the top

ax2.set_ylabel('MS2 Points\n Across Baseline',fontweight = 'bold',fontsize = 12, ha = 'center')
ax2.set_xticks([1,2,3])
ax2.set_xticklabels(['7','15','30'])
ax2.set_xlabel('Gradient Length (min)',fontweight = 'bold',fontsize = 12)


ax1.annotate(medians[0],(1,170),fontsize = 10, fontweight = 'bold', ha = 'center')
ax1.annotate(medians[1],(2,170),fontsize = 10, fontweight = 'bold', ha = 'center')
ax1.annotate(medians[2],(3,170),fontsize = 10, fontweight = 'bold', ha = 'center')

d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
plt.subplots_adjust(hspace = 0.05)

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['figure.constrained_layout.use'] = True
file_save_path= 'C:\\Users\\nlancaster\OneDrive - UW-Madison\Astral_Phosphoproteomics_Manuscript\PythonFigures\\'
# fig.savefig(file_save_path + '20231101_ms2_points.svg',dpi=2000,bbox_inches = 'tight')
