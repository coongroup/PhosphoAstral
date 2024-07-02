# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 17:31:36 2023

@author: nlancaster
"""

#Figure 2A
import pandas as pd
RT = []
mz = []
order = []

RT.append(0)
mz.append(680)
order.append(1)

RT.append(10)
mz.append(381)
order.append(2)


last_ms1_rt = 0
last_mz = 381
last_rt = 10

while last_mz<980:
    if last_rt-last_ms1_rt<600:
        RT.append(last_rt+5)
        mz.append(last_mz + 2)
        order.append(2)
        last_rt = last_rt+5
        last_mz = last_mz + 2
    if last_rt-last_ms1_rt>=600:
        RT.append(last_rt + 5)
        mz.append(680)
        order.append(1)
        last_ms1_rt = last_rt + 5
        last_rt = last_rt + 10
last_mz = 379
while last_rt-last_ms1_rt<512:
        RT.append(last_rt+5)
        mz.append(last_mz + 2)
        order.append(2)
        last_rt = last_rt+5
        last_mz = last_mz + 2
simulated_acq = pd.DataFrame({'RT':RT,'m/z':mz,'Order':order})
def create_box(center_mass, width, rt, total_time):
    box = [[rt,rt,rt + total_time,rt + total_time,rt],[center_mass-width/2,center_mass + width/2,center_mass+width/2,center_mass-width/2,center_mass - width/2]]
    return box
def create_box_MS2(center_mass, width, rt, total_time):
    box = [[rt,rt,rt,rt + total_time,rt + total_time,rt+total_time],[center_mass-width/2,center_mass + width/2,center_mass,center_mass,center_mass+width/2,center_mass - width/2]]
    return box
def generate_scan_boxes(df):
    ms1_results = []
    ms2_results = []
    ms2_masses = []
    
    ms1_total_time=512
    ms2_total_time = 5
    for i,x in enumerate(df['Order']):
        rt = df['RT'][i]
        if x == 1 :
            iso_mass = 680
            iso_wid = 600
            ms1_results.append(create_box(iso_mass,iso_wid,rt,ms1_total_time))
        elif x ==2:
            
            ms2_mass = df['m/z'][i]
            iso_wid = 2
            ms2_results.append(create_box(ms2_mass,iso_wid,rt,ms2_total_time))
            ms2_masses.append(ms2_mass)
    return ms1_results,ms2_results,ms2_masses



ms1, ms2,masses = generate_scan_boxes(simulated_acq)


import matplotlib.pyplot as plt
width_multiplier = 0.7*0.7*2
height_multiplier = 1*0.7*0.85
fig = plt.figure(figsize = (7,3))
ax = plt.subplot()

ax.set_ylim(360,1100)
ax.set_xlim(-50,1800)

ax.spines[['top','right']].set_visible(False)

    
    
ax.set_ylabel('m/z',fontweight = 'bold')
ax.set_xlabel('Time (ms)',fontweight = 'bold')


from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='#5084C4', lw=2),
                Line2D([0], [0], color='#BD222E', lw=2)]
ax.legend(custom_lines ,['MS1','MS2'],loc=(0.88,0.86),frameon=False)


axins = ax.inset_axes([0.09,0.48,0.16*1.5*1.2,0.23999*1.5*0.95])
xmax = (278-250)*1.2+ 250
ymax = (492.5263158- 475)*0.95+475
axins.set_xlim(250,xmax)
axins.set_ylim(475,ymax)
ax.indicate_inset_zoom(axins)
axins.tick_params(axis='both', labelsize=9)


for x in ms1:
    ax.plot(x[0],x[1],color = '#5084C4',linestyle = 'solid',label = 'MS1')
    axins.plot(x[0],x[1],color = '#5084C4',linestyle = 'solid',label = 'MS1')
for x in ms2:
    ax.plot(x[0],x[1],color = '#BD222E',linestyle= 'solid',linewidth = 1,label = 'MS2')
    axins.plot(x[0],x[1],color = '#BD222E',linestyle = 'solid',linewidth = 1,label = 'MS1')

axins.set_xlabel('Time (ms)',fontsize = 9,fontweight = 'bold')
axins.set_ylabel('m/z',fontsize = 9,fontweight = 'bold')



plt.rcParams['svg.fonttype'] = 'none'
fig.tight_layout(pad = 0.2)
file_save_path= 'C:\\Users\\nlancaster\OneDrive - UW-Madison\Astral_Phosphoproteomics_Manuscript\PythonFigures\\'
# fig.savefig(file_save_path + '20230928_DIA_scheme.svg')




#%%Figure 2B
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
    test_df = deepcopy(file_data)
    ptm_list = list(test_df['PTM.Group'].dropna())
    ptm_list = [x.replace('_.0','') for x in ptm_list]
    ptm_list = [x.replace('_.1','') for x in ptm_list]
    ptm_list = [x.replace('_.2','') for x in ptm_list]
    ptm_list = [x.replace('_.3','') for x in ptm_list]
    ptm_list = [x.replace('_.4','') for x in ptm_list]
    ptm_list = [x.replace('_.5','') for x in ptm_list]
    ptm_list = [x.replace('_.6','') for x in ptm_list]
    ptm_list = [x.replace('_.7','') for x in ptm_list]
    ptm_list = [x.replace('_.7','') for x in ptm_list]
    ptm_list = [x.replace('_.8','') for x in ptm_list]
    ptm_list = [x.replace('_.9','') for x in ptm_list]
    ptm_list = [x.replace('_.10','') for x in ptm_list]
    ptm_list = [x.replace('_','') for x in ptm_list]
    
    
    
    return len(list(set(ptm_list)))
for x in file_names:
    HEK_results[x] = phospho_list(phospho_report,x)
#Bin Width Results
phospho_7min_2Th = [HEK_results['20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep01'],HEK_results['20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep02'],HEK_results['20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep03']]
phospho_15min_2Th = [HEK_results['20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep01'],HEK_results['20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep02'],HEK_results['20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep03']] 
phospho_30min_2Th = [HEK_results['20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep01'],HEK_results['20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep02'],HEK_results['20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep03']] 
phospho_30min_4Th = [HEK_results['20230412_HEK239T_April_phosphobatch_250ng_30min_4Th_500AGC_rep01'],HEK_results['20230412_HEK239T_April_phosphobatch_250ng_30min_4Th_500AGC_rep02'],HEK_results['20230412_HEK239T_April_phosphobatch_250ng_30min_4Th_500AGC_rep03']] 


#%%
labels = ['2','4']#,'30min\n\n4Th']
raw_data = [phospho_30min_2Th,phospho_30min_4Th]
data = []
total_data = []
def avg_min_max(arg):
    import numpy as np
    output = [np.mean(arg),np.mean(arg)-min(arg),max(arg)-np.mean(arg)]
    return output
for x in raw_data:
    data.append(avg_min_max(x))    

# for x in total_raw_data:
#     total_data.append(avg_min_max(x))  

width = 0.5
colors = ['silver','cornflowerblue','#5084C4']
linewidth = 0.6
capsize = 5
spacing = width * 1
height_arr = np.array([data[0][0],data[1][0]])
err_array = np.array([[data[0][1],data[1][1]],[data[0][2],data[1][2]]])
fig = plt.figure(figsize = (1.6,2))
ax = plt.subplot()
x = np.arange(len(labels))
ax.bar(x,height_arr,width=width,capsize = capsize,linewidth=linewidth,color = colors[2],align='center',label = 'Localized',yerr=err_array)
# ax.bar(x,tot_height_arr,width=width,capsize = capsize, edgecolor = 'black',linewidth=linewidth,color = colors[1],align='center',label = 'Total',yerr=tot_err_array)


spacing = [[-1*width/3,0,1*width/3],[[-1*width/3+1,1,1+1*width/3]]]
for i in range(len(x)):
    # distribute scatter randomly across whole width of bar
    ax.scatter(spacing[i], raw_data[i], edgecolor='black', facecolors = 'none', s = 10 )

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_yticks([0,10000,20000,30000])
# ax.legend(fancybox = True)

# ax.set_title('HEK293T Phospho Replicates',fontweight = 'bold')


# ax.set_xlabel('Loading Mass (ng)',fontweight= 'bold')
ax.set_ylabel('Localized Phosphosites')
ax.set_xlabel('Isolation Width, Th           ')
# ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.tight_layout()
plt.rcParams['svg.fonttype'] = 'none'
file_save_path= 'C:\\Users\\nlancaster\OneDrive - UW-Madison\Astral_Phosphoproteomics_Manuscript\\NatCommsSubmission\Review_Round2\DraftResponseDocument\\UpdatedManuscriptPythonFigures'
fig.savefig(file_save_path + '\\20240702_Fig2B_isowidths.svg')


#%%Figure 2C
import pandas as pd
path= 'P:/Projects/NML_TPC_2023_ThermoFisher_Phosphoproteomics_q-OT-Astral/SpectronautSearchOutputs/2023715_dil_seires_phospho_30min_0p75library/2023715_dil_seires_phospho_30min_0p75library_-PTMSiteReport.tsv'
dil_series = pd.read_csv(path,sep = '\t')

#%%
file_names = dil_series['R.FileName'].unique()
dil_results = {}

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
    dil_results[x] = phospho_list(dil_series,x)
import matplotlib.pyplot as plt
labels = [50,100,250,500,1000]
width_multiplier = 0.7*0.7*5/4*1.1
height_multiplier = 1*0.7*0.95
height_multiplier = 3.192/4.8
fig = plt.figure(figsize = (2.5,2))
ax = plt.subplot()
dil_data = [dil_results['20221213_Coon_phospho2co_50ng_7to47B_30min_HCD27_2Th_360bins_mz380-1100_MS2_150-2k_3pt5MaxIT'],
            dil_results['20221213_Coon_phospho2co_100ng_7to47B_30min_HCD27_2Th_360bins_mz380-1100_MS2_150-2k_3pt5MaxIT_20221214185646'],
            dil_results['20221216_Coon_phospho1co_250ng_DIA_30min_2Th_380-1100_DefineMZrange150-2000_HCD27_EOWrep'],
            dil_results['20221216_OLEP07_Coon_phospho2co_DIA_500ng_7to47B_30min_HCD27_2Th_0ThOvlp_360bins_mz380-1100_MS2_150-2k_3pt5MaxIT_EOWrep2'],
            dil_results['20221216_OLEP07_Coon_phospho2co_DIA_1ug_7to47B_30min_HCD27_2Th_0ThOvlp_360bins_mz380-1100_MS2_150-2k_3pt5MaxIT_EOWrep1']]

relative_dil_data = [x/max(dil_data)*100 for x in dil_data]

ax.plot(labels,relative_dil_data,marker = 'o',color = '#5084C4')
ax.set_ylim(0,110)
ax.set_xlim(0,)


ax.spines[['top','right']].set_visible(False)
ax.set_xticks([0,250,500,750,1000])
ax.set_yticks([0,20,40,60,80,100])
# ax.set_yticks([0,5000,10000,15000,20000,25000,30000,35000,40000])
# ax.tick_params(axis='x', labelsize=12)
# ax.tick_params(axis = 'y',labelsize = 12)
ax.set_xlabel('Loading Mass (ng)',fontweight = 'bold')
ax.set_ylabel('Relative Depth (%)',fontweight = 'bold')
# fig.tight_layout(pad= 0.2)
plt.rcParams['svg.fonttype'] = 'none'
# ax2.set_ylabel('% Increase in Localized Phosphosites',fontweight = 'bold',color = color)
# fig.savefig('C:\\Users\\nlancaster\OneDrive - UW-Madison\ASMS2023_oral_presentation\FinalFigures\\20230531_dil_series.svg')
fig.tight_layout(pad = 0.2)
file_save_path= 'C:\\Users\\nlancaster\OneDrive - UW-Madison\Astral_Phosphoproteomics_Manuscript\PythonFigures\\'
# fig.savefig(file_save_path + '20231208_dil_series_relative.svg')






#%%Figure 2D
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
rel_data = normalize_data(raw_data)
data = []
def avg_min_max(arg):
    import numpy as np
    output = [np.mean(arg),np.mean(arg)-min(arg),max(arg)-np.mean(arg)]
    return output
for x in rel_data:
    data.append(avg_min_max(x))    



width = 0.5
colors = ['silver','cornflowerblue','#5084C4']
linewidth = 0.6
capsize = 5
spacing = width * 1
height_arr = np.array([data[0][0],data[1][0],data[2][0]])
err_array = np.array([[data[0][1],data[1][1],data[2][1]],[data[0][2],data[1][2],data[2][2]]])
width_multiplier = 0.7*0.7
height_multiplier = 1*0.7
fig = plt.figure(figsize = (2.3333,2))
ax = plt.subplot() 
x = np.arange(len(labels))

ax.plot(labels,height_arr,marker = 'o',color = '#5084C4',label = 'Localized')
ax.errorbar(labels,height_arr,yerr = err_array,capsize = 4,ecolor = 'black',linestyle = '')


ax.set_ylabel('Localized Phosphosites',fontweight = 'bold')
ax.set_xlabel('Gradient Length (min)',fontweight = 'bold')


ax.set_xticks([0,10,20,30])
ax.set_yticks([0,20,40,60,80,100])
ax.spines[['top','right']].set_visible(False)



ax.set_ylabel('Relative Depth (%)',fontweight = 'bold')
fig.tight_layout(pad = 0.2)
ax.set_xlim(0,)
ax.set_ylim(0,110)
file_save_path= 'C:\\Users\\nlancaster\OneDrive - UW-Madison\Astral_Phosphoproteomics_Manuscript\PythonFigures\\'
# fig.savefig(file_save_path + '20231030_gradient_length_relative.svg')

#%%Figure 2E
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

phospho_report = pd.read_csv('P:/Projects/NML_TPC_2023_ThermoFisher_Phosphoproteomics_q-OT-Astral/SpectronautSearchOutputs/20230523_0p75_Pulsar_cutoff_searches/20230412_HEK239T_April_phosphobatch_250ng_0p75PulasrCutoff_30minDIA_DecemberLibS_-PTMSiteReport.tsv',sep='\t')
phospho_report = phospho_report[phospho_report['PTM.ModificationTitle']=='Phospho (STY)']
phospho_report = phospho_report[phospho_report['PTM.SiteProbability']>=0.75]


file_names = {'7min_2Th':[    
    '20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep01',
    '20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep02',
    '20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep03'],
    '15min_2Th':[
    '20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep01',
    '20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep02',
    '20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep03'],
    '30min_2Th':[
    '20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep01',
    '20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep02',
    '20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep03']}
from copy import deepcopy
phospho_7min  = phospho_report[phospho_report['R.FileName'].isin(file_names['7min_2Th'])]
phospho_15min = phospho_report[phospho_report['R.FileName'].isin(file_names['15min_2Th'])]
phospho_30min = phospho_report[phospho_report['R.FileName'].isin(file_names['30min_2Th'])]

    
    
def phosphosite_overlap(data_df):
    file_names = [x for x in data_df['R.FileName'].unique()]
    df = deepcopy(data_df)
    df = df.sort_values('PTM.CollapseKey')
    df = df.drop_duplicates('PTM.Group')
    df = df.drop_duplicates('PTM.CollapseKey')
    total_unique = len(df['PTM.Group'])
    replicate_list = []
    for x in file_names:
        temp_df = deepcopy(data_df[data_df['R.FileName']==x])
        temp_df = temp_df.sort_values('PTM.CollapseKey')
        temp_df = temp_df.drop_duplicates('PTM.Group')
        temp_df = temp_df.drop_duplicates('PTM.CollapseKey')
        replicate_list.append(temp_df)
    total_unique_list = df['PTM.CollapseKey']
    count_list = []
    rep1 = list(replicate_list[0]['PTM.CollapseKey'])
    rep2 = list(replicate_list[1]['PTM.CollapseKey'])
    rep3 = list(replicate_list[2]['PTM.CollapseKey'])
    for i,x in enumerate(total_unique_list):
        temp_val = 0
        if x in rep1:
            temp_val = temp_val +1
        if x in rep2:
            temp_val = temp_val +1
        if x in rep3:
            temp_val = temp_val +1
        count_list.append(temp_val)
        print(str(i+1 )+ ' of ' + str(total_unique)+ ' complete.')
    site_overlap = pd.DataFrame({'Site':total_unique_list,'File Count':count_list})
    return site_overlap

    
ovlp_7 = phosphosite_overlap(phospho_7min)

ovlp_15 = phosphosite_overlap(phospho_15min)

ovlp_30 = phosphosite_overlap(phospho_30min)
val7 = ovlp_7['File Count'].value_counts().sort_index()
val15 = ovlp_15['File Count'].value_counts().sort_index()
val30 = ovlp_30['File Count'].value_counts().sort_index()

#make cluster bar chart
fig, ax = plt.subplots(figsize = (3.5,2))

in3 = [val7[3],val15[3],val30[3]]
in2 = [val7[2],val15[2],val30[2]]
in1 = [val7[1],val15[1],val30[1]]

total_unique = []
for i,x in enumerate(in3): total_unique.append(in3[i]+in2[i]+in1[i])

in1_bottom = []
for i,x in enumerate(in3):
    in1_bottom.append(x + in2[i])
width = 0.83
linewidth = 1
x = ['7','15','30']


ax.bar(x,in1, width = width, edgecolor = 'black',linewidth = linewidth,bottom = in1_bottom, label = '1 of 3',color = '#C0BFBF')
ax.bar(x,in2, width = width, edgecolor = 'black',linewidth = linewidth,bottom = in3, label = '2 of 3',color = '#88A2C2' )
ax.bar(x,in3, width = width, edgecolor = 'black',linewidth = linewidth,label = '3 of 3',color = '#5084C4' )


textloc_in3 = in3
textloc_in2 = []
for i,x in enumerate(in2): textloc_in2.append(x + in3[i])
textloc_in1 = [] 
for i,x in enumerate(in1): textloc_in1.append(in1[i] + in2[i]+in3[i])

text_in3 = []
for i,x in enumerate(in3): text_in3.append(str(round(x/total_unique[i]*100,1))+'%')
text_in2 = []
for i,x in enumerate(in2): text_in2.append(str(round(x/total_unique[i]*100,1))+'%')
text_in1 = []
for i,x in enumerate(in1): text_in1.append(str(round(x/total_unique[i]*100,1))+'%')


text_x = [0,1,2]
for i,x in enumerate(text_in1):
    ax.annotate(x,(text_x[i],textloc_in1[i]-1000),ha= 'center',va = 'top',fontsize = 8)
for i,x in enumerate(text_in3):
    ax.annotate(x,(text_x[i],textloc_in3[i]-1000),ha= 'center',va = 'top',fontsize = 8)
for i,x in enumerate(text_in2):
    ax.annotate(x,(text_x[i],textloc_in2[i]-1000),ha= 'center',va = 'top',fontsize = 8)

ax.spines[['top','right']].set_visible(False)
ax.set_xlabel('Gradient Length (min)',fontweight = 'bold')
ax.set_ylabel('Localized Phosphosites',fontweight = 'bold')
ax.legend(title = 'Data Completeness\nAcross Triplicates',frameon = False,loc = (0.95,0.3))


fig.tight_layout(pad = 0.2)
plt.rcParams['svg.fonttype'] = 'none'

file_save_path= 'C:\\Users\\nlancaster\OneDrive - UW-Madison\Astral_Phosphoproteomics_Manuscript\PythonFigures\\'
# fig.savefig(file_save_path + '20231205_phosphositeOvlpBarChart.svg')


#%%Figure 2F
import pandas as pd
hekrep_pivot = pd.read_csv('P:/Projects/NML_TPC_2023_ThermoFisher_Phosphoproteomics_q-OT-Astral/SpectronautSearchOutputs/20230523_0p75_Pulsar_cutoff_searches/20230412_HEK239T_April_phosphobatch_250ng_0p75PulasrCutoff_30minDIA_DecemberLibSearch_PTMsite_pivot.TSV',sep='\t')

site_probabilities = {'7min_2Th':['[1] 20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep01.raw.PTM.SiteProbability','[2] 20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep02.raw.PTM.SiteProbability','[3] 20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep03.raw.PTM.SiteProbability'],
                  '15min_2Th':['[7] 20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep01.raw.PTM.SiteProbability','[8] 20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep02.raw.PTM.SiteProbability','[9] 20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep03.raw.PTM.SiteProbability'],
                  '30min_2Th':['[13] 20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep01.raw.PTM.SiteProbability','[14] 20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep02.raw.PTM.SiteProbability','[15] 20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep03.raw.PTM.SiteProbability']}

ptm_quantities = {'7min_2Th':['[1] 20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep01.raw.PTM.Quantity','[2] 20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep02.raw.PTM.Quantity','[3] 20230412_HEK239T_April_phosphobatch_250ng_7min_2Th_500AGC_rep03.raw.PTM.Quantity'],
                  '15min_2Th':['[7] 20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep01.raw.PTM.Quantity','[8] 20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep02.raw.PTM.Quantity','[9] 20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep03.raw.PTM.Quantity'],
                  '30min_2Th':['[13] 20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep01.raw.PTM.Quantity','[14] 20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep02.raw.PTM.Quantity','[15] 20230412_HEK239T_April_phosphobatch_250ng_30min_2Th_500AGC_rep03.raw.PTM.Quantity']}#,
quant_dict = {}
import numpy as np
from copy import deepcopy
def calculate_rsd(row):
    mean = np.mean(row)
    std = np.std(row,ddof = len(row)-1)
    rsd = (std/mean) * 100
    return rsd

value_to_drop = 'Filtered'
for x in ptm_quantities.keys():
    temp_df = deepcopy(hekrep_pivot)
    temp_df = temp_df.dropna(subset = ptm_quantities[x])
    temp_df = temp_df.dropna(subset = site_probabilities[x])
    temp_df = temp_df[~temp_df[ptm_quantities[x]].isin([value_to_drop]).any(axis=1)]
    temp_df = temp_df[~temp_df[site_probabilities[x]].isin([value_to_drop]).any(axis=1)]
    
    temp_df[ptm_quantities[x]] = temp_df[ptm_quantities[x]].astype(float)
    temp_df[site_probabilities[x]] = temp_df[site_probabilities[x]].astype(float)
    temp_df = temp_df[temp_df[site_probabilities[x][0]]>=0.75]
    temp_df = temp_df[temp_df[site_probabilities[x][1]]>=0.75]
    temp_df = temp_df[temp_df[site_probabilities[x][2]]>=0.75]
    filtered_df = temp_df
    filtered_df['RSD'] = filtered_df[ptm_quantities[x]].apply(calculate_rsd, axis = 1)
    quant_dict[x] = filtered_df['RSD']
# Violin plots 
import matplotlib.pyplot as plt
labels = [x for x in quant_dict.keys()]
data = [x for x in quant_dict.values()]
width_multiplier = 0.7*0.7
height_multiplier = 1*0.7
fig = plt.figure(figsize = (3.5,2))
ax = plt.subplot() 


import seaborn as sb
import numpy as np

ax = sb.violinplot(data=data,color = '#5084C4',linewidth=1,inner = None,cut =0)

plt.setp(ax.collections, edgecolor="k")

ax.axhline(20,linestyle = '--',color = 'black')
ax.spines[['top','right']].set_visible(False)
ax.set_ylabel('% RSD (n=3)',fontweight = 'bold')
ax.set_ylim(0,)

ax.set_xticklabels(['7','15','30'])
ax.set_xlabel('Gradient Length (min)',fontweight = 'bold')
for i,x in enumerate(data):
    ax.annotate(str(round(np.median(x),1))+'%',xy = (i,255),ha='center',weight = 'bold')
ax.annotate('20%',xy=(2.5,20),va = 'center',fontweight = 'bold')




ax.set_yticks([0,50,100,150,200,250])

plt.rcParams['svg.fonttype'] = 'none'

fig.tight_layout(pad = 0.2)
file_save_path= 'C:\\Users\\nlancaster\OneDrive - UW-Madison\Astral_Phosphoproteomics_Manuscript\PythonFigures\\'
# fig.savefig(file_save_path + '20230928_quant_rsd.svg')
