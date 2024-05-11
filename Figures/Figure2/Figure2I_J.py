# -*- coding: utf-8 -*-
"""
Created on Fri May 10 09:29:08 2024

@author: nlancaster
"""
#% Combined FLR and Quant Eval Plot Data Prep beginning
import pandas as pd
import numpy as np
standards = pd.read_csv('P:/Projects/NML_TPC_2023_ThermoFisher_Phosphoproteomics_q-OT-Astral/SyntheticPhospho/SyntheticPhosphoStds_ReferenceSheet.csv')

def parse_accuracy_results(df,standards,sequence_repeats =1 ,probability_cutoff = 0.75,min_combinations = False,multiplicity = False,missed_cleavages =False):
    standards = standards[standards['In Yeast']==False]
    if min_combinations:
        standards = standards[standards['Possible Combinations']>=min_combinations] #1 combination would have no chance of mislocalization
    standards = standards[standards['Sequence Repeats']==1] # removes phospho positional isomers in dataset (would need to manually validate)
    if multiplicity:
        standards = standards[standards['Multiplicity']==multiplicity]
    if missed_cleavages:
        standards = standards[standards['Missed Cleavages']==missed_cleavages]
    df = df[df['PG.Database']=='nml']
    df = df[df['PEP.StrippedSequence'].isin(standards['Modified Sequences'])]
    df = df.reset_index(drop = True)
    drop_indices = []
    for i,x in enumerate(df['PEP.StrippedSequence']):
        filtered_standards = standards[standards['Modified Sequences']==x]
        filtered_standards = filtered_standards.reset_index(drop = True)
        if len(filtered_standards['Modified Sequences'])>1:
            print("Don't currently have this setup to consider standards that have positional isomers, but this should be filtered out for the sequence repeats filtering above")
        else:
            #these if statements use a continue statement to avoid the drop_indices.append() statement below if the m/z value matches
            if df['FG.Charge'][i]>7: 
                print('Charge state greater than 7 found. Need to raise the charge state range used in the reference sheet')
                break
            if df['FG.Charge'][i]==1:
                if round(df['FG.TheoreticalMz'][i],2) == round(filtered_standards['m/z (+1)'][0],2): continue
            elif df['FG.Charge'][i]==2:
                if round(df['FG.TheoreticalMz'][i],2) == round(filtered_standards['m/z (+2)'][0],2): continue
            elif df['FG.Charge'][i]==3:
                if round(df['FG.TheoreticalMz'][i],2) == round(filtered_standards['m/z (+3)'][0],2): continue
            elif df['FG.Charge'][i]==4:
                if round(df['FG.TheoreticalMz'][i],2) == round(filtered_standards['m/z (+4)'][0],2): continue
            elif df['FG.Charge'][i]==5:
                if round(df['FG.TheoreticalMz'][i],2) == round(filtered_standards['m/z (+5)'][0],2): continue
            elif df['FG.Charge'][i]==6:
                if round(df['FG.TheoreticalMz'][i],2) == round(filtered_standards['m/z (+6)'][0],2): continue
            elif df['FG.Charge'][i]==7:
                if round(df['FG.TheoreticalMz'][i],2) == round(filtered_standards['m/z (+7)'][0],2): continue
            drop_indices.append(i)
    df = df.drop(drop_indices)
    df = df[df['EG.PTMAssayProbability']>=probability_cutoff]
    df = df.reset_index(drop = True)
    
    return df



def phospho_match(df,standards):
    match_status = []
    match_compare = []
    for i,x in enumerate(df['PEP.StrippedSequence']):
        filtered_standards = standards[standards['Modified Sequences']==x]
        filtered_standards = filtered_standards.reset_index(drop=True)
        if len(filtered_standards['Modified Sequences'])>1:
            print('NEED TO FIX THIS FUNCTION TO FACILITATE MANUAL VALIDATION OF THE TRUE POSITIVES')
        else:
            #the standard_site variable contains the modified site list
            standard_site = filtered_standards['Sites (Spectronaut Notation)'][0] #This statement wouldn't work for the positional isomers where there would be multiple phospho configurations for the same amino acid sequence
            # print('Still need to figure out what EG.PTMAssayProbability to use for multiply phosphorylated' )
            data_site = df['EG.ProteinPTMLocations'][i].replace('(','').replace(')','')
            # print(standard_site,data_site)
            match_compare.append(str(standard_site) + ':' + str(data_site)) # this prints out the standard vs the data for manual confirming that the true match and false match values are correct
            if data_site == standard_site:
                match_status.append(True)
            else:
                # print('false')
                match_status.append(False)
    df['Match Status'] = match_status
    df['Match Compare'] = match_compare
    return df


def phospho_error_rate(data,standards,file_name):
    data = data[data['R.FileName']==file_name]
    data = data.reset_index(drop = True)
    loc_probs = [0.75,0.8,0.85,0.9,0.95,0.99,0.999]
    error_rate = []
    num_precursors = []
    for x in loc_probs:
        parsed_results = parse_accuracy_results(data, standards,probability_cutoff = x)
        match_results = phospho_match(parsed_results,standards)
        num_true = len(match_results[match_results['Match Status']==True]['Match Status'])
        num_false = len(match_results[match_results['Match Status']==False]['Match Status'])
        error_rate.append(round(num_false/(num_false+ num_true)*100,4))
        num_precursors.append(num_true)
    return pd.DataFrame({'File':file_name,'Localization Probability':loc_probs,'Error Rate':error_rate,'Correct Precursors':num_precursors})


coon_results = pd.read_csv('P:/Projects/NML_TPC_2023_ThermoFisher_Phosphoproteomics_q-OT-Astral/20240306_PhosphoStdYeast_Astral/20240309_20240308_NML_Yeast_PhosphoStd_dDIA_0.75lib/20240308_NML_Yeast_PhosphoStd_dDIA_0.75lib_Report_NM_FullEG_FG (Normal).tsv',sep = '\t')

coon_files = ['20240306_NML_15min_2Th_250ngYeast_39amol_Stds01',
       '20240306_NML_15min_2Th_250ngYeast_39amol_Stds02',
       '20240306_NML_15min_2Th_250ngYeast_39amol_Stds03',
       '20240306_NML_15min_2Th_250ngYeast_156amol_Stds01',
       '20240306_NML_15min_2Th_250ngYeast_156amol_Stds02',
       '20240306_NML_15min_2Th_250ngYeast_156amol_Stds03',
       '20240306_NML_15min_2Th_250ngYeast_625amol_Stds01',
       '20240306_NML_15min_2Th_250ngYeast_625amol_Stds02',
       '20240306_NML_15min_2Th_250ngYeast_625amol_Stds03',
       '20240306_NML_15min_2Th_250ngYeast_2500amol_Stds01',
       '20240306_NML_15min_2Th_250ngYeast_2500amol_Stds02',
       '20240306_NML_15min_2Th_250ngYeast_2500amol_Stds03',
       '20240306_NML_15min_2Th_250ngYeast_10000amol_Stds01',
       '20240306_NML_15min_2Th_250ngYeast_10000amol_Stds02',
       '20240306_NML_15min_2Th_250ngYeast_10000amol_Stds03']

combined_flr = pd.DataFrame()

for x in coon_files:
    temp_df = phospho_error_rate(coon_results,standards,x)
    combined_flr = pd.concat((combined_flr,temp_df),ignore_index = True)



def avg_error_data(df,files,label):
    import numpy as np
    df = df[df['File'].isin(files)]
    loc_prob = list(df['Localization Probability'].unique())
    avg_flr = []
    min_flr = []
    max_flr = []
    avg_prec = []
    min_prec = []
    max_prec = []
    for x in loc_prob:
        temp_df = df[df['Localization Probability']==x]
        avg_flr.append(np.average(temp_df['Error Rate']))
        min_flr.append(np.average(temp_df['Error Rate'])-min(temp_df['Error Rate']))
        max_flr.append(max(temp_df['Error Rate'])-np.average(temp_df['Error Rate']))
        avg_prec.append(np.average(temp_df['Correct Precursors']))
        min_prec.append(np.average(temp_df['Correct Precursors'])-min(temp_df['Correct Precursors']))
        max_prec.append(max(temp_df['Correct Precursors'])-np.average(temp_df['Correct Precursors']))
    output_df = pd.DataFrame({'Localization Probability':loc_prob,'Experiment':label,'Avg FLR':avg_flr,'Minus FLR':min_flr,'Plus FLR':max_flr,'Avg Precursors':avg_prec,'Minus Precursors':min_prec,'Plus Precursors':max_prec})
    return output_df

flr_chart = avg_error_data(combined_flr,coon_files,'Coon')

#append reported results
loc_probs = [0.75,0.8,0.85,0.9,0.95,0.99,0.999]



#%Construct Combined Figure

import matplotlib.pyplot as plt
plot_coon = flr_chart[flr_chart['Experiment']=='Coon']


fig, ax  = plt.subplots()
ax.spines[['top','right']].set_visible(False)


colors =['#5084C4','#89A2C2','#C1BFBF']

x = np.arange(len(loc_probs))
width = 0.5
linewidth = 1

ax.bar(x,plot_coon['Avg FLR'],linewidth = linewidth, width = width, edgecolor = 'black', label = 'Orbitrap Astral',color = colors[0])
ax.set_xticks(x,loc_probs)
ax.set_xlabel('Localization Probability Cutoff')
ax.set_ylabel('Average Error Rate (%)')

#%% Figure 2J
#Phosphosite Quantitation Evaluation - Astral
import pandas as pd
site_ref = pd.read_csv('P:/Projects/NML_TPC_2023_ThermoFisher_Phosphoproteomics_q-OT-Astral/SyntheticPhospho/SyntheticPhosphoStds_ReferenceSheet_SiteLevel.csv')
site_ref = site_ref[site_ref['Sequence Repeats']==1]
site_ref = site_ref[site_ref['In Yeast']==False]
site_ref = site_ref.reset_index(drop  = True)
def generate_site_id_ref(df):
    ptm_site_id =[]
    for i,x in enumerate(df['Accession Number']):
        ptm_site_id.append(str(x)+'-'+str(df['SiteAA'][i])+str(df['SiteLocation'][i]))
    df['SiteId'] = ptm_site_id
    return df

site_ref = generate_site_id_ref(site_ref)
#read in the PTM Site Report pivot table
ptm_pivot = pd.read_csv('P:/Projects/NML_TPC_2023_ThermoFisher_Phosphoproteomics_q-OT-Astral/20240306_PhosphoStdYeast_Astral/20240309_20240308_NML_Yeast_PhosphoStd_dDIA_0.75lib/20240315_114035_20240308_NML_Yeast_PhosphoStd_dDIA_0.75lib_filtered_0p75cutoff.tsv',sep='\t')

def generate_site_id_data(df):
    ptm_site_id =[]
    for i,x in enumerate(df['PTM.ProteinId']):
        ptm_site_id.append(str(x)+'-'+str(df['PTM.SiteAA'][i])+str(df['PTM.SiteLocation'][i]))
    df['SiteId'] = ptm_site_id
    return df

def filter_pivot_report(site_report,site_ref):
    site_report = site_report[site_report['PG.Database']=='nml']
    site_report = site_report[site_report['PTM.ModificationTitle']=='Phospho (STY)']
    site_report = site_report[site_report['SiteId'].isin(site_ref['SiteId'])]
    site_report = site_report.reset_index(drop = True)
    return site_report
ptm_pivot = generate_site_id_data(ptm_pivot)

filtered_pivot = filter_pivot_report(ptm_pivot,site_ref)


concentration_dict = {
    '39amol':0.0390625,
    '156amol':0.15625,
    '625amol':0.625,
    '2500amol':2.5,
    '10000amol':10}


quantity_dict = {
    '39amol':['[1] 20240306_NML_15min_2Th_250ngYeast_39amol_Stds01.htrms.PTM.Quantity',
    '[2] 20240306_NML_15min_2Th_250ngYeast_39amol_Stds02.htrms.PTM.Quantity',
    '[3] 20240306_NML_15min_2Th_250ngYeast_39amol_Stds03.htrms.PTM.Quantity'],
    '156amol': ['[4] 20240306_NML_15min_2Th_250ngYeast_156amol_Stds01.htrms.PTM.Quantity',
    '[5] 20240306_NML_15min_2Th_250ngYeast_156amol_Stds02.htrms.PTM.Quantity',
    '[6] 20240306_NML_15min_2Th_250ngYeast_156amol_Stds03.htrms.PTM.Quantity'],
    '625amol':['[7] 20240306_NML_15min_2Th_250ngYeast_625amol_Stds01.htrms.PTM.Quantity',
    '[8] 20240306_NML_15min_2Th_250ngYeast_625amol_Stds02.htrms.PTM.Quantity',
    '[9] 20240306_NML_15min_2Th_250ngYeast_625amol_Stds03.htrms.PTM.Quantity'],
    '2500amol':['[10] 20240306_NML_15min_2Th_250ngYeast_2500amol_Stds01.htrms.PTM.Quantity',
    '[11] 20240306_NML_15min_2Th_250ngYeast_2500amol_Stds02.htrms.PTM.Quantity',
    '[12] 20240306_NML_15min_2Th_250ngYeast_2500amol_Stds03.htrms.PTM.Quantity'],
    '10000amol':['[13] 20240306_NML_15min_2Th_250ngYeast_10000amol_Stds01.htrms.PTM.Quantity',
    '[14] 20240306_NML_15min_2Th_250ngYeast_10000amol_Stds02.htrms.PTM.Quantity',
    '[15] 20240306_NML_15min_2Th_250ngYeast_10000amol_Stds03.htrms.PTM.Quantity']
    }


import numpy as np
concentration_list =[0.0390625,0.0390625,0.0390625,0.15625,0.15625,0.15625,0.625,0.625,0.625,2.5,2.5,2.5,10,10,10]
concentration_list = [np.log2(x) for x in concentration_list]

def assemble_curve_df(ptm_pivot,concentration_list,quantity_dict):
    import pandas as pd
    import numpy as np
    curve_dict = {'Column Load (fmol)':concentration_list}
    for x in ptm_pivot['SiteId']:
        quantity_list = []
        site_pivot = ptm_pivot[ptm_pivot['SiteId']==x]
        site_pivot = site_pivot.reset_index(drop = True)
        for y in quantity_dict.keys():
            for z in quantity_dict[y]:
                if site_pivot[z][0] == 'Filtered':
                    quantity_list.append(np.nan)
                else:
                    quantity_list.append(float(site_pivot[z][0]))
        curve_dict[x] = [np.log2(x) for x in quantity_list]
    curve_df = pd.DataFrame(curve_dict)
    return curve_df

calibration_curve = assemble_curve_df(filtered_pivot,concentration_list,quantity_dict)

from scipy.stats import pearsonr
r2= []
for x in filtered_pivot['SiteId']:
    temp_df = calibration_curve[['Column Load (fmol)',x]].dropna(subset = [x])
    if len(temp_df['Column Load (fmol)'].unique())>=3:
        r2.append(pearsonr(temp_df['Column Load (fmol)'],temp_df[x])[0]**2)
print(len(r2))

import matplotlib.pyplot as plt
fig,ax3 = plt.subplots()

#R2 plot for Std Quant

bins = [0,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100]
bins = [x/100 for x  in bins ]
labels = ['<0.8','0.81',
 '0.82',
 '0.83',
 '0.84',
 '0.85',
 '0.86',
 '0.87',
 '0.88',
 '0.89',
 '0.9',
 '0.91',
 '0.92',
 '0.93',
 '0.94',
 '0.95',
 '0.96',
 '0.97',
 '0.98',
 '0.99',
 '1.0']
ax3.spines[['top','right']].set_visible(False)
histogram = np.histogram(r2, bins = bins)
ax3.set_xlabel('R\u00B2')
ax3.set_ylabel('Frequency')
ax3.bar(labels,histogram[0],edgecolor = 'black',width =1 ,linewidth = 0.8)
ax3.set_xticks(['<0.8','0.81',
 '0.82',
 '0.83',
 '0.84',
 '0.85',
 '0.86',
 '0.87',
 '0.88',
 '0.89',
 '0.9',
 '0.91',
 '0.92',
 '0.93',
 '0.94',
 '0.95',
 '0.96',
 '0.97',
 '0.98',
 '0.99',
 '1.0'],['<0.8','',
     '',
     '',
     '',
     '0.85',
     '',
     '',
     '',
     '',
     '0.9',
     '',
     '',
     '',
     '',
     '0.95',
     '',
     '',
     '',
     '',
     '1.0'])
