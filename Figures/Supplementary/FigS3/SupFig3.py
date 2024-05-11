# -*- coding: utf-8 -*-
"""
Created on Fri May 10 12:01:37 2024

@author: nlancaster
"""
#Sup Fig 3A-C
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
    #Now that I've filtered out the strippled sequence, I need to filter out for any peptides that are not the correct m/z (could be from phospho loss during sample prep, oxidation during sample prep, etc.)
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
olsen_results = pd.read_csv('P:/Projects/NML_TPC_2023_ThermoFisher_Phosphoproteomics_q-OT-Astral/20240319_Olsen_PhosphoDIA_StdReSearch/20240319_20240319_Olsen_phospoDIA_PhosphoStd_dDIA_ReSearch/20240328_163050_20240319_Olsen_phospoDIA_PhosphoStd_dDIA_ReSearch_FullEG_FG_ReExport.tsv',sep = '\t')
#%
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
olsen_files = ['20181116_QE3_nLC3_AH_SA_PhosTest_DIA_1xYeast-JPTstySig123_1',
       '20181116_QE3_nLC3_AH_SA_PhosTest_DIA_1xYeast-JPTstySig123_2',
       '20181116_QE3_nLC3_AH_SA_PhosTest_DIA_1xYeast-JPTstySig123_3',
       '20181116_QE3_nLC3_AH_SA_PhosTest_DIA_10xYeast-JPTstySig123_1',
       '20181116_QE3_nLC3_AH_SA_PhosTest_DIA_10xYeast-JPTstySig123_2',
       '20181116_QE3_nLC3_AH_SA_PhosTest_DIA_10xYeast-JPTstySig123_3',
       '20181116_QE3_nLC3_AH_SA_PhosTest_DIA_100xYeast-JPTstySig123_1',
       '20181116_QE3_nLC3_AH_SA_PhosTest_DIA_100xYeast-JPTstySig123_2',
       '20181116_QE3_nLC3_AH_SA_PhosTest_DIA_100xYeast-JPTstySig123_3',
       '20181116_QE3_nLC3_AH_SA_PhosTest_DIA_1000xYeast-JPTstySig123_1',
       '20181116_QE3_nLC3_AH_SA_PhosTest_DIA_1000xYeast-JPTstySig123_2',
       '20181116_QE3_nLC3_AH_SA_PhosTest_DIA_1000xYeast-JPTstySig123_3']
combined_flr = pd.DataFrame()

for x in coon_files:
    temp_df = phospho_error_rate(coon_results,standards,x)
    combined_flr = pd.concat((combined_flr,temp_df),ignore_index = True)
for x in olsen_files:
    temp_df = phospho_error_rate(olsen_results,standards,x)
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

flr_chart = pd.concat((avg_error_data(combined_flr,coon_files,'Coon'),avg_error_data(combined_flr,olsen_files,'Olsen')),ignore_index = True)

#append reported results
loc_probs = [0.75,0.8,0.85,0.9,0.95,0.99,0.999]
olsen_reported_flr = pd.DataFrame({'Localization Probability':loc_probs,'Experiment':'Olsen_reported','Avg FLR':[5.8,5.5,4.7,3.9,3.4,1.7,0.9],'Min FLR':np.nan,'Max FLR':np.nan,'Avg Precursors':np.nan,'Min Precursors':np.nan,'Max Precursors':np.nan})
flr_chart = pd.concat((flr_chart,olsen_reported_flr),ignore_index = True)


#% Quantitative Data for Plot

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
#%
def generate_site_id_data(df):
    ptm_site_id =[]
    for i,x in enumerate(df['PTM.ProteinId']):
        ptm_site_id.append(str(x)+'-'+str(df['PTM.SiteAA'][i])+str(df['PTM.SiteLocation'][i]))
    df['SiteId'] = ptm_site_id
    return df

def filter_pivot_report(site_report,site_ref):
    site_report = site_report[site_report['PG.Database']=='nml']
    site_report = site_report[site_report['PTM.ModificationTitle']=='Phospho (STY)']
    # site_report = site_report[site_report['PTM.SiteProbability']>=localization]
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


    
#% Phosphosite Quant Completeness data prep
# Plot of completeness of detection across the five injections

#construct initial list where the numbers represent [in0, in1, in2, in3]
c39amol = [0,0,0,0]
c156amol = [0,0,0,0]
c625amol = [0,0,0,0]
c2500amol = [0,0,0,0]
c10000amol = [0,0,0,0]


for x in calibration_curve.columns[1:]:
    if calibration_curve[x][0:3].isna().sum() == 0:
        c39amol[3] = c39amol[3]+1
    elif calibration_curve[x][0:3].isna().sum() == 1:
        c39amol[2] = c39amol[2]+1
    elif calibration_curve[x][0:3].isna().sum() == 2:
        c39amol[1] = c39amol[1]+1
    elif calibration_curve[x][0:3].isna().sum() == 3:
        c39amol[0] = c39amol[0]+1
    if calibration_curve[x][3:6].isna().sum() == 0:
        c156amol[3] = c156amol[3]+1
    elif calibration_curve[x][3:6].isna().sum() == 1:
        c156amol[2] = c156amol[2]+1
    elif calibration_curve[x][3:6].isna().sum() == 2:
        c156amol[1] = c156amol[1]+1
    elif calibration_curve[x][3:6].isna().sum() == 3:
        c156amol[0] = c156amol[0]+1
    if calibration_curve[x][6:9].isna().sum() == 0:
        c625amol[3] = c625amol[3]+1
    elif calibration_curve[x][6:9].isna().sum() == 1:
        c625amol[2] = c625amol[2]+1
    elif calibration_curve[x][6:9].isna().sum() == 2:
        c625amol[1] = c625amol[1]+1
    elif calibration_curve[x][6:9].isna().sum() == 3:
        c625amol[0] = c625amol[0]+1
    if calibration_curve[x][9:12].isna().sum() == 0:
        c2500amol[3] = c2500amol[3]+1
    elif calibration_curve[x][9:12].isna().sum() == 1:
        c2500amol[2] = c2500amol[2]+1
    elif calibration_curve[x][9:12].isna().sum() == 2:
        c2500amol[1] = c2500amol[1]+1
    elif calibration_curve[x][9:12].isna().sum() == 3:
        c2500amol[0] = c2500amol[0]+1
    if calibration_curve[x][12:15].isna().sum() == 0:
        c10000amol[3] = c10000amol[3]+1
    elif calibration_curve[x][12:15].isna().sum() == 1:
        c10000amol[2] = c10000amol[2]+1
    elif calibration_curve[x][12:15].isna().sum() == 2:
        c10000amol[1] = c10000amol[1]+1
    elif calibration_curve[x][12:15].isna().sum() == 3:
        c10000amol[0] = c10000amol[0]+1
    

completeness_df = pd.DataFrame({'39amol':c39amol,'156amol':c156amol,'625amol':c625amol,'2500amol':c2500amol,'10000amol':c10000amol})




#%Construct Combined Figure
import numpy as np
import matplotlib.pyplot as plt
plot_coon = flr_chart[flr_chart['Experiment']=='Coon']
plot_olsen = flr_chart[flr_chart['Experiment']=='Olsen']
plot_olsen_reported = flr_chart[flr_chart['Experiment']=='Olsen_reported']



fig, ((ax,ax2),(ax3,ax4))  = plt.subplots(2,2,figsize = (7,6))
ax.spines[['top','right']].set_visible(False)


colors =['#5084C4','#89A2C2','#C1BFBF']

x = np.arange(len(loc_probs))
width = 0.2
linewidth = 1
ax.bar(x-width,plot_olsen_reported['Avg FLR'],linewidth = linewidth, width = width, edgecolor = 'black',label = 'Olsen Reported', color = colors[2])
ax.bar(x,plot_olsen['Avg FLR'],linewidth = linewidth, width = width, edgecolor = 'black',label = 'Olsen Re-Analysis ',color = colors[1])
ax.bar(x+width,plot_coon['Avg FLR'],linewidth = linewidth, width = width, edgecolor = 'black', label = 'Orbitrap Astral',color = colors[0])
ax.set_xticks(x,loc_probs)
ax.set_xlabel('Localization Probability Cutoff')
ax.set_ylabel('Average Error Rate (%)')
ax.legend(handlelength = 0.75)




#Correct Precursors vs Loc Prob
width = 0.5
ax2.bar(x,plot_coon['Avg Precursors'],linewidth = linewidth, width = width, edgecolor = 'black', label = 'Coon Astral',color = colors[0])

ax2.set_xticks(x,loc_probs)
ax2.set_xlabel('Localization Probability Cutoff')
ax2.set_ylabel('Average Correct Precursors')
ax2.spines[['top','right']].set_visible(False)






#Phosphosite Quant Completeness

labels = list(completeness_df.columns)
x = np.arange(5)
ax3.spines[['top','right']].set_visible(False)
width = 0.05

colors =['#5084C4','#7698C2','#9BABC1','#C1BFBF']
edge_width = 0.8
bottom1 = completeness_df.loc[3] + completeness_df.loc[2]
bottom0 =bottom1 + completeness_df.loc[1]
ax3.bar(x,completeness_df.loc[0].to_list(),label='0 of 3',bottom = bottom0.to_list(),color = colors[3],edgecolor = 'black',linewidth = edge_width)
ax3.bar(x,completeness_df.loc[1].to_list(),label='1 of 3',bottom =bottom1.to_list(),color = colors[2],edgecolor = 'black',linewidth = edge_width)
ax3.bar(x,completeness_df.loc[2].to_list(),label='2 of 3',bottom = completeness_df.loc[3].to_list(),color = colors[1],edgecolor = 'black',linewidth = edge_width)
ax3.bar(x,completeness_df.loc[3].to_list(),label='3 of 3',color = colors[0],edgecolor = 'black',linewidth = edge_width)


ax3.set_xticks(x,[x.replace('amol','') for x in labels])
ax3.set_xlabel('Column Load (amol)')
ax3.set_ylabel('Phosphorylation Sites')
ax3.legend(loc = (0.98,0.5),handlelength = 0.85)





ax4.axis('off')


plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
fig.tight_layout()

#%% Sup Fig 3D
#Calibration Curve Examples

from scipy.stats import pearsonr
from scipy.stats import linregress

fig,ax = plt.subplots(figsize = (5,2))
ax.spines[['top','right']].set_visible(False)

concentration_list =[0.0390625,0.15625,0.625,2.5,10]
log2_conc = [np.log2(x) for x in concentration_list]

site = 'S109-S3'


ax.set_xlabel('Log2 Column Load (fmol)')
ax.set_ylabel('Log2 PTM.Quantity')
curve_df = pd.DataFrame({'Conc':calibration_curve['Column Load (fmol)'],'Signal':calibration_curve[site]})
curve_df = curve_df.dropna(subset = 'Signal')
lin_regress = linregress(curve_df['Conc'],curve_df['Signal'])
ax.annotate( 'y = ' + str(round(lin_regress[0],3))+ 'x + ' + str(round(lin_regress[1],2))+ '\nR\u00B2 = ' + str(round(lin_regress[2]*lin_regress[2],3)),(-5,14) )

conc = [min(calibration_curve['Column Load (fmol)']),max(calibration_curve['Column Load (fmol)'])]

def linreg_fxn(x,lin_regress):
    y = lin_regress[0]*x + lin_regress[1]
    return y
signal = [linreg_fxn(x,lin_regress) for x in conc]

ax.plot(conc,signal,linestyle ='--',color = 'black',linewidth = 1)
ax.scatter(calibration_curve['Column Load (fmol)'],calibration_curve[site],marker = 's',s = 4)
ax.set_title(site)
ax.set_xticks(log2_conc)


plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
fig.tight_layout()


#%%



fig,ax = plt.subplots(figsize = (5,2))
ax.spines[['top','right']].set_visible(False)

concentration_list =[0.0390625,0.15625,0.625,2.5,10]
log2_conc = [np.log2(x) for x in concentration_list]


site = 'S121-S3'


ax.set_xlabel('Log2 Column Load (fmol)')
ax.set_ylabel('Log2 PTM.Quantity')
curve_df = pd.DataFrame({'Conc':calibration_curve['Column Load (fmol)'],'Signal':calibration_curve[site]})
curve_df = curve_df.dropna(subset = 'Signal')
lin_regress = linregress(curve_df['Conc'],curve_df['Signal'])
ax.annotate( 'y = ' + str(round(lin_regress[0],3))+ 'x + ' + str(round(lin_regress[1],2))+ '\nR\u00B2 = ' + str(round(lin_regress[2]*lin_regress[2],3)),(-5,14) )

conc = [min(calibration_curve['Column Load (fmol)']),max(calibration_curve['Column Load (fmol)'])]

def linreg_fxn(x,lin_regress):
    y = lin_regress[0]*x + lin_regress[1]
    return y
signal = [linreg_fxn(x,lin_regress) for x in conc]

ax.plot(conc,signal,linestyle ='--',color = 'black',linewidth = 1)
ax.scatter(calibration_curve['Column Load (fmol)'],calibration_curve[site],marker = 's',s = 4)
ax.set_title(site)
ax.set_xticks(log2_conc)


plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
fig.tight_layout()
#%%



fig,ax = plt.subplots(figsize = (5,2))
ax.spines[['top','right']].set_visible(False)

concentration_list =[0.0390625,0.15625,0.625,2.5,10]
log2_conc = [np.log2(x) for x in concentration_list]



site = 'S29-Y4'

ax.set_xlabel('Log2 Column Load (fmol)')
ax.set_ylabel('Log2 PTM.Quantity')
curve_df = pd.DataFrame({'Conc':calibration_curve['Column Load (fmol)'],'Signal':calibration_curve[site]})
curve_df = curve_df.dropna(subset = 'Signal')
lin_regress = linregress(curve_df['Conc'],curve_df['Signal'])
ax.annotate( 'y = ' + str(round(lin_regress[0],3))+ 'x + ' + str(round(lin_regress[1],2))+ '\nR\u00B2 = ' + str(round(lin_regress[2]*lin_regress[2],3)),(-5,15.3) )

conc = [min(calibration_curve['Column Load (fmol)']),max(calibration_curve['Column Load (fmol)'])]

def linreg_fxn(x,lin_regress):
    y = lin_regress[0]*x + lin_regress[1]
    return y
signal = [linreg_fxn(x,lin_regress) for x in conc]

ax.plot(conc,signal,linestyle ='--',color = 'black',linewidth = 1)
ax.scatter(calibration_curve['Column Load (fmol)'],calibration_curve[site],marker = 's',s = 4)
ax.set_title(site)
ax.set_xticks(log2_conc)


plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
fig.tight_layout()

#%%Sup Fig 3E
# Intensity distributions for yeast and phosphopeptide standards

import pandas as pd


ptm_pivot = pd.read_csv('P:/Projects/NML_TPC_2023_ThermoFisher_Phosphoproteomics_q-OT-Astral/20240306_PhosphoStdYeast_Astral/20240309_20240308_NML_Yeast_PhosphoStd_dDIA_0.75lib/20240315_114035_20240308_NML_Yeast_PhosphoStd_dDIA_0.75lib_filtered_0p75cutoff.tsv',sep='\t')



def generate_site_id_data(df):
    ptm_site_id =[]
    for i,x in enumerate(df['PTM.ProteinId']):
        ptm_site_id.append(str(x)+'-'+str(df['PTM.SiteAA'][i])+str(df['PTM.SiteLocation'][i]))
    df['SiteId'] = ptm_site_id
    return df

ptm_pivot = generate_site_id_data(ptm_pivot)
ptm_pivot = ptm_pivot[ptm_pivot['PTM.ModificationTitle']=='Phospho (STY)']
ptm_pivot = ptm_pivot.reset_index(drop = True)



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


condition_dict ={
    '39amol':['[1] 20240306_NML_15min_2Th_250ngYeast_39amol_Stds01.htrms.PTM.Quantity',
    '[2] 20240306_NML_15min_2Th_250ngYeast_39amol_Stds02.htrms.PTM.Quantity',
    '[3] 20240306_NML_15min_2Th_250ngYeast_39amol_Stds03.htrms.PTM.Quantity'],
    '156amol':['[4] 20240306_NML_15min_2Th_250ngYeast_156amol_Stds01.htrms.PTM.Quantity',
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
    '[15] 20240306_NML_15min_2Th_250ngYeast_10000amol_Stds03.htrms.PTM.Quantity']}




def filter_pivot_report(site_report,site_ref):
    site_report = site_report[site_report['PG.Database']=='nml']
    site_report = site_report[site_report['PTM.ModificationTitle']=='Phospho (STY)']

    site_report = site_report[site_report['SiteId'].isin(site_ref['SiteId'])]
    site_report = site_report.reset_index(drop = True)
    return site_report

yeast_df = ptm_pivot[ptm_pivot['PG.Organisms']=='Saccharomyces cerevisiae (strain ATCC 204508 / S288c)']

def filter_pivot_report(site_report,site_ref):
    site_report = site_report[site_report['PG.Database']=='nml']
    site_report = site_report[site_report['PTM.ModificationTitle']=='Phospho (STY)']
    site_report = site_report[site_report['SiteId'].isin(site_ref['SiteId'])]
    site_report = site_report.reset_index(drop = True)
    return site_report

filtered_pivot = filter_pivot_report(ptm_pivot,site_ref)


yeast_int ={
    '39amol':[],
    '156amol':[],
    '625amol':[],
    '2500amol':[],
    '10000amol':[]}

std_int ={
    '39amol':[],
    '156amol':[],
    '625amol':[],
    '2500amol':[],
    '10000amol':[]}


import numpy as np
for i,x in enumerate(yeast_df['PTM.CollapseKey']):
    for condition in condition_dict.keys():
        file_list = condition_dict[condition]
        for file in file_list:
            value = yeast_df[file][i]
            if value != 'Filtered':
                yeast_int[condition].append(np.log2(float(yeast_df[file][i])))


for i,x in enumerate(filtered_pivot['PTM.CollapseKey']):
    for condition in condition_dict.keys():
        file_list = condition_dict[condition]
        for file in file_list:
            value = filtered_pivot[file][i]
            if value != 'Filtered':
                std_int[condition].append(np.log2(float(filtered_pivot[file][i])))


import matplotlib.pyplot as plt
import numpy as np

fig, ([ax1,ax2],[ax3,ax4],[ax5,ax6]) = plt.subplots(3,2,figsize = (7,9))

#10000amol
ax_1 = ax1.twinx()
ax1.spines[['top','right']].set_visible(False)
ax_1.spines[['top','left']].set_visible(False)

condition = '10000amol'

ax1.set_title(condition,fontweight = 'bold')


ax1.set_ylabel('Yeast Frequency',fontweight = 'bold')
ax_1.spines['right'].set_color('#5084C4')
ax_1.tick_params(axis = 'y', colors = '#5084C4')


std_plot = ax_1.hist(std_int[condition],bins = 50, histtype = 'step', color = '#5084C4',label = 'Standards')[2]
ax_1.axvline(np.median(std_int[condition]),color = '#5084C4',linewidth = 1, linestyle = '--')
yeast_plot = ax1.hist(yeast_int[condition], bins = 50,histtype = 'step',color = 'black',label = 'Yeast')[2]
ax1.axvline(np.median(yeast_int[condition]),color = 'black',linewidth = 1, linestyle = '--')
legend = yeast_plot + std_plot
leg_labels = ['Yeast','Standards']

#2500amol
ax_2 = ax2.twinx()
ax2.spines[['top','right']].set_visible(False)
ax_2.spines[['top','left']].set_visible(False)

condition = '2500amol'

ax2.set_title(condition,fontweight = 'bold')


ax_2.set_ylabel('Standard Frequency',fontweight = 'bold',color = '#5084C4')
ax_2.spines['right'].set_color('#5084C4')
ax_2.tick_params(axis = 'y', colors = '#5084C4')


std_plot = ax_2.hist(std_int[condition],bins = 50, histtype = 'step', color = '#5084C4',label = 'Standards')[2]
ax_2.axvline(np.median(std_int[condition]),color = '#5084C4',linewidth = 1, linestyle = '--')
yeast_plot = ax2.hist(yeast_int[condition], bins = 50,histtype = 'step',color = 'black',label = 'Yeast')[2]
ax2.axvline(np.median(yeast_int[condition]),color = 'black',linewidth = 1, linestyle = '--')
legend = yeast_plot + std_plot
leg_labels = ['Yeast','Standards']


#625amol
ax_3 = ax3.twinx()
ax3.spines[['top','right']].set_visible(False)
ax_3.spines[['top','left']].set_visible(False)

condition = '625amol'

ax3.set_title(condition,fontweight = 'bold')


ax3.set_ylabel('Yeast Frequency',fontweight = 'bold')
ax_3.spines['right'].set_color('#5084C4')
ax_3.tick_params(axis = 'y', colors = '#5084C4')


std_plot = ax_3.hist(std_int[condition],bins = 50, histtype = 'step', color = '#5084C4',label = 'Standards')[2]
ax_3.axvline(np.median(std_int[condition]),color = '#5084C4',linewidth = 1, linestyle = '--')
yeast_plot = ax3.hist(yeast_int[condition], bins = 50,histtype = 'step',color = 'black',label = 'Yeast')[2]
ax3.axvline(np.median(yeast_int[condition]),color = 'black',linewidth = 1, linestyle = '--')
legend = yeast_plot + std_plot
leg_labels = ['Yeast','Standards']

#156amol
ax_4 = ax4.twinx()
ax4.spines[['top','right']].set_visible(False)
ax_4.spines[['top','left']].set_visible(False)
ax4.set_xlabel('log2 PTM.Quantity',fontweight = 'bold')

condition = '156amol'

ax4.set_title(condition,fontweight = 'bold')


ax_4.set_ylabel('Standard Frequency',fontweight = 'bold',color = '#5084C4')
ax_4.spines['right'].set_color('#5084C4')
ax_4.tick_params(axis = 'y', colors = '#5084C4')


std_plot = ax_4.hist(std_int[condition],bins = 50, histtype = 'step', color = '#5084C4',label = 'Standards')[2]
ax_4.axvline(np.median(std_int[condition]),color = '#5084C4',linewidth = 1, linestyle = '--')
yeast_plot = ax4.hist(yeast_int[condition], bins = 50,histtype = 'step',color = 'black',label = 'Yeast')[2]
ax4.axvline(np.median(yeast_int[condition]),color = 'black',linewidth = 1, linestyle = '--')
legend = yeast_plot + std_plot
leg_labels = ['Yeast','Standards']



#39amol
ax_5 = ax5.twinx()
ax5.spines[['top','right']].set_visible(False)
ax_5.spines[['top','left']].set_visible(False)
ax5.set_xlabel('log2 PTM.Quantity',fontweight = 'bold')

condition = '39amol'

ax5.set_title(condition,fontweight = 'bold')


ax5.set_ylabel('Yeast Frequency',fontweight = 'bold')
ax_5.set_ylabel('Standard Frequency',fontweight = 'bold',color = '#5084C4')
ax_5.spines['right'].set_color('#5084C4')
ax_5.tick_params(axis = 'y', colors = '#5084C4')


std_plot = ax_5.hist(std_int[condition],bins = 50, histtype = 'step', color = '#5084C4',label = 'Standards')[2]
ax_5.axvline(np.median(std_int[condition]),color = '#5084C4',linewidth = 1, linestyle = '--')
yeast_plot = ax5.hist(yeast_int[condition], bins = 50,histtype = 'step',color = 'black',label = 'Yeast')[2]
ax5.axvline(np.median(yeast_int[condition]),color = 'black',linewidth = 1, linestyle = '--')
legend = yeast_plot + std_plot
leg_labels = ['Yeast','Standards']




ax6.axis('off')
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
fig.tight_layout()


