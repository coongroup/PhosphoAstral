# -*- coding: utf-8 -*-
"""
Created on Fri May 10 12:37:58 2024

@author: nlancaster
"""

raw_file = "P:/Projects/NML_TPC_2023_ThermoFisher_Phosphoproteomics_q-OT-Astral/20230412_HEK293T_phospho_rep/20230412_HEK239T_April_phosphobatch_250ng_15min_2Th_500AGC_rep01.raw"

def RT_MS1_MS2(rawfile,ms2_mass):#add rawfile object
    ms1_RT = []
    ms2_RT = []
    total_scans = rawfile.GetNumSpectra()
    for x in range(1,total_scans):
        if rawfile.GetMSOrderForScanNum(x) == 1:
            ms1_RT.append(rawfile.RTFromScanNum(x))
            print('RT Extraction: '+ str(x) + ' of ' + str(total_scans))
            
        elif rawfile.GetMSOrderForScanNum(x) == 2:
            if rawfile.GetPrecursorMassForScanNum(x,2)==ms2_mass:
                ms2_RT.append(rawfile.RTFromScanNum(x))
                print('RT Extraction: '+ str(x) + ' of ' + str(total_scans))

    return ms1_RT, ms2_RT 

def _p_cyc_times(rawfile,ms2_mass):
    ms1_RT,ms2_RT = RT_MS1_MS2(rawfile, ms2_mass)
    ms1 = []
    ms2 = []
    length = len(ms1_RT)
    for i,x in enumerate(ms1_RT):
        if i+1 == length: break
        ms1.append((ms1_RT[i+1]-ms1_RT[i])*60*1000)
        print('MS1 Cycle Time: ' +str(i) + ' of ' + str(length))
    length = len(ms2_RT)
    for i,x in enumerate(ms2_RT):
        if i+1 == length: break
        ms2.append((ms2_RT[i+1]-ms2_RT[i])*60*1000)
        print('MS2 Cycle Time: ' +str(i) + ' of ' + str(length))
    return (ms1_RT[1:],ms1),(ms2_RT[1:],ms2)
def cyc_times(rawfile):
    from pymsfilereader import MSFileReader
    raw_obj = MSFileReader(rawfile)
    ms2_mass = 0
    while ms2_mass == 0:
        scan_num = 2
        if raw_obj.GetMSOrderForScanNum(scan_num) == 2:
            ms2_mass = raw_obj.GetPrecursorMassForScanNum(scan_num,2)
        else: scan_num = scan_num+1
    ms1,ms2 = _p_cyc_times(raw_obj,ms2_mass)
    return ms1,ms2
    

cyc_time_dict = {}


cyc_time_dict['Result'] = cyc_times(raw_file)



def make_cyc_time_plot(dictionary,key,save_file= False):
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize = (3.333333,2.5))

    ax= plt.subplot()
    ax.scatter(dictionary[key][0][0],dictionary[key][0][1],marker = '.',label = 'MS1',color = 'black',s= 6)
    ax.scatter(dictionary[key][1][0],dictionary[key][1][1],marker = '.',label ='MS2',color = 'royalblue',s=6)
    ax.set_ylim(0,2200)
    ax.spines[['top','right']].set_visible(False)
    ax.set_xlabel('Retention Time (min)',fontweight = 'bold',fontsize = 12)
    ax.set_ylabel('Cycle Time (ms)',fontweight = 'bold',fontsize = 12)
    ax.legend(fontsize = 10,loc = 'upper right',markerscale = 2)
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['font.family'] = 'Arial'
    fig.tight_layout()
    
    if save_file:
        fig.savefig(save_file,dpi=800)
    else: 
        fig.show()
        

make_cyc_time_plot(cyc_time_dict,'Result')