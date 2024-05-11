# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 17:33:43 2023

@author: nlancaster
"""
import matplot.pyplot as plt
#Supplemental Figure 2
def concave_grad(Vf,Vt,Tf,Tt,Te,curve_type):
    #pass in the curve_type number and k is determined from the dictionary
    if curve_type>=6:
        k_dict = {6:0.75,7:0.5,8:0.25,9:0}
        k = k_dict[curve_type]
        Ve = Vf + (1-k)*(Vt-Vf)*(2**((-10*(Tt-Te))/(Tt-Tf)))/(1-2**(-10)) +((k*(Vt-Vf)*(Te-Tf))/(Tt-Tf))
        return Ve
    else:
        k_dict = {1:0,2:0.25,3:0.5,4:0.75,5:1}
        k = k_dict[curve_type]
        Ve = Vf + (1-k)*(Vt-Vf)*(1-2**((-10*(Te-Tf))/(Tt-Tf)))/(1-2**(-10)) +((k*(Vt-Vf)*(Te-Tf))/(Tt-Tf))
        return Ve




import pandas as pd


chromatograms = pd.read_excel("C:/Users/nlancaster/OneDrive - UW-Madison/Astral_Phosphoproteomics_Manuscript/PythonFigures/ChromatogramExamples/PhosphoChromatogramsExamples_BPC.xlsx")
gradients = pd.read_excel("C:/Users/nlancaster/OneDrive - UW-Madison/Astral_Phosphoproteomics_Manuscript/PythonFigures/ChromatogramExamples/PhosphoChromatogramsExamples_Gradients.xlsx")


fig = plt.figure(figsize = (7, 8))
ax1 = plt.subplot(311)
ax2 = plt.subplot(312)
ax3 = plt.subplot(313)
ax1.spines[['top','right']].set_visible(False)
ax2.spines[['top','right']].set_visible(False)
ax3.spines[['top','right']].set_visible(False)




ax1g = ax1.twinx()
ax1g.spines['top'].set_visible(False)
ax2g = ax2.twinx()
ax2g.spines['top'].set_visible(False)
ax3g = ax3.twinx()
ax3g.spines['top'].set_visible(False)

times = list(gradients['Grad_7'])
mpb = list(gradients['MPB_7'])

t = [5]
times.append(5)
while t[-1]<=12.1:
    t.append(t[-1]+0.1)
    times.append(t[-1]+0.1)
for x in t:
    mpb.append(concave_grad(14,54,5,12.2,x,6))
times.append(max(times))
mpb.append(0)
gradient_df = pd.DataFrame({'t':times,'%B':mpb})
gradient_df = gradient_df.sort_values('t')    
ax1g.plot(gradient_df['t'],gradient_df['%B'],color = '#5084C4')
ax1g.set_ylim(0,)
ax1g.set_ylabel('%B Composition',fontweight = 'bold')


times = list(gradients['Grad_15'])
mpb = list(gradients['MPB_15'])
t = [5.2]
times.append(5.2)
while t[-1]<=20.1:
    t.append(t[-1]+0.1)
    times.append(t[-1]+0.1)
for x in t:
    mpb.append(concave_grad(12,52,5.2,20.2,x,6))
times.append(max(times))
mpb.append(0)
gradient_df = pd.DataFrame({'t':times,'%B':mpb})
gradient_df = gradient_df.sort_values('t')    
ax2g.plot(gradient_df['t'],gradient_df['%B'],color = '#5084C4')
ax2g.set_ylim(0,)
ax2g.set_ylabel('%B Composition',fontweight = 'bold')


times = list(gradients['Grad_30'])
mpb = list(gradients['MPB_30'])
t = [5]
times.append(5)
while t[-1]<=35.1:
    t.append(t[-1]+0.1)
    times.append(t[-1]+0.1)
for x in t:
    mpb.append(concave_grad(11,51,5,35.2,x,6))
times.append(max(times))
mpb.append(0)
gradient_df = pd.DataFrame({'t':times,'%B':mpb})
gradient_df = gradient_df.sort_values('t')    
ax3g.plot(gradient_df['t'],gradient_df['%B'],color = '#5084C4')
ax3g.set_ylim(0,)
ax3g.set_ylabel('%B Composition',fontweight = 'bold')




ax1.plot(chromatograms['RT_7'],chromatograms['BPC_7'],color = 'black',linewidth =1)
ax2.plot(chromatograms['RT_15'],chromatograms['BPC_15'],color = 'black',linewidth =1)
ax3.plot(chromatograms['RT_30'],chromatograms['BPC_30'],color = 'black',linewidth =1)


ax1.set_xlim(0,41)
ax2.set_xlim(0,41)
ax3.set_xlim(0,41)

ax1.set_ylim(0,)
ax2.set_ylim(0,)
ax3.set_ylim(0,)

ax1.set_ylabel('Base Peak Chromatogram',fontweight = 'bold')
ax2.set_ylabel('Base Peak Chromatogram',fontweight = 'bold')
ax3.set_ylabel('Base Peak Chromatogram',fontweight = 'bold')
ax3.set_xlabel('Retention Time (min)',fontweight = 'bold')


ax1.set_yticks([0,1e9,2e9,3e9,4e9,5e9,6e9])
ax2.set_yticks([0,1e9,2e9,3e9,4e9,5e9])
ax3.set_yticks([0,1e9,2e9,3e9,4e9,5e9])

ax1g.yaxis.label.set_color('#5084C4')
ax1g.tick_params(axis='y', colors='#5084C4')
ax1g.spines['right'].set_color('#5084C4')


ax2g.yaxis.label.set_color('#5084C4')
ax2g.tick_params(axis='y', colors='#5084C4')
ax2g.spines['right'].set_color('#5084C4')


ax3g.yaxis.label.set_color('#5084C4')
ax3g.tick_params(axis='y', colors='#5084C4')
ax3g.spines['right'].set_color('#5084C4')


plt.subplots_adjust(hspace =0.6)





plt.rcParams['svg.fonttype'] = 'none'

file_save_path= 'C:\\Users\\nlancaster\OneDrive - UW-Madison\Astral_Phosphoproteomics_Manuscript\PythonFigures\\'
# fig.savefig(file_save_path + '20230927_ChromatogramExamplesWGradient_BPC.svg')

