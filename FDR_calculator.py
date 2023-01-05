#!/usr/bin/python
#Extract pvalue corresponding to FDR
#based on https://en.wikipedia.org/wiki/False_discovery_rate#Benjamini.E2.80.93Hochberg_procedure
#Doesn't work for methylases that only target few sites e.g. mettl4.
#Requires output of m6A_quant_v7.py or one of its subsequent output files.

import read
import math
import numpy as np
import scipy.stats as stats
import math

FDR=0.1
cutoff=2.0#this is log2fc
cutoff=2.0**cutoff
#1 for mono-directional (condition1>condition2), -1 for mono-directional (condition1<condition2), 2 for bi-directional (condition1>condition2 or condition1<condition2)
directional=2

z=open('../../Downloads/Empty-E2_Mettl3-134B4_consolidate_output_annotated_contexted_appended.txt','r')
##z=open('../../Downloads/Empty-E2_Mettl4-63A3_consolidate_output_annotated_contexted_appended_without+1.txt','r')
##z=open('../../Downloads/Mettl4-63A3-OE-WT-Mettl4_Mettl4-63A3-OE-CD-Mettl4_consolidate_output_annotated_contexted_appended.txt','r')
##z=open('../../Downloads/Empty-E2_Pcif1-183B2_consolidate_output_annotated_contexted_appended.txt','r')
##z=open('../../Downloads/Empty-E2-2_Pcif1-183B2_consolidate_output_annotated_contexted_appended.txt','r')
##z=open('../../Downloads/Empty-E2+IAV_Pcif1-183B2+IAV_consolidate_output_annotated_contexted_appended.txt','r')
##z=open('../../Downloads/WT-S2_Mettl3-KO-S2_consolidate_output.txt','r')
##z=open('../../Downloads/W1118-ovary_Mettl3-KO-ovary_consolidate_output.txt','r')


#Actual running of script
if directional==1:
    Ps=[]
    for line in z:
        ar=read.tab(line)
        value1=[]
        for i in range(4,7):
            value1.append(float(ar[i]))
        value2=[]
        for i in range(7,10):
            value2.append(float(ar[i]))
        if ( float( np.mean(value1) ) > 0.0 and float( np.mean(value2) ) ==0.0 ) or \
           ( ( np.mean(value1)/np.mean(value2) ) >= cutoff ):
            Ps += [stats.ttest_ind(value1,value2)[1]/2,]#This is divided by 2 because I'm looking for mono-directional shift.

elif directional==-1:
    Ps=[]
    for line in z:
        ar=read.tab(line)
        value1=[]
        for i in range(4,7):
            value1.append(float(ar[i]))
        value2=[]
        for i in range(7,10):
            value2.append(float(ar[i]))
        if ( float( np.mean(value2) ) > 0.0 and float( np.mean(value1) ) ==0.0 ) or \
           ( ( np.mean(value2)/np.mean(value1) ) >= cutoff ):
            Ps += [stats.ttest_ind(value1,value2)[1]/2,]#This is divided by 2 because I'm looking for mono-directional shift.

elif directional==2:
    Ps=[]
    for line in z:
        ar=read.tab(line)
        value1=[]
        for i in range(4,7):
            value1.append(float(ar[i]))
        value2=[]
        for i in range(7,10):
            value2.append(float(ar[i]))
        if ( float( np.mean(value1) ) > 0.0 and float( np.mean(value2) ) ==0.0 ) or \
           ( ( np.mean(value1)/np.mean(value2) ) >= cutoff ) or \
           ( float( np.mean(value2) ) > 0.0 and float( np.mean(value1) ) ==0.0 ) or \
           ( ( np.mean(value2)/np.mean(value1) ) >= cutoff ):
            Ps += [stats.ttest_ind(value1,value2)[1],]#Not divided by 2 because I'm looking for bi-directional shift.
            
#Actual calculation of FDR-based p-val cutoff.
m=len(Ps)#m is total number of called peaks to be assessed.
k='NA'#k is to eventually be final number of peaks that pass FDR   
for prob in enumerate(sorted(Ps),start=1):
    if ( float(prob[1]) ) > ( float(prob[0])*FDR/m ):
        k=int(prob[0])-1.0
        break
if k=='NA':#'NA' means all sites passed FDR. 0 means no sites passed.
    k=m
if k==0.0:
    cutoff_report='NA'
else:
    cutoff_report=sorted(Ps)[int(k)-1]
print ('k',int(k),'m',int(m),'FDR',FDR,'pvalue',cutoff_report)

