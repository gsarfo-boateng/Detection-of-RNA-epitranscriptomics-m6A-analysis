#!/usr/bin/python
#Make bedfile of writer-dependent m6A/m6Am sites with different logfold_cutoff.
#Need to do FDR_calculator.py to determine corresponeding pval_cutoff for each condition and each logfold_cutoff.

import read
import numpy as np
import math
import scipy.stats as stats




logfold_cutoff=2.0
pval_cutoff=0.09495834812337565#FDR 0.1





logfold_cutoff=2.0**logfold_cutoff
deseq2padj_cutoff=0.05#0.05 always
filename='Empty-E2_Mettl3-134B4'


c=open('Mettl3-dependent.bed','w')

a=open('store/quanti/'+filename+'_consolidate_output_annotated_contexted_appended.txt','r')
count=0
for line in a:
    ar=read.tab(line)
    filterPass=False
    for j in range(13,15):
        if ar[j]!='NA':
            if float(ar[j]) < deseq2padj_cutoff:
                filterPass=True
                break
    if filterPass==True:
        value1=[]
        for i in range(4,7):
            value1.append(float(ar[i]))
        value2=[]
        for i in range(7,10):
            value2.append(float(ar[i]))

        b=2.0#Just an arbitrary impossible pvalue.
        if ( float( np.mean(value1) ) >0.0 and float( np.mean(value2) ) ==0.0 ):
            b=stats.ttest_ind(value1,value2)[1]/2                
        elif ( ( np.mean(value1)/np.mean(value2) ) >= logfold_cutoff ):
            b=stats.ttest_ind(value1,value2)[1]/2

        #If PCIF1-dependent m6Am, write append coordinates.
        if b<pval_cutoff:
            count+=1
            score = math.log2(np.mean(value1)/np.mean(value2))
            c.write(ar[0]+'\t'+str(ar[1])+'\t'+str(ar[2])+'\t'+str(count)+'\t'+str(score)+'\t'+ar[3]+'\n')

                    
