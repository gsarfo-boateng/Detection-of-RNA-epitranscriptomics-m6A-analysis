#!/usr/bin/python
#Given significant differential methylation sites
#output fasta file for meme-chip and gene list for GO analysis
#Execute on aquila.
#Make sure to first mkdir 'metagene' in /store/

import read
import numpy as np
from scipy import stats
import motif_generation
import math

#Either 'over', 'under' or 'all'.
#Note for 'over' and 'under', anything without triplicate non-'NA' values in m6A levels will be filtered out.

#Mohua: Use 'all' if you want to see both directions. Use 'over' if you want to see only the hits that have positive RML fold change and vice versa (depending on the direction you defined)
#Technically speaking, you shoudl alter the pval_cutoff if you are only interested in over or under. If you're confused about this, we can do a quick call about it.
filter_type='all'

deseq2padj_cutoff=0.05#0.05 always
pval_cutoff=0.08606434852035316#Mohua, input pval_cutoff calculated from FDR_calculator.py. 
logfold_cutoff=1.0#Mohua, make sure this value matches the cutoff you used for FDR_calculator.py.
logfold_cutoff=2.0**logfold_cutoff

ref_genome='hg38'
upstream=250
seqlen=250


#This part is for Pcif1-dependent m6Am 
filename='Empty-E2_Pcif1-183B2'#Mohua, input your filename for each library here.

z=open('store/quanti/'+filename+'_gene_list.txt','w')#For GO analysis.
a=open('store/quanti/'+filename+'_consolidate_output_annotated_contexted_appended.txt','r')
coors={}
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
        if filter_type=='over':
            if ( float( np.mean(value1) ) >0.0 and float( np.mean(value2) ) ==0.0 ):
                b=stats.ttest_ind(value1,value2)[1]/2                
            elif ( ( np.mean(value1)/np.mean(value2) ) >= logfold_cutoff ):
                b=stats.ttest_ind(value1,value2)[1]/2
        elif filter_type=='under':
            if ( float( np.mean(value2) ) >0.0 and float( np.mean(value1) ) ==0.0 ):
                b=stats.ttest_ind(value1,value2)[1]/2                
            elif ( ( np.mean(value2)/np.mean(value1) ) >= logfold_cutoff ):
                b=stats.ttest_ind(value1,value2)[1]/2

        #Write into gene_name output file if satisfy fold change and pval cutoffs.
        #Also generate coordinates for motif_generation.
        chro=ar[0]
        strand=ar[3]
        if chro not in coors.keys():
            coors[chro]={}
            coors[chro]['+']=[]
            coors[chro]['-']=[]
        if strand=='+':
            coor=int(ar[1])
        elif strand=='-':
            coor=int(ar[2])
            
        if filter_type=='over' or filter_type=='under' :
           if b<pval_cutoff:               
                z.write(ar[10]+'\n')
                coors[chro][strand].append(coor)
        elif filter_type=='all':
            z.write(ar[10]+'\n')
            coors[chro][strand].append(coor)

a.close()
z.close()

motif_generation.fasta(ref_genome,upstream,seqlen,coors,filename)

#This part is for CAGE_subset sites not within +/5kb of m6Am-lfc1 sites.
filename='CAGE_subset'
c=open('store/metagene/CAGE_subset.txt','r')
coors={}
for line in c:
    ar=read.tab(line)
    chro=ar[0]
    if chro not in coors.keys():
        coors[chro]={}
        coors[chro]['+']=[]
        coors[chro]['-']=[]
    strand=ar[2]
    coor=int(ar[1])
    if strand=='-':
        coor-=1
    coors[chro][strand].append(coor)

motif_generation.fasta(ref_genome,upstream,seqlen,coors,filename)


