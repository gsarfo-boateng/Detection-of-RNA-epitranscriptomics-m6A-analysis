#!/usr/bin/python
#Add motif context to annoteted m6Aquant outputs.

import read
import numpy
from scipy import stats
import motif_generation

ref_genome='dm6'
upstream=2
seqlen=3
##filename='Empty-E2_Mettl3-134B4'
##filename='Empty-E2_FTO-168C4'
##filename='Empty-E2_Alkbh5-7G4'
##filename='WT-nucl_WT-cyto'
##filename='Empty-E2_Mettl4-63A3'
##filename='Empty-E2_Pcif1-A5'
##filename='Scr-B9_Mettl16-A2'
##filename='Scr-B9_Pcif1-A5'
##filename='Hek-cyto_Hek-nucl'
##filename='Empty-E2_Mettl16-A2'
##filename='siRNA-screen'
##filename='crispr-screen'
##filename='Empty-E2_Pcif1-183B2'
##filename='Mettl4-63A3-OE-WT-Mettl4_Mettl4-63A3-OE-CD-Mettl4'
##filename='Empty-E2_Empty-E2'#Strictly for plotting metagene and memechip. RML values are fake fill-ins.
##filename='Empty-E2-2_Mettl3-134B4-2'
##filename='Quantitative-Mix'#Can't actually be used. Just do empty-e2-2 veruss mettl3-134b2-2 then change name afterwards.
##filename='Empty-E2_Empty-E2-2'
##filename='Empty-E2_Empty-E2-OE-FL-Fto'
##filename='Empty-E2_Mettl5-71C2'
##filename='Empty-E2_Empty-E2-abcam'
##filename='Empty-E2_Empty-E2-sysy202111'
##filename='Empty-E2_Hemk2-165D2'
##filename='Empty-E2_Empty-E2'
##filename='HCT116_HCT116'
##filename='Empty-E2-2_Pcif1-183B2'
##filename='9DIV_19DIV'
##filename='Empty-E2+IAV_Pcif1-183B2+IAV'
filename='W1118-ovary_Mettl3-KO-ovary'
filename='WT-S2_Mettl3-KO-S2'



a=open('store/quanti/'+filename+'_consolidate_output_annotated.txt','r')
coors={}
for line in a:
    ar=read.tab(line)
    chro=ar[0]
    if chro not in coors.keys():
        coors[chro]={}
        coors[chro]['+']=[]                
        coors[chro]['-']=[]                
    strand=ar[3]
    if strand=='+':
        coor=int(ar[1])
    elif strand=='-':
        coor=int(ar[2])
    coors[chro][strand].append(coor)
a.close()

contexts=motif_generation.append(ref_genome,upstream,seqlen,coors,filename)


a=open('store/quanti/'+filename+'_consolidate_output_annotated.txt','r')
b=open('store/quanti/'+filename+'_consolidate_output_annotated_contexted.txt','w')
for line in a:
    ar=read.tab(line)
    chro=ar[0]
    strand=ar[3]
    if strand=='+':
        coor=int(ar[1])
    elif strand=='-':
        coor=int(ar[2])
    context=contexts[chro][strand][coor]
    b.write(line.strip()+'\t'+context+'\n')
    
#Remove input file to save space.
import subprocess as sp               
prog='rm '
command=prog+\
         'store/quanti/'+filename+'_consolidate_output_annotated.txt'
sp.call(command,shell=True)

print ('tag done')
