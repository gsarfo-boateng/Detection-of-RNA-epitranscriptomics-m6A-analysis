#!/usr/bin/python
#Filter out all sites that do do not look like single sharp peaks.
#Single sharp peak denoted as DESeq2 padj of the site being at least 1,000-fold
#lower than the sites at -2,-1,+1,+2 of the site being checked. 

import read
import complimentary
import math
import numpy

sigFoldCutoff=1000.0

##sample='Mettl3-134B4'
##sample='Empty-E2'
##sample='FTO-168C4'
##sample='Alkbh5-7G4'
##sample='WT-nucl'
##sample='Mettl4-63A3'
##sample='Pcif1-A5'
##sample='Scr-B9'
##sample='Mettl16-A2'
##sample='Hek-cyto'
##sample='Hek-nucl'
##sample='siRNA-screen'
##sample='Pcif1-183B2'
##sample='Mettl4-63A3-OE-WT-Mettl4'
##sample='Mettl4-63A3-OE-CD-Mettl4'
##sample='Empty-E2-2'
##sample='Mettl3-134B4-2'
##sample='Empty-E2-OE-FL-Fto'
##sample='Mettl5-71C2'
##sample='Empty-E2-abcam'
##sample='Empty-E2-sysy202111'
##sample='Hemk2-165D2'
##sample='HCT116'
##sample='9DIV'
##sample='19DIV'
##sample='Empty-E2+IAV'
##sample='Pcif1-183B2+IAV'
sample='W1118-ovary'
##sample='Mettl3-KO-ovary'
##sample='WT-S2'
##sample='Mettl3-KO-S2'



#Get list of m6As to be queried from condition1
querys={}
y=open('store/DESeq2/'+sample+'_condition1_m6A_peaks.bed','r')
for i in range(3):
    next(y)
for line in y:
    ar=read.tab(line)
    chro=ar[0]
    end=int(ar[2])
    strand=ar[5]
    querys[chro,end,strand]=line
print ('before',len(querys.keys()))
y.close()

#Get list of padj from DESeq2 output.
b=open('store/DESeq2/'+sample+'_NoVsYes_DESeq2_results_condition_and_batch_all.txt','r')
next(b)
padjs={}
for line in b:
    ar=read.tab(line)
    if ar[6]!='NA':
        strand=ar[0][-1]
        ar[0]=ar[0][0:-1]
        chro=ar[0].split(':')[0]
        end=int(ar[0].split('-')[1])
        if chro not in padjs.keys():
            padjs[chro]={}
            for k in ['+','-']:
                padjs[chro][k]={}
        padjs[chro][strand][end]=float(ar[6])            
b.close()

to_writes=[]
adjusts=[-2,-1,1,2]    
for query in querys.keys():
    chro=query[0]
    end=query[1]
    strand=query[2]
    centerpadj=padjs[chro][strand][end]
    sides=[]
    for adjust in adjusts:
        side=end+adjust
        sides.append(side)
    passed=1
    for side in sides:
        if centerpadj==0.0:
            #If centerpadj is really so significant, just pass it.
            passed=1
        else:
            try:
                if ( padjs[chro][strand][side]/centerpadj ) < sigFoldCutoff:
                    passed=0
                    break
            except KeyError:
                nothinghappens=1
    if passed==1:
        to_writes.append(querys[query])            
print ('after',len(to_writes))

w=open('store/DESeq2/'+sample+'_good_m6A_peaks.bed','w')
for to_write in to_writes:
    w.write(to_write)

