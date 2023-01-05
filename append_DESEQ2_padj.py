#!/usr/bin/python
#After tag_motif.py, append the DESeq2 padj for condition 1 and 2 respectively
#as 2 new columns.

import read

version_to_use='Aonly'#'all' or 'Aonly'

#Create dictionary of padj from the 2 paired DESeq2 outputs.
##filename='Empty-E2_FTO-168C4'
##filename='Empty-E2_Mettl3-134B4'
##filename='Empty-E2_Mettl4-63A3'
##filename='Empty-E2_Pcif1-A5'
##filename='Empty-E2_Alkbh5-7G4'
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

samples=filename.split('_')
deseq={}
for sample in samples:
    deseq[sample]={}
    #Use output of DESeq2
    a=open('store/DESeq2/'+sample+'_NoVsYes_DESeq2_results_condition_and_batch_'+version_to_use+'.txt','r')
    next(a)
    for line in a:
        ar=read.tab(line)
        if ar[6]=='NA':
            deseq[sample][ar[0]]='NA'
        else:
            deseq[sample][ar[0]]=ar[6]
    a.close()

c=open('store/quanti/'+filename+'_consolidate_output_annotated_contexted_appended.txt','w')
b=open('store/quanti/'+filename+'_consolidate_output_annotated_contexted.txt','r')
for line in b:
    ar=read.tab(line)
    chro=ar[0]
    start=ar[1]
    end=ar[2]
    strand=ar[3]
    key=chro+':'+str(start)+'-'+str(end)+strand
    try:
        cond1=deseq[samples[0]][key]
    except KeyError:
        cond1='NA'
    try:
        cond2=deseq[samples[1]][key]
    except KeyError:
        cond2='NA'
    c.write(line.strip()+'\t'+str(cond1)+'\t'+str(cond2)+'\n')


#Remove input file to save space.
import subprocess as sp               
prog='rm '
command=prog+\
         'store/quanti/'+filename+'_consolidate_output_annotated_contexted.txt'
sp.call(command,shell=True)

print ('append done')
    

