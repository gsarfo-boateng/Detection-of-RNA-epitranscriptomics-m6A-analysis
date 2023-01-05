#!/usr/bin/python
#Generate separate bed files for m6A sites depending on context
#Find more m6A sites from non-canonical peaks
#Similar to meme_from_deseq.py
#Note that current shifting uses GENOMIC sequences, not transcriptomic for conditions
#2 to 5
#Same as v2 but uses a less stringent initial maxpadj=0.1 to define "SIGNIFICANT"
#for the initial categorization of conditions and shiftings.
#Afterwards, then applies secondajd=0.05 to filter out the post-shifted sites by padj.
#This effectively removes any sites that are upstream of a not-so-significant Rm6AC site since they'd
#otherwise contribute to false-positives.

import read
import complimentary
import math
import numpy

minLFC=1.0
maxpadj=0.1#Might need to do more stringent maxpadj= for mettl4-related libs. 
secondadj=0.05
version_to_use='Aonly'
ref_genome='dm6'#######
upstream=10
seqlen=11
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


#Get list of chromosomes. Currently set at all chromosomes.
chromosomes=['all',]#Input here. Either specific chromosomes or 'all'.
chros=[]
y=open('store/STARdir/genomes/'+ref_genome+'/'+ref_genome+'_chromosomes.txt','r')
for line in y:
    chros+=[line.strip(),]
y.close()
if chromosomes[0] != 'all':
    chros = chromosomes


b=open('store/DESeq2/'+sample+'_temp_DESeq2_conditions_and_batch.txt','w')
total=0#For counting total number of starting peaks.

#Use output of DESeq2
a=open('store/DESeq2/'+sample+'_NoVsYes_DESeq2_results_condition_and_batch_'+version_to_use+'.txt','r')
next(a)
#Isolate lines that pass LFC and padj filters.
for line in a:
    ar=read.tab(line)
    padj=ar[6]
    if str(padj)!='NA':
        if (float(ar[2])>=minLFC) and (float(padj)<maxpadj):
            b.write(line)
            total+=1
a.close()
print ('total',total)

#Center is 10:11
#Condition 1: Anything that's not condition 2, 3, 4, 5.
#Condition 2: AAC at position 10:13; m6A at position 12 is SIGNIFICANT
#Condition 3: ARAC at position 10:14; m6A at position 13 is SIGNIFICANT
#Condition 4: Not RAC at position 9:12; A at position 10:11; RAC at position 12:15; m6A at position 14 is SIGNIFICANT
#Condition 5: Not RAC at position 9:12; A at position 10:11; NRAC at position 12:16; m6A at position 15 is SIGNIFICANT
conditions={}
#Coordinate shifts to be implemented for each condition.
adjust={}
for i in range(1,6):
    conditions[i]=[]
    adjust[i]=i-1

#Parse lines based on conditions.
for chro in chros:
    z=open('store/STARdir/genomes/genomes_parsed_by_chr/'+ref_genome+'/'+chro+'.fa','r')
    for line in z:
        genome=line.strip()
    z.close()
    passed=False
    b=open('store/DESeq2/'+sample+'_temp_DESeq2_conditions_and_batch.txt','r')
    for line in b:
        ar=read.tab(line)
        if chro == str(ar[0]).split(':')[0]:
            passed=True
            strand=str(ar[0])[-1]
            if strand == '+':
                start=int(str(ar[0])[0:-1].split(':')[1].split('-')[0])
                sequence=str.upper(genome[start-upstream:start+seqlen])           
            elif strand == '-':
                end=int(str(ar[0])[0:-1].split(':')[1].split('-')[1])
                sequence=str.upper(genome[end-seqlen:end+upstream])
                sequence=complimentary.seq(sequence)

            if sequence[10:13]=='AAC':
                conditions[2].append(line.strip())
            elif (sequence[10:14]=='AGAC' or sequence[10:14]=='AAAC'):
                conditions[3].append(line.strip())
            elif sequence[10:11]=='A' and sequence[9:12] not in ['GAC','AAC'] and sequence[12:15] in ['GAC','AAC']:
                conditions[4].append(line.strip())
            elif sequence[10:11]=='A' and sequence[9:12] not in ['GAC','AAC'] and sequence[13:16] in ['GAC','AAC']:
                conditions[5].append(line.strip())
            elif sequence[10:11] =='A':
                conditions[1].append(line.strip())
                
        elif ( chro != str(ar[0]).split(':')[0] ) and ( passed==True ) :
            break
    b.close()


#Shift lines from conditions2/3/4/5 to condition 1 if sites in the former don't have a corresponding significant counterpart in
#condition 1.
for i in range(2,6):
    to_delete=[]
    for count,line in enumerate(conditions[i]):
        ar=read.tab(line)
        chro = str(ar[0]).split(':')[0]
        strand=str(ar[0])[-1]
        start=int(str(ar[0])[0:-1].split(':')[1].split('-')[0])
        if strand=='+':
            start=start+adjust[i]
        else:
            start=start-adjust[i]
        end=start+1
        to_check=chro+':'+str(start)+'-'+str(end)+strand
        shift=True
        for q in range(len(conditions[1])):
            if to_check in conditions[1][q]:
                shift=False
                break
        if shift==True:
            conditions[1].append(line)
            to_delete.append(count)
    for c in sorted(to_delete,reverse=True):
        del conditions[i][c]

#Create bed files for each condition.
y=open('store/DESeq2/'+sample+'_condition1_m6A_peaks.bed','w')

#Make trackline for bed file.
if ref_genome == 'hg38':
    y.write('browser position chr6:31829816-31829965\n')
elif ref_genome == 'mm10':
    y.write('browser position chr6:31829816-31829965\n')
y.write('browser hide all\n'\
    +'track type=bed name="m6A condition'+str(i)+'" '\
    +'visibility=2 colorByStrand="255,0,0 0,0,255"\n')

for line in conditions[1]:
    ar=read.tab(line)
    pval_write=float(ar[6])
    #Following is when padj is so low that python registers it as '0.0'.
    if str(pval_write)=='0.0':
        print (line)
        pval_write=float(numpy.finfo(float).tiny)
    #
    if pval_write < secondadj:
        chro = str(ar[0]).split(':')[0]
        strand=str(ar[0])[-1]
        start=int(str(ar[0])[0:-1].split(':')[1].split('-')[0])
        end=start+1    
        y.write(chro+'\t'+str(start)+'\t'+str(end)+'\t'+\
        str(pval_write)+'\t1000\t'+\
        strand+'\n')

print (sample)
for i in range(1,6):
    print ('condition',i,len(conditions[i]))






    

