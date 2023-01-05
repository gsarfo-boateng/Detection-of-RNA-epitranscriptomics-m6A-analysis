#!/usr/bin/python

#Annotate bed file reads with 
#NCBI RefSeq table of from UCSC:
#Go to UCSC genome browser > Tools > Table Browser > select correct genome and "NCBI Refseq" in track dropdown menu.
#0-based start, 1-based end.

#NOTE:
##For mRNA, only use transcripts with NM_ prefix in name column
##5_prime_UTR is anything upstream of CDS
##3_prime_UTR is anything downstream of CDS
##intron is anything outside of CDS, 5_prime_UTR and 3_prime_UTR but within txStart and tsEnd.
##In cases where there's more than 1 NM_ prefix transcript, use all but when annotating, give priority
##in this order: CDS>5'UTR,3'UTR,intron.
##e.g. if a site is found in CDS of transcript1 but intron in transcript2 of same gene, denote as CDS.

##For ncRNA, only use transcripts with NR_ prefix in name column
##v3 same as v2 but annotates unannotated sites that are X nucleotides upstream of 5'UTR.


def format_gff3(directory,ref_genome):
    #Format gff3 file to define 5'UTR, 3'UTR, introns and promoters    

    import read

    #First pass: for each gene, find most confident gene prefix.
    gffs={}    
    a=open(directory+ref_genome+'_NCBI_RefSeq.txt','r')
    next(a)
    for line in a:
        ar=read.tab(line)
        gene=ar[12]
        prefix=ar[1][0:2]
        
        try:
            gffs[gene].append(prefix)
        except KeyError:
            gffs[gene]=[prefix,]
    for gene in gffs.keys():
        #For mRNA
        if 'NM' in gffs[gene]:
            gffs[gene]='NM'
        #For ncRNA
        elif 'NR' in gffs[gene]:
            gffs[gene]='NR'
        #Not sure what YP is but it only tags mitochondrial-encoded mRNA.
        elif 'YP' in gffs[gene]:
            gffs[gene]='YP'
        #For predicted mRNA
        elif 'XM' in gffs[gene]:
            gffs[gene]='XM'
        #For predicted ncRNA
        elif 'XR' in gffs[gene]:
            gffs[gene]='XR'

        
    #Second pass, collect all lines that have prefix corresponding to highest confidence gene prefix for each gene, and collect under dictionary "temp".     
    temps={}
    a=open(directory+ref_genome+'_NCBI_RefSeq.txt','r')
    next(a)
    for line in a:
        ar=read.tab(line)
        gene=ar[12]
        prefix=ar[1][0:2]
        if prefix == gffs[gene]:
            try:
                temps[gene].append(line.strip())
            except KeyError:
                temps[gene]=[line.strip(),]
    a.close()

    #Define all regions and collect under "annotations" dictionary
    annotations={}
    #Generation sub-dictionaries.
    for gene in temps.keys():
        for line in temps[gene]:
            ar=read.tab(line)
            chro=ar[2]
            strand=ar[3]
            try:
                annotations[chro+strand][gene]={}                
            except KeyError:
                annotations[chro+strand]={}
                annotations[chro+strand][gene]={}
            for region in ['cds','3utr','5utr','intron','ncRNA','intron_ncRNA']:
                annotations[chro+strand][gene][region]=[]
    
    #Actual defining of regions.            
    for gene in temps.keys():
        for line in temps[gene]:
            ar=read.tab(line)
            chro=ar[2]
            strand=ar[3]
            if strand=='+':
                x='5utr'
                y='3utr'
            else:
                y='5utr'
                x='3utr'                
            txstart=int(ar[4])
            txend=int(ar[5])
            exonstarts=[]
            for br in ar[9].split(',')[0:-1]:
                exonstarts.append(int(br))
            exonends=[]
            for br in ar[10].split(',')[0:-1]:            
                exonends.append(int(br))

            #For mRNAs
            if gffs[gene] in ['NM','XM']:
                cdsstart=int(ar[6])
                cdsend=int(ar[7])
                #Single exon.
                if len(exonstarts)==1:
                    annotations[chro+strand][gene]['cds'].append([cdsstart,cdsend])
                    annotations[chro+strand][gene][x].append([txstart,cdsstart])
                    annotations[chro+strand][gene][y].append([cdsend,txend])
                #Multiple exons.
                else:
                    #Find out if any 5utr or 3utr is interrupted by introns; also denote introns.
                    for i in range(len(exonstarts)):
                        if exonstarts[i] <= cdsstart < exonends[i]:
                            leftutr=i
                        if exonstarts[i] < cdsend <= exonends[i]:
                            rightutr=i
                        if i>1:
                            annotations[chro+strand][gene]['intron'].append([exonends[i-1],exonstarts[i]])                            
                    #Denote cds.
                    for i in range(leftutr,rightutr+1):
                        if i==leftutr:
                            annotations[chro+strand][gene]['cds'].append([cdsstart,exonends[i]])
                        elif i==rightutr:
                            annotations[chro+strand][gene]['cds'].append([exonstarts[i],cdsend])
                        else:
                            annotations[chro+strand][gene]['cds'].append([exonstarts[i],exonends[i]])
                    #Denote left-most UTR.
                    for i in range(0,leftutr+1):
                        if i<leftutr:
                            annotations[chro+strand][gene][x].append([exonstarts[i],exonends[i]])
                        elif i==leftutr:
                            annotations[chro+strand][gene][x].append([exonstarts[i],cdsstart])
                    #Denote right-most UTR.
                    for i in range(rightutr,len(exonstarts)):
                        if i==rightutr:
                            annotations[chro+strand][gene][y].append([cdsend,exonends[i]])
                        if i>rightutr:
                            annotations[chro+strand][gene][y].append([exonstarts[i],exonends[i]])
                                
            #For ncRNAs.
            elif gffs[gene] in ['NR','XR']:
                #Single exon.
                if len(exonstarts)==1:
                    annotations[chro+strand][gene]['ncRNA'].append([txstart,txend])
                #Multiple exons.
                else:
                    for i in range(len(exonstarts)):
                        annotations[chro+strand][gene]['ncRNA'].append([exonstarts[i],exonends[i]])
                        if i>1:
                            annotations[chro+strand][gene]['intron_ncRNA'].append([exonends[i-1],exonstarts[i]])                            

            #For YP_ prefixes.                
            elif gffs[gene]=='YP':
                annotations[chro+strand][gene]['cds'].append([txstart,txend])

    del gffs            
    return annotations

    
###########START HERE
import read
ref_genome='hg38'#########

directory='store/annotations/'+ref_genome+'/'
priorities=['cds','3utr','5utr','intron','ncRNA','intron_ncRNA']
annotations=format_gff3(directory,ref_genome)
up=200#200nt upstream seems to cover most sites that are wrongly annotated as being upstream of TSS/5'UTR.

#####For single sample type.
##sample='Empty-E2'
##a=open('store/DESeq2/'+sample+'_condition'+condition+'_m6A_peaks.bed','r')
##for i in range(3):
##    next(a)
##condition='1'
##z=5#For strand index
##d=open('store/DESeq2/'+sample+'_condition'+condition+'_m6A_peaks_annotated.txt','w')
#####

###For paired sample quantitative comparisons.
##filename='Empty-E2_FTO-168C4'
##filename='Empty-E2_Mettl3-134B4'
##filename='Empty-E2_Mettl4-63A3'
##filename='Empty-E2_Alkbh5-7G4'
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
filename='Empty-E2+IAV_Pcif1-183B2+IAV'

a=open('store/quanti/'+filename+'_consolidate_output.txt','r')
z=3#For strand index
d=open('store/quanti/'+filename+'_consolidate_output_annotated.txt','w')
###


#Actual annotating.
sites={}
sitecount=0
for line in a:
    sitecount+=1
    sites[sitecount]=read.tab(line)
a.close()


for site in sites.keys():
    ar=sites[site]
    chro=ar[0]
    strand=ar[3]
    if chro+strand not in annotations.keys():
        sites[site].append('NA')
        sites[site].append('unannotated')
    else:
        end=int(ar[2])
        unannotated=True    
        for priority in priorities:
            for gene in annotations[chro+strand].keys():
                try:
                    for region in annotations[chro+strand][gene][priority]:
                        if region[0] < end <= region[1]:
                            sites[site].append(gene)
                            sites[site].append(priority)
                            unannotated=False
                            break
                except KeyError:
                    nothinghappens=True
                if unannotated==False:
                    break
            if unannotated==False:
                break
            
        if unannotated==True:
            for gene in annotations[chro+strand].keys():
                try:
                    for region in annotations[chro+strand][gene]['5utr']:
                        if ( strand=='+' and (region[0] -up) < end <= region[0] ) \
                           or \
                           ( strand=='-' and region[1] < end <= (region[1] +up) ):                        
                            sites[site].append(gene)
                            sites[site].append('5utr')
                            unannotated=False
                            break
                except KeyError:
                    nothinghappens=True
                if unannotated==False:
                    break
            
            if unannotated==True:
                sites[site].append('NA')
                sites[site].append('unannotated')
del annotations


#To output frequency of each annotation.
represent={}
represent['unannotated']=0
for priority in priorities:
    represent[priority]=0
for annotat in sites.values():
    if annotat[-1]=='unannotated':
        represent['unannotated']+=1
    else:        
        represent[annotat[-1]]+=1
    for ar in annotat:
        d.write(str(ar)+'\t')
    d.write('\n')

        

for key in represent.keys():
    print (str(key)+'\t'+str(represent[key]))

#Remove input file to save space.
import subprocess as sp               
prog='rm '
command=prog+\
         'store/quanti/'+filename+'_consolidate_output.txt'
sp.call(command,shell=True)
                
                                           
               
    
            
    

            
        
        
