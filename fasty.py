#!/usr/bin/python
#Scripts for fastq/fasta manipulation.

import subprocess as sp
import read
import numpy as np
from operator import itemgetter
import complimentary

def unzip_fastq(lib_names,origin,read_end,paired):
    #Unzip zipped fastq file in 'fastqgz/' directory into BT2/.

    for lib_name in lib_names:

        prog='gzip -dc '

        if paired==False:
            print ('Unzipping',lib_name+read_end[0],'...')
            
            sp.call(prog+' '+origin+lib_name+read_end[0]+'.fastq.gz > '+origin+lib_name+read_end[0]+\
                    '.fastq ', \
                     shell=True)

            print ('... '+lib_name+read_end[0]+' unzipped.')
            
        elif paired==True:
            print ('Unzipping',lib_name+read_end[0],' and ',lib_name+read_end[1],'...')
            
            sp.call(prog+' '+origin+lib_name+read_end[0]+'.fastq.gz > '+origin+lib_name+read_end[0]+\
                    '.fastq ', \
                     shell=True)
            sp.call(prog+' '+origin+lib_name+read_end[1]+'.fastq.gz > '+origin+lib_name+read_end[1]+\
                    '.fastq ', \
                     shell=True)
            
            print ('... '+lib_name+read_end[0]+' and '+lib_name+read_end[1]+' unzipped.')

        
def fastq_trimming(lib_names,origin,first,last,read_end,paired,modifier):
    #Trim fastq.clipped sequences

    for lib_name in lib_names:
        if paired==False:
            h=1
        elif paired==True:
            h=2
        for q in range(h):

            print ('Trimming ',lib_name+read_end[q],\
                  ' fastq and keeping from '+str(first)+' to '+str(last)+'nt...')

            b=open(origin+lib_name+read_end[q]+'_PE'+modifier+'_trimmed.fastq','w')
            a=open(origin+lib_name+read_end[q]+'_PE'+modifier+'.fastq','r')
            linecount=0
            for line in a:
                line=line.strip()
                linecount += 1
                if ( linecount == 2 ) or ( linecount == 4 ):
                    line=line[first-1:last]                    
                b.write(line+'\n')
                if linecount == 4:
                    linecount=0
            a.close()
            b.close()

            print ('... '+lib_name+read_end[q]+' trimmed.')


def fastq_clipping(lib_names,origin,clip_length,adapter1,adapter2,paired,read_end,polya,atail,ttail,qual,threads):
    
    #Clip fastq sequences

    #Set quality threshold and denote whether this is nextseq run or not.
    if qual[1]==True:
        qual_trim='--nextseq-trim '+str(qual[0])+' '
    elif qual[1]==False:
        qual_trim='-q '+str(qual[0])+' '

    for lib_name in lib_names:
        #Actual clipping.

        if paired==False:
            prog='miniconda/bin/cutadapt '
            print ('Clipping '+adapter1+' from '+str(lib_name)+read_end[0]+'.fastq with minimum length of '+\
                  str(clip_length)+'nt and '+qual_trim+'... ')

            command=prog+\
                     qual_trim+\
                     '-m '+\
                     str(clip_length)+' '+\
                     '-j '+\
                     str(threads)+' '+\
                     '-a '+\
                     adapter1+' '+\
                     '-o '+\
                     origin+lib_name+read_end[0]+'_clipped.fastq '+\
                     origin+lib_name+read_end[0]+'.fastq '
            sp.call(command,shell=True)

            #Delete unclipped fastq file to save space.
            prog='rm '
            command=prog+\
                     origin+lib_name+read_end[0]+'.fastq '
            sp.call(command,shell=True)

        elif paired==True:
            print ('Clipping '+adapter1+' from '+str(lib_name)+read_end[0]+' and '+\
                  adapter2+' from '+str(lib_name)+read_end[1]+' '+\
                  'fastqs with minimum length of '+\
                  str(clip_length)+'nt and '+qual_trim+'... ')
            
            prog='miniconda/bin/cutadapt '
            if polya == True:
                #First remove regular adapters
                command=prog+\
                         qual_trim+\
                         '-m '+\
                         str(clip_length)+' '+\
                         '-j '+\
                         str(threads)+' '+\
                         '-a '+\
                         adapter1+' '+\
                         '-A '+\
                         adapter2+' '+\
                         '-o '+\
                         origin+lib_name+read_end[0]+'_PE_temp.fastq '+\
                         '-p '+\
                         origin+lib_name+read_end[1]+'_PE_temp.fastq '+\
                         origin+lib_name+read_end[0]+'.fastq '+\
                         origin+lib_name+read_end[1]+'.fastq '
                sp.call(command,shell=True)
                
                #Then remove polyA tail in R1 3' and polyT in R2 5'.
                command=prog+\
                         '-m '+\
                         str(clip_length)+' '+\
                         '-j '+\
                         str(threads)+' '+\
                         '-a '+\
                         atail+' '+\
                         '-G '+\
                         ttail+' '+\
                         '-o '+\
                         origin+lib_name+read_end[0]+'_PE_clipped.fastq '+\
                         '-p '+\
                         origin+lib_name+read_end[1]+'_PE_clipped.fastq '+\
                         origin+lib_name+read_end[0]+'_PE_temp.fastq '+\
                         origin+lib_name+read_end[1]+'_PE_temp.fastq '
                sp.call(command,shell=True)

                #Delete unclipped fastq file to save space.
                prog='rm '
                command=prog+\
                         origin+lib_name+read_end[0]+'_PE_temp.fastq '
                sp.call(command,shell=True)
                command=prog+\
                         origin+lib_name+read_end[1]+'_PE_temp.fastq '
                sp.call(command,shell=True)


            else:
                command=prog+\
                         qual_trim+\
                         '-m '+\
                         str(clip_length)+' '+\
                         '-j '+\
                         str(threads)+' '+\
                         '-a '+\
                         adapter1+' '+\
                         '-A '+\
                         adapter2+' '+\
                         '-o '+\
                         origin+lib_name+read_end[0]+'_PE_clipped.fastq '+\
                         '-p '+\
                         origin+lib_name+read_end[1]+'_PE_clipped.fastq '+\
                         origin+lib_name+read_end[0]+'.fastq '+\
                         origin+lib_name+read_end[1]+'.fastq '
                sp.call(command,shell=True)
            
            #Delete unclipped fastq file to save space.
            prog='rm '
            command=prog+\
                     origin+lib_name+read_end[0]+'.fastq '
            sp.call(command,shell=True)
            command=prog+\
                     origin+lib_name+read_end[1]+'.fastq '
            sp.call(command,shell=True)

        print ('... clipping done.')



def fastq_UMI(lib_names,origin,paired,readwithUMI,readwithoutUMI,UMIposition,modifier,read_end):
    #Remove UMIs and adds the UMI to the sequence name in the fastq.
    #Note that there's an unusually higher percentage of 'GGGGGGGG' UMIs but
    #these are likely due to sequencing machine base-calling errors. Regardless,
    # they are removed as these don't map well to human genome (so problem solves itself).

    for lib_name in lib_names:

        if paired==False:
            print ('Transfering UMI at position'+\
                  str(UMIposition[0]),'to',\
                  str(UMIposition[1])+',',\
                  'to sequence name for',lib_name+readwithUMI)
            end=''
        if paired==True:
            print ('Transfering UMI at position'+\
                  str(UMIposition[0]),'to',\
                  str(UMIposition[1])+',',\
                  'to sequence name for',lib_name+readwithUMI,'and processing',\
                  readwithoutUMI)
            end='_PE'

        UMIs=[]
        b=open(origin+lib_name+readwithUMI+end+modifier+'_identified.fastq','w')
        a=open(origin+lib_name+readwithUMI+end+modifier+'.fastq','r')
        linecount=0
        seqs=[]
        for line in a:
            line=line.strip()
            linecount += 1
            seqs+=[line,]
            if linecount == 4:
                UMI=seqs[1][UMIposition[0]-1:UMIposition[1]]
                UMIs.append(UMI)
                seqs[0]=seqs[0].split(' ')[0]+':'+UMI+' '+seqs[0].split(' ')[1]
                seqs[1]=seqs[1][UMIposition[1]::]
                seqs[3]=seqs[3][UMIposition[1]::]
                for seq in seqs:
                    b.write(seq+'\n')
                linecount=0
                seqs=[]
        a.close()
        b.close()

        c=open(origin+lib_name+'_UMIstats.txt','w')
        UMIfreq={}
        for UMI in UMIs:
            try:
                UMIfreq[UMI]+=1
            except KeyError:
                UMIfreq[UMI]=1
        for value in sorted(UMIfreq.items(),key=itemgetter(1),reverse=True):
            c.write(value[0]+'\t'+str(value[1])+'\n')
        c.close()

        #Trim 3'end of read2 accordingly
        #See how much of 3' end matches complement of UMI and trim.
        if paired==True:
            d=open(origin+lib_name+readwithoutUMI+end+modifier+'_identified.fastq','w')
            e=open(origin+lib_name+readwithoutUMI+end+modifier+'.fastq','r')
            linecount=0
            UMIcount=0
            seqs=[]
            for line in e:
                line=line.strip()
                linecount += 1
                seqs+=[line,]
                if linecount == 4:
                    seqs[0]=seqs[0].split(' ')[0]+':'+UMIs[UMIcount]+' '+seqs[0].split(' ')[1]
                    umi_rc=complimentary.seq(UMIs[UMIcount])
                    for h in [8,7,6,5,4,3,2,1]:
                        if seqs[1][-h::]==umi_rc[0:h]:
                            seqs[1]=seqs[1][0:-h]
                            seqs[3]=seqs[3][0:-h]
                            break                        
                    for seq in seqs:
                        d.write(seq+'\n')
                    linecount=0
                    seqs=[]
                    UMIcount+=1
            d.close()
            e.close()

            #Delete clipped but not identified fastq file to save space.
            if paired==True:
                prog='rm '
                command=prog+\
                         origin+lib_name+read_end[0]+'_PE_clipped.fastq '
                sp.call(command,shell=True)
                command=prog+\
                         origin+lib_name+read_end[1]+'_PE_clipped.fastq '
                sp.call(command,shell=True)

            

        print ('... '+lib_name+' UMI transferred.')






    
def fastq_nucl_dist(lib_names,origin,modifier,maxlength,read_end,paired):
    #Extract nucleotide sequences from fastq files and plot
    #nucleotide distribution.
    
    for lib_name in lib_names:
        if paired==False:
            h=1
        elif paired==True:
            h=2
        for q in range(h):
            print ('Extracting sequence of '+lib_name+read_end[q]+modifier+'.fastq...')

            #Extracting sequences in fastq.
            if paired==False:
                a=open(origin+lib_name+read_end[q]+modifier+'.fastq','r')
            elif paired==True:
                a=open(origin+lib_name+read_end[q]+'_PE'+modifier+'.fastq','r')                
            seqs=[]
            linecount=0
            for line in a:
                linecount +=1
                if linecount ==2:
                    line=line.strip()[0:maxlength]
                    seqs += [line, ]
                if linecount == 4:
                    linecount= 0
            a.close()
            seqs=sorted(seqs)
            collapseds={}
            prev=''
            for seq in seqs:
                if seq==prev:
                    collapseds[seq] += 1.0
                else:
                    collapseds[seq] = 1.0
                prev=seq
            del seqs

            print ('... sequences extracted.')
            print ('Plotting nucleotide distribution of '+lib_name+read_end[q]+modifier+'.fastq...')

            #Calculating nucleotide distribution.
            dist={}
            for i in ['A','C','G','T','N','all']:
                dist[i]={}
                for j in range(1,maxlength+1):
                    dist[i][j]=0.0
            
            for collapsed in collapseds.keys():
                for i in list(enumerate(collapsed,start=1)):
                    dist['all'][i[0]] += collapseds[collapsed]
                    dist[i[1]][i[0]] += collapseds[collapsed]
            del collapseds
            
            for_plot={}
            alls=[]
            for j in range(1,maxlength+1):
                alls+= [dist['all'][j], ]
            for i in ['A','C','G','T','N']:
                to_div=[]
                for j in range(1,maxlength+1):
                    to_div+= [dist[i][j], ]
                for_plot[i]=list( np.divide(to_div,alls) )
            del dist
            
            #Actual plotting.
            import matplotlib
            matplotlib.use('Agg')#Agg (low quality) or SVG (high quality).
            import matplotlib.pyplot as plt
            N=maxlength
            ind=np.arange(N) #The x locations for the groups.
            width=0.85 #The width of the bars; can also be len(x) sequence.
            pA=plt.bar(ind,for_plot['A'],width,color='r')
            pC=plt.bar(ind,for_plot['C'],width,color='b',bottom=for_plot['A'])
            nextbottom=list( np.add(for_plot['A'],for_plot['C']) )
            pG=plt.bar(ind,for_plot['G'],width,color='y',bottom=nextbottom)
            nextbottom=list( np.add(nextbottom,for_plot['G']) )
            pT=plt.bar(ind,for_plot['T'],width,color='g',bottom=nextbottom)
            nextbottom=list( np.add(nextbottom,for_plot['T']) )
            pN=plt.bar(ind,for_plot['N'],width,color='m',bottom=nextbottom)
            plt.ylabel('Fraction')
            ind=np.arange(0,N+1,10)
            xlabeling=[]
            for j in range(0,N+1,10):
                xlabeling += [str(j), ]
            xlabeling=tuple(xlabeling)
            plt.xticks(ind-width/2., xlabeling )
            plt.yticks(np.arange(0.0,1.01,0.25))
            plt.legend( (pA[0],pC[0],pG[0],pT[0],pN[0]), ('A','C','G','T','N') )
            plt.savefig(origin+lib_name+read_end[q]+modifier+'_nucldist')        
            del for_plot

            print ('...nucleotide distribution plotted.')


def fastq_length(lib_names,origin,modifier,read_end,paired):
    #Determine length distribution of fastq sequences and provide raw
    #data for plotting histogram.
    
    for lib_name in lib_names:
        if paired==False:
            h=1
        elif paired==True:
            h=2
        for q in range(h):

            print ('Extracting sequence of '+lib_name+read_end[q]+modifier+'.fastq...')

            #Extracting sequences in fastq.
            a=open(origin+lib_name+read_end[q]+modifier+'.fastq','r')
            lengths={}
            linecount=0
            for line in a:
                linecount +=1
                if linecount ==2:
                    length=len(line.strip())
                    try:
                        lengths[length]+=1
                    except KeyError:
                        lengths[length]=1
                if linecount == 4:
                    linecount= 0
            a.close()

            #This is just to calculate percentage of total clipped reads that have lengths >=max read length.
            k=0.0
            for i in lengths.values():
                k+=float(i)
            print (lib_name, str(float(max(sorted(lengths.values())) ) / k))

            for ii in sorted(lengths.keys()):
                print (str(ii)+'\t'+str(lengths[ii]))

    


