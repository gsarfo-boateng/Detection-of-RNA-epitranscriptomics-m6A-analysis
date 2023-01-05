#!/usr/bin/python
#Scripts for bowtie2 mapping, then converts sam output to
#sorted bam file.

import subprocess as sp
import read
import complimentary


def pairing(paired,lib_name,destination,modifier,read_end,origin):
    #For paired end mapping.
    if paired==True:
        to_run=['-1 '+origin+lib_name+read_end[0]+'_PE'+modifier+'.fastq '+\
                '-2 '+origin+lib_name+read_end[1]+'_PE'+modifier+'.fastq ',\
                '_PE']
    #For single end mapping.
    elif paired==False:
        to_run=[origin+lib_name+read_end[0]+modifier+'.fastq ',\
                '_SE']
    return to_run


def run_bowtie2(lib_names,destination,\
                modifier,ref_genome,\
                mismatch,paired,threads,read_end,origin):
    #Running bowtie2.

    for lib_name in lib_names:
    
        prog='bowtie2 '

        to_run=pairing(paired,lib_name,destination,modifier,read_end,origin)

        if paired==True:
            #For PE.
            print ('Running bowtie2 on '+lib_name+read_end[0]+'_PE'+modifier+'.fastq and '+\
                  lib_name+read_end[1]+'_PE'+modifier+'.fastq...')


            command=prog+\
                     '-q '+\
                     '-N '+str(mismatch)+' '+\
                     '--no-mixed '+\
                     '--no-discordant '+\
                     '-p '+str(threads)+' '+\
                     '-x store/BT2/bt_files/'+ref_genome+' '+\
                     to_run[0]+\
                     '-S '+destination+lib_name+modifier+to_run[1]+'.sam '
        elif paired==False:
            #For SE.
            print ('Running bowtie2 on '+lib_name+read_end[0]+modifier+'.fastq...')

            command=prog+\
                     '-q '+\
                     '-N '+str(mismatch)+' '+\
                     '--dovetail '+\
                     '--no-mixed '+\
                     '--no-discordant '+\
                     '-p '+str(threads)+' '+\
                     '-x store/BT2/bt_files/'+ref_genome+' '+\
                     to_run[0]+\
                     '-S '+destination+lib_name+modifier+to_run[1]+'.sam '
        sp.call(command,shell=True)

        print ('...bowtie2 done.')


def sam_to_sorted_bam(lib_names,destination,\
                      modifier,ref_genome,\
                      paired,read_end,origin):
    #See title.

    for lib_name in lib_names:

        print ('Converting '+lib_name+modifier+'.sam to sorted bam...')

        to_run=pairing(paired,lib_name,destination,modifier,read_end,origin)

        #Convert sam to bam here.
        prog='samtools view -Sb'
        command=prog+' '+\
                 destination+lib_name+modifier+to_run[1]+'.sam '+\
                 '-o '+\
                 destination+lib_name+modifier+to_run[1]+'.bam '
        sp.call(command,shell=True)   

        #Delete sam file to save space.
        prog='rm '
        command=prog+\
                 destination+lib_name+modifier+to_run[1]+'.sam '
        sp.call(command,shell=True)
        
        #Sort bam.
        prog='samtools sort'
        command=prog+' '+\
                 destination+lib_name+modifier+to_run[1]+'.bam '+\
                 destination+lib_name+modifier+to_run[1]+'_sorted '
        sp.call(command,shell=True)
        
        #Delete unsorted bam file to save space.
        prog='rm '
        command=prog+\
                 destination+lib_name+modifier+to_run[1]+'.bam '
        sp.call(command,shell=True)

        print ('... conversion done.')
        


def remove_PCR_duplicates(lib_names,destination,\
                          modifier,ref_genome,\
                          paired):

    #NOT USED SINCE RMDUP WITHOUT USING UMI RESULTS IN OVER-REMOVAL
    
    #Remove PCR duplicates. Needs sorted bam.
    #THIS VERSION IS ONLY FOR PAIR-END READS!
    #This version does NOT use UMIs.
    
    for lib_name in lib_names:

        print ('Removing PCR duplicates for '+lib_name+'...')
    
        to_run=pairing(paired,lib_name,destination,modifier)

        #Removing duplicates here.
        prog='samtools rmdup '
        command=prog+' '+\
                 destination+lib_name+modifier+to_run[1]+'_sorted.bam '+\
                 destination+lib_name+modifier+to_run[1]+'_rmdup.bam '
        sp.call(command,shell=True)

        print ('...duplicates removed.')


def indexing_bam(lib_names,destination,\
                 modifier,ref_genome,\
                 paired):
    #Use samtools to index sorted.bam.

    for lib_name in lib_names:

        print ('Indexing sorted.bam of '+lib_name+'...')

        to_run=pairing(paired,lib_name,destination,modifier)

        #Generate sorted.bam index to be used in samtools idxstats or
        #preprocessing of MACE.
        prog='samtools index '
        command=prog+destination+lib_name+modifier+to_run[1]+'_sorted.bam '
        sp.call(command,shell=True)

    
def total_mapped_reads(lib_name,destination,modifier,paired):
    #NOT USED SINCE COUNTING TOTAL READS WITHOUT USING UMI RESULTS IN OVER-COUNTING
    
    #Calculate total number of million mapped reads.
    #DOING NON-RMDUP VERSION FOR NOW.

    print ('Calculating total number of mapped reads for '+lib_name+'...')

    to_run=pairing(paired,lib_name,destination,modifier)

    #Generate sorted.bam index to be used in samtools idxstats.
    prog='samtools index '
    command=prog+destination+lib_name+modifier+to_run[1]+'_sorted.bam '
    sp.call(command,shell=True)
    
    #Generate stats of sorted.bam.
    prog='samtools idxstats '
    command=prog+destination+lib_name+modifier+to_run[1]+'_sorted.bam '+\
             '> '+\
             destination+lib_name+modifier+to_run[1]+'.stats '
    sp.call(command,shell=True)

    #Get total mapped reads.
    total=0.0
    a=open(destination+lib_name+modifier+to_run[1]+'.stats','r')
    for line in a:
        ar=read.tab(line)
        total += float(ar[2])
    a.close()
    if paired==True:
        total=total/2.0/1000000.0

    print ('... reads counted.')
    return total



def bam_to_bedgraph(lib_names,destination,modifier,ref_genome,paired,maxlength):
    #STOP using this and instead do bedgraph()
    #Custom code to generate bedgraph from paired end mapping
    #where each bed fragment is the entire fragment from 1 mapped end to its
    #paired end.
    #USE SORTED AND RMDUP BAM.

    for lib_name in lib_names:
        
        print ('Convert '+lib_name+modifier+'_rmdup.bam to bedgraph...')

        to_run=pairing(paired,lib_name,destination,modifier)

        chrolists={}
        #Use pre-generated list of chromosome names.
        b=open('store/BT2/genomes/'+ref_genome+'/'+ref_genome+'_chromosomes.txt','r')
        for line in b:
            chrolists[line.strip()]=[3000000000,1]
        b.close()

        #Get total mapped read counts for RPM normalization.
        total=total_mapped_reads(lib_name,destination,modifier,paired)

        #Convert _rmdup.bam to sam here.
        prog='samtools view '
        command=prog+' '+\
                 '-h '+\
                 destination+lib_name+modifier+to_run[1]+'_rmdup.bam '+\
                 '-o '+\
                 destination+lib_name+modifier+to_run[1]+'_rmdup.sam '
        sp.call(command,shell=True)   

        print ('Filling intervals...')

        #Fill in interval between mapped paired-ends.
        #SAM FILE HAS TO BE ARRANGED IN PAIRS OF LINES THAT
        #CORRESPOND TO THE 2 PAIRS IN PAIR-END READ.
        a=open(destination+lib_name+modifier+to_run[1]+'_rmdup.sam','r')
        b=open(destination+lib_name+modifier+'_intervals.txt','w')

        for i in range(68):
            next(a)
        for line in a:
            ar=read.tab(line)
            chro=str(ar[2])
            test=int(ar[8])
            if abs(test) != -test:
                start=int(ar[3])-1
                end=int(ar[7])+74
                if start < chrolists[chro][0]:
                    chrolists[chro][0] = start
                if end > chrolists[chro][1]:
                    chrolists[chro][1] = end
            b.write(chro+'\t'+str(start)+'\t'+str(end)+'\n')
        a.close()
        b.close()
                
        print ('Making trackline...')       
        #Make trackline.
        d=open(destination+lib_name+modifier+'.bg','w')
        if ref_genome == 'hg38':
            d.write('browser position chr6:33751017-33954459\n')
        elif ref_genome == 'mm10':
            d.write('browser position chr17:27250609-27728218\n')
        d.write('browser hide all\n'\
                +'track type=bedGraph name="'+lib_name+'" '\
                +'visibility=full color=255,0,0 altColor=0,10,255\n')
        
        #FOR NOW, BECAUSE NOT BOTHERING WITH STRANDED INFO, I'M DISPLAYING
        #EVERYTHING AS IF THEY ARE ALL +STRAND. MIGHT NEED TO DEAL WITH THIS
        #IN THE FUTURE!!!!!
        
        #Calculate score at each coordinate.
        for chrolist in chrolists.keys():            

            print ('Calculating scores at each coordinate for '+chrolist+'...')

            scores={}
            scores[chrolist]={}
            for j in range(chrolists[chrolist][0],chrolists[chrolist][1]+1):
                scores[chrolist][j]=0.0
            b=open(destination+lib_name+modifier+'_intervals.txt','r')
            for line in b:
                ar=read.tab(line)
                chro=str(ar[0])
                if chro == chrolist:
                    start=int(ar[1])
                    end=int(ar[2])
                    for i in range(start,end):
                        scores[chro][i] += 1.0
            b.close()

            print ('Deleting zero-scores for '+chrolist+'...')

            for i in scores[chrolist].keys():
                if scores[chrolist][i] == 0.0:
                    del scores[chrolist][i]
            
            print ('Merging interavals for '+chrolist+'...')

            #Merge intervals to get continuous bedgraph.
            prevcoor=1
            prevscore=0
            d.write(chrolist+'\t'+str(prevcoor)+'\t')            
            for coor in sorted( scores[chrolist].keys() ):
                score= scores[chrolist][coor]
                if ( score != prevscore ) or \
                   ( coor != (prevcoor+1) ):
                    d.write(str(prevcoor+1)+'\t')
                    reads=str( round(float(prevscore)/total,2) )
                    d.write(reads+'\n')
                    d.write(chrolist+'\t'+str(coor)+'\t')
                prevscore=score
                prevcoor=coor
            d.write(str(prevcoor+1)+'\t')
            reads=str( round(float(prevscore)/total,2) )
            d.write(reads+'\n')
            del scores
        d.close()

        #Delete this sam file to save space.
        prog='rm '
        command=prog+\
                 destination+lib_name+modifier+to_run[1]+'_rmdup.sam '
        sp.call(command,shell=True)



def bam_to_bed_for_Exo_CLIP(lib_names,destination,modifier,ref_genome,paired,maxlength):
    #Custom code to generate bed file from single end mapping.
    #Retains stranded info.
    #USE SORTED BAM.
    #UNSCALED FOR NOW, BUT CHANGED TO SCALED LATER.

    for lib_name in lib_names:
        
        print ('Convert '+lib_name+modifier+'_sorted.bam to stranded bed...')

        to_run=pairing(paired,lib_name,destination,modifier)

        #Convert _sorted.bam to sam here.
        prog='samtools view '
        command=prog+' '+\
                 '-h '+\
                 destination+lib_name+modifier+to_run[1]+'_sorted.bam '+\
                 '-o '+\
                 destination+lib_name+modifier+to_run[1]+'_sorted.sam '
        sp.call(command,shell=True)   

        print ('Making trackline...')        
        #Make trackline.
        b=open(destination+lib_name+modifier+'.bed','w')
        if ref_genome == 'hg38':
            b.write('browser position chrM:1-16569\n')
            header=457
        elif ref_genome == 'mm10':
            b.write('browser position chrM:1-16,299\n')
            header=68
        elif ref_genome == 'K12':
            b.write('browser position chr:10,001-35,000\n')
            header=3
        elif ref_genome == 'UTI89':
            b.write('browser position chr:10,001-35,000\n')
            header=4

        b.write('browser hide all\n'\
                +'track type=bed name="'+lib_name+'" '\
                +'visibility=2 colorByStrand="255,0,0 0,0,255"\n')

        #Use pre-generated list of chromosome names.
        #Note that these lack the tricky chromosomes like chrUn and * (unmapped).
        w=open('store/BT2/genomes/'+ref_genome+'/'+ref_genome+'_chromosomes.txt','r')
        chromos=[]
        for line in w:
            chro=line.strip()
            chromos += [chro,]
        w.close()

        a=open(destination+lib_name+modifier+to_run[1]+'_sorted.sam','r')
        for i in range(header):
            next(a)
        bedcount=0
        for line in a:
            ar=read.tab(line)
            chro=str(ar[2])
            if chro in chromos:
                bedcount+=1
                flag=int(ar[1])
                start=int(ar[3])-1
                end=start+len(ar[9])
                if flag==0:
                    strand='+'
                    score='1'
                elif flag==16:
                    strand='-'
                    score='-1'
                b.write(chro+'\t'+str(start)+'\t'+str(end)+'\t'+\
                        str(bedcount)+'\t'+score+'\t'+strand+'\n')
        a.close()
        b.close()
                
##        #Delete this sam file to save space.
##        prog='rm '
##        command=prog+\
##                 destination+lib_name+modifier+to_run[1]+'_sorted.sam '
##        sp.call(command,shell=True)
##
        print (lib_name+modifier+' bed made.')



def sorted_bam_to_sam(lib_names,destination,modifier,ref_genome,paired,read_end,origin):

    for lib_name in lib_names:

        to_run=pairing(paired,lib_name,destination,modifier,read_end,origin)

        print ('Convert '+lib_name+modifier+to_run[1]+'_sorted.bam to sorted.sam...')

        #Convert _sorted.bam to sam here.
        prog='samtools view '
        command=prog+' '+\
                 '-h '+\
                 destination+lib_name+modifier+to_run[1]+'_sorted.bam '+\
                 '-o '+\
                 destination+lib_name+modifier+to_run[1]+'_sorted.sam '
        
        sp.call(command,shell=True)   



        
def rmdup_via_UMI(lib_names,destination,modifier,ref_genome,paired,read_end,origin):
    #Custom code to remove duplicates via common UMIs and mapped coordinates.
    #Need to sort sam file via mapped coordinates first for input here.
    #Then remove duplicates via same mapped coordinates and UMI.
    #Requires paired-end mapping.
    #Only outputs the lines that match the 2 paired-end flags 99,83,147,163
    
    for lib_name in lib_names:

        to_run=pairing(paired,lib_name,destination,modifier,read_end,origin)

        print ('Removing PCR duplicates for '+lib_name+'...')

        if paired==True:
            flag1=99#For pair-end
            flag2=83#For pair-end
        aligns={}
        aligns[flag1]={}
        aligns[flag2]={}
            
        if ref_genome == 'pg435':
            header=3
        elif ref_genome == 'hg38':
            header=457

        total_count=0
        a=open(destination+lib_name+modifier+to_run[1]+'_sorted.sam','r')
        for i in range(header):
            next(a)        
        for line in a:
            ar=read.tab(line)
            flag=int(ar[1])
            if (flag==flag1) or (flag==flag2):
                total_count+=1
                chro = str(ar[2])
                start=int(ar[3])-1 
                leng=abs(int(ar[8]))
                name=str(ar[0])
                try:
                    aligns[flag][chro,start,leng].append(name)
                except KeyError:
                    aligns[flag][chro,start,leng]=[name,]
        a.close()

        #Actual rmdup
        non_dup=[]
        duplicate_count=0
        for flag in [flag1,flag2]:
            for key in aligns[flag].keys():
                if len( aligns[flag][key] ) ==1:
                    non_dup += [aligns[flag][key][0],]
                else:
                    UMIs_so_far=[]
                    for name in aligns[flag][key]:                        
                        UMI=name.split(':')[-1]
                        if UMI not in UMIs_so_far:
                            UMIs_so_far.append(UMI)
                            non_dup += [name,]
                        else:
                            duplicate_count+=1
            del aligns[flag]
            

        b=open(destination+lib_name+modifier+to_run[1]+'_sorted_rmdup.sam','w')
        #Write in header into rmdup.sam
        a=open(destination+lib_name+modifier+to_run[1]+'_sorted.sam','r')
        linecount=0
        for line in a:
            b.write(line.strip()+'\n')
            linecount+=1
            if linecount==header:
                break
        sams={}
        a=open(destination+lib_name+modifier+to_run[1]+'_sorted.sam','r')
        for i in range(header):
            next(a)
        for line in a:
            ar=read.tab(line)
            try:
                sams[ar[0]] += [line,]
            except KeyError:
                sams[ar[0]] = [line,]
        a.close()
        for name in non_dup:
            for line in sams[name]:
                b.write(line)
        del sams
        print ('total',total_count)
        print ('non_dup',len(non_dup))
        print ('percent_dup',float(len(non_dup))/float(total_count)*100.0)
            
        c=open(destination+lib_name+modifier+to_run[1]+'_rmdup-stats.txt','w')
        dup_percent=round(float(duplicate_count)/float(total_count)*100.0,1)
        non_dup_percent=100.0 - dup_percent
        c.write(lib_name+modifier+to_run[1]+'\n'+\
                str(total_count)+' mapped alignments\n'+\
                str(duplicate_count)+' duplicate alignments '+str(dup_percent)+'%\n'+\
                str(total_count-duplicate_count)+' non-duplicate alignments '+str(non_dup_percent)+'%')
                
            
        #Delete this sam file to save space.
        prog='rm '
        command=prog+\
                 destination+lib_name+modifier+to_run[1]+'_sorted.sam '
        sp.call(command,shell=True)
    


def five_prime_bg(lib_names,destination,modifier,ref_genome,paired,read_end,origin,UMIng):
    #Custom code to generate bedgraph from single end mapping
    #Only maps first nucleotide.
    #USE SORTED and clipped.

    from collections import defaultdict

    for lib_name in lib_names:
            
        to_run=pairing(paired,lib_name,destination,modifier,read_end,origin)

        if UMIng==True:
            samversion='_rmdup'
        else:
            samversion=''
            
        print ('Convert '+lib_name+modifier+to_run[1]+'_sorted_rmdup.sam to bedgraph...')

        print ('Making trackline...')       
        #Make trackline.
        b=open(destination+lib_name+modifier+to_run[1]+'_5prime.bg','w')
        if ref_genome == 'hg38':
            b.write('browser position chrM:1-16569\n')
            header=457
        elif ref_genome == 'mm10':
            b.write('browser position chrM:1-16,299\n')
            header=68
        elif ref_genome == 'K12':
            b.write('browser position chr:10,001-35,000\n')
            header=3
        elif ref_genome == 'UTI89':
            b.write('browser position chr:10,001-35,000\n')
            header=4
        elif ref_genome == '14028S':
            b.write('browser position chr:10,001-35,000\n')
            header=3
        elif ref_genome == 'pg482':
            b.write('browser position pg482:1-200\n')
            header=3
        elif ref_genome == 'pg435':
            b.write('browser position pg435:1-57\n')
            header=3
        elif ref_genome == 'pg144':
            b.write('browser position pg144:1-20\n')
            header=3
        elif ref_genome == 'pg494':
            b.write('browser position pg494:1-200\n')
            header=3
            
        b.write('browser hide all\n'\
                +'track type=bedGraph name="'+lib_name+to_run[1]+'_5prime" '\
                +'visibility=full color=255,0,0 altColor=0,10,255\n')

        #Calculate reads at each +strand start or -strand end position.
        print ('Calculating reads at each start position...')
        #Use pre-generated list of chromosome names.
        #Note that these lack the tricky chromosomes like chrUn and * (unmapped).


        if ref_genome=='hg38':
            chrMend=16569
        elif ref_genome=='mm10':
            chrMend=16299

        if paired == False:
            flag1=0
            flag2=16
        elif paired == True:
            flag1=99
            flag2=83
        
        w=open('store/BT2/genomes/'+ref_genome+'/'+ref_genome+'_chromosomes.txt','r')
        for line in w:
            chro=line.strip()
            beds_pos = defaultdict(int)
            beds_neg = defaultdict(int)
            
            a=open(destination+lib_name+modifier+to_run[1]+'_sorted'+samversion+'.sam','r')
            for i in range(header):
                next(a)
            for line in a:
                ar=read.tab(line)
                if chro == str(ar[2]):
                    start=int(ar[3])-1 
                    flag=int(ar[1])
                    if flag==flag1:                       
                        beds_pos[start] += 1                       
                    elif flag==flag2:
                        end=start+len(ar[9])
                        beds_neg[end] += 1                       
            a.close()
            
            #Write reads into .bg.
            for start in beds_pos.keys():
                b.write(chro+'\t'+str(start)+'\t'+str(start+1)+'\t'+ \
                        str(beds_pos[start])+'\n')
            for end in beds_neg.keys():
                
                #Have to do this for all circular genomes.
                if ( chro=='chrM' ) and ( end > chrMend ):
                    end=end-chrMend
                    
                b.write(chro+'\t'+str(end-1)+'\t'+str(end)+'\t'+ \
                        str(-beds_neg[end])+'\n')
            del beds_pos
            del beds_neg
        b.close()

        print ('.bg written.')


   
def bedgraph(lib_names,destination,modifier,ref_genome,paired,read_end,origin):
    #NOT REALLY USED BECAUSE GOT .bigwig MODULE.

    #This version is to use sorted_rmdup.sam (which ironically is no longer sorted because of my manipulation)
    # to make RPM.bedgraph.
    #Also converts sorted_rmdup.sam to sorted_rmdup.bam, which is required input for bedtools genomecov

    from collections import defaultdict

    for lib_name in lib_names:
        to_run=pairing(paired,lib_name,destination,modifier,read_end,origin)

        print ('Generate bedgraph for '+lib_name+to_run[1]+'...')

        #Get total mapped read count (all unique counts + one of each multi count).
        flag1=99#For paired-end
        flag2=83#For paired-end
        total=0
        b=open(destination+lib_name+modifier+to_run[1]+'_sorted_rmdup.sam','r')
        if ref_genome == 'hg38':
            header=457
        for i in range(header):
            next(b)        
        for line in b:
            ar=read.tab(line)
            flag=int(ar[1])
            if (flag==flag1) or (flag==flag2):
                total+=1
        b.close()        
        print (lib_name+to_run[1],'total',total)
        total=float(total)/1000000.0#Million mapped reads.

        #Convert sorted_rmdup.sam to rmdup.bam. Note that sorted_rmdup.sam is not really sorted at this point.
        prog='samtools view -S -b '        
        sp.call(prog+' '+destination+lib_name+modifier+to_run[1]+'_sorted_rmdup.sam > '+\
                destination+lib_name+modifier+to_run[1]+'_rmdup.bam ',\
                shell=True)

        #Sort rmdup.bam.
        prog='samtools sort '        
        sp.call(prog+' '+destination+lib_name+modifier+to_run[1]+'_rmdup.bam '+\
                destination+lib_name+modifier+to_run[1]+'_sorted_rmdup ',\
                shell=True)
        
        #Delete rmdup.bam to save space.
        prog='rm '        
        sp.call(prog+' '+destination+lib_name+modifier+to_run[1]+'_rmdup.bam ',\
                shell=True)

        #Generate bedgraphs.
        prog='bedtools genomecov -bg -split '        
        sp.call(prog+' -strand + -ibam '+destination+lib_name+modifier+to_run[1]+'_sorted.bam '+\
                '-g '+destination+'/genomes/chrom_sizes/'+ref_genome+'.chrom.sizes > '+\
                destination+lib_name+modifier+to_run[1]+'_plus.bedgraph ',\
                shell=True)
        sp.call(prog+' -strand - -ibam '+destination+lib_name+modifier+to_run[1]+'_sorted.bam '+\
                '-g '+destination+'/genomes/chrom_sizes/'+ref_genome+'.chrom.sizes > '+\
                destination+lib_name+modifier+to_run[1]+'_minus.bedgraph ',\
                shell=True)
        

        #Combine plus and minus bedgraphs and output as RPM
        c=open(destination+lib_name+modifier+to_run[1]+'_RPM.bedgraph','w')
        b=open(destination+lib_name+modifier+to_run[1]+'_plus.bedgraph','r')
        for line in b:
            ar=read.tab(line)
            rpm=round(float(ar[3])/total,3)
            c.write(ar[0]+'\t'+ar[1]+'\t'+ar[2]+'\t'+str(rpm)+'\n')
        b.close()
        b=open(destination+lib_name+modifier+to_run[1]+'_minus.bedgraph','r')
        for line in b:
            ar=read.tab(line)
            rpm=-round(float(ar[3])/total,3)
            c.write(ar[0]+'\t'+ar[1]+'\t'+ar[2]+'\t'+str(rpm)+'\n')
        b.close()
        c.close()

        #Remove plus and minus bedgraphs to save space.
        prog='rm '
        sp.call(prog+' '+destination+lib_name+modifier+to_run[1]+'_plus.bedgraph ',\
                shell=True)
        sp.call(prog+' '+destination+lib_name+modifier+to_run[1]+'_minus.bedgraph ',\
                shell=True)
        
        #gzip new .bg file.
        prog='gzip '
        sp.call(prog+' '+destination+lib_name+modifier+to_run[1]+'_RPM.bedgraph ',\
                shell=True)

        print (lib_name+to_run[1]+'_RPM.bedgraph normalized.')




def bam_to_bigwig(lib_names,destination,modifier,ref_genome,paired,read_end,threads,egz,origin):

    #This version is to use sorted_rmdup.sam (which ironically is no longer sorted because of my manipulation)
    #Also converts sorted_rmdup.sam to sorted_rmdup.bam, which is required input for bedtools genomecov

    from collections import defaultdict

    for lib_name in lib_names:
        to_run=pairing(paired,lib_name,destination,modifier,read_end,origin)

        print ('Generate bigwig for '+lib_name+to_run[1]+'...')

##        ###FOR RMDUP###
##        #Convert sorted_rmdup.sam to rmdup.bam. Note that rmdup.sam is not really sorted at this point.
##        prog='samtools view -S -b '        
##        sp.call(prog+' '+destination+lib_name+modifier+to_run[1]+'_sorted_rmdup.sam > '+\
##                destination+lib_name+modifier+to_run[1]+'_rmdup.bam ',\
##                shell=True)
##        #Sort rmdup.bam.
##        prog='samtools sort '        
##        sp.call(prog+' '+destination+lib_name+modifier+to_run[1]+'_rmdup.bam '+\
##                destination+lib_name+modifier+to_run[1]+'_sorted_rmdup ',\
##                shell=True)
##        #Delete rmdup.bam to save space.
##        prog='rm '        
##        sp.call(prog+' '+destination+lib_name+modifier+to_run[1]+'_rmdup.bam ',\
##                shell=True)
##        ###

        #Index bam file so it can be used for bamCoverage
        prog='samtools index '+destination+lib_name+modifier+to_run[1]+'_sorted.bam '
        sp.call(prog,shell=True)

        #Actual bamCoverage
        prog='bamCoverage '
        prog += ('-b '+destination+lib_name+modifier+to_run[1]+'_sorted.bam '+\
                 '--effectiveGenomeSize '+egz+\
                 '--normalizeUsing CPM '+\
                 '-p 30 '+\
                 '-e '+\
                 '-o '+destination+lib_name+modifier+to_run[1]+'.bw ')

        sp.call(prog,shell=True)
        
        print (lib_name+to_run[1]+'_bigwig generated.')
