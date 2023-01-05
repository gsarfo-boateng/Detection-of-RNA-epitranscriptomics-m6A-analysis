#!/usr/bin/python
#Scripts for STAR mapping, then converts to bedgraph or browser file.
#STAR can take both PE and SE but PE is preferred.
#See https://groups.google.com/forum/#!topic/rna-star/GuUxYI6RHJw for multiple-mapping issues.

#hg38_transcriptome.gtf is gencode.v28.chr_patch_hapl_scaff.annotation.gtf downloaded from https://www.gencodegenes.org/releases/current.html
#mm10_transcriptome.gtf is gencode.vM24.chr_patch_hapl_scaff.annotation.gtf downloaded from https://www.gencodegenes.org/mouse/
#Drosophila_melanogaster.BDGP6.94_chr_formatted.gtf is from Brian from Eric Lai Lab.

import subprocess as sp
import read
import complimentary


def pairing(paired,lib_name,origin,modifier,read_end):
    #For paired end mapping.
    if paired==True:
        to_run=[origin+lib_name+read_end[0]+'_PE'+modifier+'.fastq '+\
                origin+lib_name+read_end[1]+'_PE'+modifier+'.fastq ',\
                '_PE']
    #For single end mapping.
    elif paired==False:
        to_run=[origin+lib_name+read_end[0]+modifier+'.fastq ',\
                '_SE']
    return to_run

def sam_header(ref_genome,destination):
    #Denote number of sam header lines present.
    a=open(destination+'genomes/headers.txt','r')
    for line in a:
        ar=read.tab(line)
        if ar[0]==ref_genome:
            b=int(ar[1])
            break
    return b


def run_star(lib_names,origin,destination,\
                modifier,ref_genome,\
                paired,threads,read_end):
    #Running STAR.

    for lib_name in lib_names:
    
        to_run=pairing(paired,lib_name,origin,modifier,read_end)

        #Make subdirectory for STAR output.
        prog='mkdir '
        command=prog+\
                 destination+lib_name+to_run[1]
        sp.call(command,shell=True)

        if paired==True:
            #For PE.
            print ('Running STAR on '+\
                  lib_name+read_end[0]+modifier+'.fastq and '+\
                  lib_name+read_end[1]+modifier+'.fastq...')
        elif paired==False:
            #For SE.
            print ('Running STAR on '+\
                  lib_name+read_end[0]+modifier+'.fastq...')

        #Run STAR here.
##        prog='bin/STAR-2.6.0c/bin/Linux_x86_64_static/STAR '
        prog='STAR '
        #Immediately outputs sorted.bam file
        #For multiple-mappers with multiple top-scoring alignments, randomly chooses 1 top scorer and reports only that.
        command=prog+\
                 '--runThreadN '+str(threads)+' '+\
                 '--outSAMtype BAM SortedByCoordinate '+\
                 '--outSAMattributes All '+\
                 '--outMultimapperOrder Random --outSAMmultNmax 1 '+\
                 '--outFileNamePrefix '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+' '+\
                 '--genomeDir '+destination+'index_files/'+ref_genome+' '+\
                 '--readFilesIn '+to_run[0]

        #Only add following line to command for genomes with exons:
        if ref_genome in ['hg38','mm10','hg38_pr8','mm10_mhvA59','dm6']:
            command+='--sjdbGTFfile '+destination+'genomes/'+ref_genome+'/'+ref_genome+'_transcriptome.gtf '

##        #Only add following line to command for genomes with IAV infection:
##        #NO LONGER NECESSARY
##        if ref_genome in ['hg38_pr8',]:
##            command+='--outReadsUnmapped Fastx '

        #Following is for mapping of cap-snatched IAV sequences because
        #when mapping to pr8, the default RAM set to sort BAM is tied to pr8 genome size,
        #so have to add this option to increase the RAM available:
        if ref_genome in ['pr8',]:
            command+='--limitBAMsortRAM 2000000000 '
            
        sp.call(command,shell=True)
        print (lib_name+modifier,'...STAR done.')



def sorted_bam_to_sam(lib_names,destination,modifier,ref_genome,paired,read_end,origin):

    for lib_name in lib_names:

        to_run=pairing(paired,lib_name,origin,modifier,read_end)

        print ('Convert '+lib_name+to_run[1]+'_sorted.bam to sorted.sam...')

        #Convert _sorted.bam to sam here.
        prog='samtools view '
        command=prog+' '+\
                 '-h '+\
                 destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+\
                 'Aligned.sortedByCoord.out.bam '+\
                 '-o '+\
                 destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+\
                 '_sorted.sam '
        sp.call(command,shell=True)   


        
def rmdup_via_UMI(lib_names,destination,modifier,ref_genome,paired,read_end,origin):
    #Custom code to remove duplicates via common UMIs and mapped coordinates.
    #Need to sort sam file via mapped coordinates first for input here.
    #Then remove duplicates via same mapped coordinates and UMI.
    #Requires paired-end mapping.
    #Outputs the lines that match the 2 paired-end flags 99,83, 
    #as well as the corresponding flags 147 and 163 (the 2nd of each pair).
    
    for lib_name in lib_names:

        to_run=pairing(paired,lib_name,origin,modifier,read_end)

        print ('Removing PCR duplicates for '+lib_name+'...')

        flag1=99#For paired-end
        flag2=83#For paired-end
        aligns={}
        aligns[flag1]={}
        aligns[flag2]={}

        header=sam_header(ref_genome,destination)
        
        total_count=0
        a=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted.sam','r')
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
            

        b=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted_rmdup.sam','w')
        #Write in header into rmdup.sam
        a=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted.sam','r')
        linecount=0
        for line in a:
            b.write(line.strip()+'\n')
            linecount+=1
            if linecount==header:
                break
        sams={}
        a=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted.sam','r')
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
        print ('non_dup',len(non_dup))
            
        c=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_rmdup-stats.txt','w')
        dup_percent=round(float(duplicate_count)/float(total_count)*100.0,1)
        non_dup_percent=100.0 - dup_percent
        c.write(lib_name+modifier+to_run[1]+'\n'+\
                str(total_count)+' mapped alignments\n'+\
                str(duplicate_count)+' duplicate alignments '+str(dup_percent)+'%\n'+\
                str(total_count-duplicate_count)+' non-duplicate alignments '+str(non_dup_percent)+'%')
                
            
        #Delete this sam file to save space.
        prog='rm '
        command=prog+\
                 destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted.sam '
        sp.call(command,shell=True)
    


def rm_multimappers(lib_names,destination,modifier,ref_genome,paired,read_end,origin):
    #Custom code to remove duplicates via common UMIs and mapped coordinates.
    #Need to sort sam file via mapped coordinates first for input here.
    #Then remove duplicates via same mapped coordinates and UMI.
    #Requires paired-end mapping.
    #Outputs the lines that match the 2 paired-end flags 99,83, 
    #as well as the corresponding flags 147 and 163 (the 2nd of each pair).
    
    for lib_name in lib_names:

        to_run=pairing(paired,lib_name,origin,modifier,read_end)

        print ('Removing multimappers for '+lib_name+'...')

        flags=[99,83,147,163]#For paired-end

        header=sam_header(ref_genome,destination)

        total_count=0
        unique_count=0
        a=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted_rmdup.sam','r')
        uniques=[]
        for i in range(header):
            next(a)        
        for line in a:
            ar=read.tab(line)
            flag=int(ar[1])
            if flag in flags:
                total_count+=1
                mapped= int(str(ar[11])[5::])
                if mapped==1:
                    uniques+=[line,]
                    unique_count+=1
        a.close()
        total_count=total_count/2#Since I counted both forward and reverse reads.
        unique_count=unique_count/2#Since I counted both forward and reverse reads.

        a=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted_rmdup.sam','r')
        #Write in header into rmdup.sam
        b=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted_rmdup_unique.sam','w')
        linecount=0
        for line in a:
            b.write(line.strip()+'\n')
            linecount+=1
            if linecount==header:
                break
        for line in uniques:
            b.write(line)
            
        c=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_rm_multimapper-stats.txt','w')
        c.write(lib_name+modifier+to_run[1]+'\n'+\
                str(int(total_count))+' mapped\n'+\
                str(int(total_count-unique_count))+' multimappers '+str(round(float(total_count-unique_count)/float(total_count)*100.0,1))+'%\n'+\
                str(int(unique_count))+' unique mappers '+str(round(float(unique_count)/float(total_count)*100.0,1))+'%')
                
            
        #Delete this sam file to save space.
        prog='rm '
        command=prog+\
                 destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted_rmdup.sam '
        sp.call(command,shell=True)

        
    
    
def five_prime_bg(lib_names,destination,modifier,ref_genome,paired,read_end,origin,UMIng):
    #Custom code to generate bedgraph of only 5'end of R1.
    #Only maps first nucleotide.
    #USE SORTED and clipped.

    from collections import defaultdict

    for lib_name in lib_names:
        to_run=pairing(paired,lib_name,origin,modifier,read_end)

        if UMIng==True:
            samversion='_rmdup'
        else:
            samversion=''

        print ('Convert '+lib_name+to_run[1]+'_sorted'+samversion+'.sam to bedgraph...')

        b=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_5prime.bg','w')

        print ('Making trackline...')       
        #Make trackline.
        header=sam_header(ref_genome,destination)
        if ref_genome == 'hg38':
            b.write('browser position chr11:65500306-65500516\n')
        elif ref_genome == 'mm10':
            b.write('browser position chrM:1-16,299\n')
        elif ref_genome == 'mm10_mhvA59':
            b.write('browser position chrM:1-16,299\n')
        elif ref_genome == 'hg38_pr8':
            b.write('browser position chr11:65500306-65500516\n')
        elif ref_genome == 'dm6':
            b.write('browser position chr2L:100-1000\n')
            
        b.write('browser hide all\n'\
                +'track type=bedGraph name="'+lib_name+to_run[1]+'_5prime" '\
                +'visibility=full color=255,0,0 altColor=0,10,255\n')

        #Calculate reads at each +strand start or -strand end position.
        print ('Calculating reads at each start position...')


        if ref_genome=='hg38':
            chrMend=16569
        elif ref_genome=='hg38_pr8':
            chrMend=16569
        elif ref_genome=='mm10':
            chrMend=16299
        elif ref_genome=='mm10_mhvA59':
            chrMend=16299
        elif ref_genome=='dm6':
            chrMend=19524

        if paired == False:
            flag1=0
            flag2=16
        elif paired == True:
            flag1=99
            flag2=83

        beds_pos={}
        beds_neg={}
        #Use pre-generated list of chromosome names.
        w=open(destination+'genomes/'+ref_genome+'/'+ref_genome+'_chromosomes.txt','r')
        for line in w:
            chro=line.strip()
            beds_pos[chro] = defaultdict(int)
            beds_neg[chro] = defaultdict(int)
        w.close()
        
        a=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted'+samversion+'.sam','r')
        for i in range(header):
            next(a)
        for line in a:
            ar=read.tab(line)
            chro = str(ar[2])
            start=int(ar[3])-1 
            flag=int(ar[1])
            leng=len(ar[9])
            if flag==flag1:                       
                beds_pos[chro][start] += 1
            elif flag==flag2:
                #Following are to factor in intron distances in calculating actual start
                #site for a negative strand read that spans exon junction(s).
                JI=ar[18].split(',')
                if JI[1]=='-1':
                    end=start+leng
                else:
                    del JI[0]
                    intron=0
                    add=False
                    for ji in JI:
                        if add==False:
                            intron=intron-(int(ji)-1)
                            add=True
                        elif add==True:
                            intron=intron+int(ji)
                            add=False
                    end=start+leng+intron                
                beds_neg[chro][end] += 1
        a.close()

        #Remove any sites that exceed maximum chrom_size of each chr.
        chrosizes={}
        x=open('store/STARdir/genomes/chrom_sizes/'+ref_genome+'.chrom.sizes','r')
        for line in x:
            ar=read.tab(line)
            chrosizes[ar[0]]=int(ar[1])
        x.close()
        chrosizes['chrM']=100000#Doing this because chrM is circular chromosome so special case.
        
        w=open('store/STARdir/genomes/'+ref_genome+'/'+ref_genome+'_chromosomes.txt','r')
        for line in w:
            chro=line.strip()
            maxcoor=chrosizes[chro]
            for start in sorted(beds_pos[chro].keys()):
                if start>=maxcoor:
                    del beds_pos[chro][start]
            for end in sorted(beds_neg[chro].keys()):
                if end>=maxcoor:
                    del beds_neg[chro][end]
        del chrosizes
        
        #Write reads into .bg.
        w=open('store/STARdir/genomes/'+ref_genome+'/'+ref_genome+'_chromosomes.txt','r')
        for line in w:
            chro=line.strip()
                
            for start in sorted(beds_pos[chro].keys()):
                score= beds_pos[chro][start]
                b.write(chro+'\t'+str(start)+'\t'+str(start+1)+'\t'+ \
                        str(score)+'\n')
            for end in sorted(beds_neg[chro].keys()):
                #Have to do this for all circular genomes.
                if ( chro=='chrM' ) and ( end > chrMend ):
                    end=end-chrMend
                score= -beds_neg[chro][end]
                b.write(chro+'\t'+str(end-1)+'\t'+str(end)+'\t'+ \
                        str(score)+'\n')
            del beds_pos[chro]
            del beds_neg[chro]                
        del beds_pos
        del beds_neg
        b.close()
        print ('.bg written.')


def normalize_bg(lib_names,destination,modifier,ref_genome,paired,read_end,origin):

    #This version is only for 5prime.bg files from m6ACE.

    from collections import defaultdict

    for lib_name in lib_names:
        to_run=pairing(paired,lib_name,origin,modifier,read_end)

        print ('Normalize '+lib_name+to_run[1]+'_5prime.bg and gzip it...')

        #Get total mapped read count (all unique counts + one of each multi count).
        c=open('store/STARdir/'+lib_name+'_PE/'+lib_name+'_PE_rmdup-stats.txt','r')
        for j in range(3):
            next(c)
        for line in c:
            ar=read.space(line)
            total=float(ar[0])
        print (lib_name+to_run[1],'total',total)
        total=total/1e6#for RPM.
        c.close()

        #Output normalized .bg in RPM
        c=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_5prime_RPM.bedgraph','w')
        b=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_5prime.bg','r')
        linecount=1
        for line in b:
            if linecount>3:
                ar=read.tab(line)
                rpm=round(float(ar[3])/total,3)
                c.write(ar[0]+'\t'+ar[1]+'\t'+ar[2]+'\t'+str(rpm)+'\n')
            else:
                linecount+=1
                c.write(line)
        b.close()
        c.close()

        #gzip new .bg file.
        prog='gzip '
        print ('Zipping',lib_name+to_run[1]+'_5prime_RPM.bg','...')
        
        sp.call(prog+' '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_5prime_RPM.bedgraph ',\
                shell=True)

        print (lib_name+to_run[1]+'_5prime.bg normalized.')





def bedgraph(lib_names,destination,modifier,ref_genome,paired,read_end,origin):

    #CURRENTLY NOT WORKING. switched to using alternative method via bamCoverage tool.
    #This version is to use sorted_rmdup.sam (which ironically is no longer sorted because of my manipulation)
    # to make RPM.bedgraph, currently specialized for eCLIP but doesn't really output nice looking bedgraph at all.
    #Also converts sorted_rmdup.sam to sorted_rmdup.bam, which is required input for bedtools genomecov

    from collections import defaultdict

    for lib_name in lib_names:
        to_run=pairing(paired,lib_name,origin,modifier,read_end)

        print ('Generate bedgraph for '+lib_name+to_run[1]+'...')

        #Get total mapped read count (all unique counts + one of each multi count).
        flag1=99#For paired-end
        flag2=83#For paired-end
        total=0
        b=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted_rmdup.sam','r')

        header=sam_header(ref_genome,destination)            
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
        sp.call(prog+' '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted_rmdup.sam > '+\
                destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_rmdup.bam ',\
                shell=True)

        #Sort rmdup.bam.
        prog='samtools sort '        
        sp.call(prog+' '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_rmdup.bam '+\
                destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted_rmdup ',\
                shell=True)
        
        #Delete rmdup.bam to save space.
        prog='rm '        
        sp.call(prog+' '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_rmdup.bam ',\
                shell=True)

        #Generate bedgraphs.
        prog='bedtools genomecov -bg -split '        
        sp.call(prog+' -strand + -ibam '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted_rmdup.bam '+\
                '-g '+destination+'/genomes/chrom_sizes/'+ref_genome+'.chrom.sizes > '+\
                destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_plus.bedgraph ',\
                shell=True)
        sp.call(prog+' -strand - -ibam '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted_rmdup.bam '+\
                '-g '+destination+'/genomes/chrom_sizes/'+ref_genome+'.chrom.sizes > '+\
                destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_minus.bedgraph ',\
                shell=True)
        

        #Combine plus and minus bedgraphs and output as RPM
        c=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_RPM.bedgraph','w')
        b=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_plus.bedgraph','r')
        for line in b:
            ar=read.tab(line)
            rpm=round(float(ar[3])/total,3)
            c.write(ar[0]+'\t'+ar[1]+'\t'+ar[2]+'\t'+str(rpm)+'\n')
        b.close()
        b=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_minus.bedgraph','r')
        for line in b:
            ar=read.tab(line)
            rpm=-round(float(ar[3])/total,3)
            c.write(ar[0]+'\t'+ar[1]+'\t'+ar[2]+'\t'+str(rpm)+'\n')
        b.close()
        c.close()

        #Remove plus and minus bedgraphs to save space.
        prog='rm '
        sp.call(prog+' '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_plus.bedgraph ',\
                shell=True)
        sp.call(prog+' '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_minus.bedgraph ',\
                shell=True)
        
        #gzip new .bg file.
        prog='gzip '
        sp.call(prog+' '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_RPM.bedgraph ',\
                shell=True)

        print (lib_name+to_run[1]+'_RPM.bedgraph normalized.')



def bam_to_bigwig(lib_names,destination,modifier,ref_genome,paired,read_end,threads,egz,skiplength,origin,kit):
    #Use bamCoverage to generate stranded bigwig from paired end mapping
    #where each fragment is the entire fragment from 1 mapped end to its
    #paired end.
    #Input is STAR output.
    #Currently does not use rmdup or UMIs (mainly for Scriptseq V2 or other kinds of non-UMI RNA-seq).
    #Currently will bigwig over the whole intron (pseudo-fix is filter out intervals bigger than 300bp)
    #Also generates stranded bedgraph for other applications.

    for lib_name in lib_names:
        
        print ('Convert '+lib_name+modifier+' sorted bam to stranded bigwig...')

        to_run=pairing(paired,lib_name,origin,modifier,read_end)

        #Index bam file so it can be used for bamCoverage
        prog='samtools index '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'Aligned.sortedByCoord.out.bam '
        sp.call(prog,shell=True)

        #Actual bamCoverage to generate .bw
        #OUTPUT OF THIS IS MEANT FOR IGV ONLY.
        prog='bamCoverage '
        prog += ('-b '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'Aligned.sortedByCoord.out.bam '+\
                 '--effectiveGenomeSize '+egz+\
                 '--normalizeUsing CPM '+\
                 '-p 30 '+\
                 '-e '+\
                 '--maxFragmentLength '+skiplength)
        if kit==1:
            plus='--samFlagInclude 99 --scaleFactor 1.0 -o '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_plus.bw '
            minus='--samFlagInclude 83 --scaleFactor -1.0 -o '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_minus.bw '
        elif kit==2:
            plus='--samFlagInclude 83 --scaleFactor 1.0 -o '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_plus.bw '
            minus='--samFlagInclude 99 --scaleFactor -1.0 -o '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_minus.bw '
        sp.call(prog+plus,shell=True)
        sp.call(prog+minus,shell=True)

        #Actual bamCoverage (pre-normalized via RPM) to generate .bedgraph
        #Decided not to use --normalizeUsing CPM option of bamCoverage and instead use my own way to normalize via RPM for consistency
        #with eclip_sam_to_bamCoverage()
        #OUTPUT OF THIS IS MEANT FOR bedgraph_to_DESeq2.py
        prog='bamCoverage '
        prog += ('-b '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'Aligned.sortedByCoord.out.bam '+\
                 '--effectiveGenomeSize '+egz+\
                 '-p 30 '+\
                 '-e '+\
                 '-of bedgraph '+\
                 '--maxFragmentLength '+skiplength)
        if kit==1:
            plus='--samFlagInclude 99 -o '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_plus.temp '
            minus='--samFlagInclude 83 -o '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_minus.temp '
        elif kit==2:
            plus='--samFlagInclude 83 -o '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_plus.temp '
            minus='--samFlagInclude 99 -o '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_minus.temp '
        sp.call(prog+plus,shell=True)
        sp.call(prog+minus,shell=True)

        print (lib_name+modifier+' stranded bigwig and stranded temp done.')                 

                    
def eclip_sam_to_bamCoverage(lib_names,destination,modifier,ref_genome,paired,read_end,threads,egz,skiplength,origin,kit):
    #Use bamCoverage to generate stranded bedgraph from paired end mapping
    #where each fragment is the entire fragment from 1 mapped end to its
    #paired end.
    #Input is rmdup_sorted_unique.sam
    #Really the same as bam_to_bigwig just that this is adapted ONLY for eCLIP sorted_rmdup_unique.sam, is not normalized using bamCoverage but is normalized using my method.


    for lib_name in lib_names:
        
        print ('Convert '+lib_name+modifier+' sorted bam to stranded bigwig...')

        to_run=pairing(paired,lib_name,origin,modifier,read_end)


        #Convert sorted_rmdup_unique.sam to rmdup_unique.bam. Note that sorted_rmdup.sam is not really sorted at this point.
        prog='samtools view -S -b '        
        sp.call(prog+' '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted_rmdup_unique.sam > '+\
                destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_rmdup_unique.bam ',\
                shell=True)

        #Sort rmdup_unique.bam.
        prog='samtools sort '        
        sp.call(prog+' '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_rmdup_unique.bam '+\
                destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted_rmdup_unique ',\
                shell=True)
        
        #Delete rmdup_unique.bam to save space.
        prog='rm '        
        sp.call(prog+' '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_rmdup_unique.bam ',\
                shell=True)

        #Index bam file so it can be used for bamCoverage
        prog='samtools index '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted_rmdup_unique.bam '
        sp.call(prog,shell=True)

        
        #Actual bamCoverage (pre-normalized via RPM).
        prog='bamCoverage '
        prog += ('-b '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_sorted_rmdup_unique.bam '+\
                 '--effectiveGenomeSize '+egz+\
                 '-p 30 '+\
                 '-e '+\
                 '-of bedgraph '+\
                 '--maxFragmentLength '+skiplength)
        if kit==1:
            plus='--samFlagInclude 99 -o '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_plus.temp '
            minus='--samFlagInclude 83 -o '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_minus.temp '
        elif kit==2:
            plus='--samFlagInclude 83 -o '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_plus.temp '
            minus='--samFlagInclude 99 -o '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_minus.temp '
        sp.call(prog+plus,shell=True)
        sp.call(prog+minus,shell=True)

        #Get total mapped read count (ONLY unique counts).
        c=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_rm_multimapper-stats.txt','r')
        for i in range(3):
            next(c)
        for line in c:
            ar=read.space(line)
            total=float(ar[0])/1000000.0#Million mapped reads.
        print(lib_name+to_run[1],'total',total)

        #Normalize bedgraphs via total mapped reads i.e. RPM.
        c=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_plus.temp','r')
        d=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_plus.bedgraph','w')
        for line in c:
            ar=read.tab(line)
            if str(ar[3]) != '0':
                d.write(ar[0]+'\t'+ar[1]+'\t'+ar[2]+'\t'+str(round(float(ar[3])/total,3))+'\n')
        c=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_minus.temp','r')
        d=open(destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_minus.bedgraph','w')
        for line in c:
            ar=read.tab(line)
            if str(ar[3]) != '0':
                d.write(ar[0]+'\t'+ar[1]+'\t'+ar[2]+'\t'+str(-round(float(ar[3])/total,3))+'\n')

        #gzip new .bedgraph file.
        prog='gzip '
        for strand in ['plus','minus']:        
            sp.call(prog+' '+destination+lib_name+to_run[1]+'/'+lib_name+to_run[1]+'_'+strand+'.bedgraph ',\
                    shell=True)

        print(lib_name+' bedgraph.gz generated.')
