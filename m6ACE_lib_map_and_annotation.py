#!/usr/bin/python
#Pipeline for mapping and annotating m6ACE-seq libraries.
#Spike-in mapping is done with a separate pipeline.
# for fastq processing.
#NOTE STAR uses ~80Gb for 80M input reads.

import read
import fasty
import browser_from_star as bfs

#INPUT Parameters here.

lib_names=[]
for i in range(907,911):
    lib_names += ['LSG'+str(i),]
    
ORIGIN='store/fastqgz/'
DESTINATION='store/STARdir/'
REF_GENOME='hg38_pr8'
THREADS=30
PAIRED=True
POLYA=True
to_UMI=True
READ_END=['_R1','_R2']#If want to do SE on R2, just flip the order in this list, and adapters.
RA3='TGGAATTCTCGGGTGCCAAGG'#For m6ACE
RA5='GATCGTCGGACTGTAGAACTCTGAAC'#For m6ACE

QUALITY=[20,False]#[1] denotes whether nextseq (True) or not nextseq run (False).
ATAIL='A{100}'
TTAIL='T{100}'
CLIP_LENGTH=20#20 for RNA.

READWITHUMI='_R1'
READWITHOUTUMI='_R2'
UMIPOSITION=[1,8]
MODIFIER='_clipped'




#START PIPELINE HERE:
##
######Unzip fastq.gz.
fasty.unzip_fastq(lib_names,origin=ORIGIN,read_end=READ_END,paired=PAIRED)


#Clip fastq.
fasty.fastq_clipping(lib_names,origin=ORIGIN,clip_length=CLIP_LENGTH,\
                     adapter1=RA3,adapter2=RA5,paired=PAIRED,read_end=READ_END,\
                     polya=POLYA,atail=ATAIL,ttail=TTAIL,qual=QUALITY,threads=THREADS)

##Trim off UMI and add it to the sequence name in the fastq.
if to_UMI==True:
    fasty.fastq_UMI(lib_names,origin=ORIGIN,paired=PAIRED,readwithUMI=READWITHUMI,\
                    readwithoutUMI=READWITHOUTUMI,UMIposition=UMIPOSITION,modifier=MODIFIER,read_end=READ_END)
    MODIFIER=MODIFIER+'_identified'

##
####
if (to_UMI==True) and ('_identified' not in MODIFIER):
    MODIFIER=MODIFIER+'_identified'#ALWAYS HAVE THIS ON
    
##########Run STAR.
bfs.run_star(lib_names,origin=ORIGIN,destination=DESTINATION,\
                 modifier=MODIFIER,ref_genome=REF_GENOME,\
                 paired=PAIRED,threads=THREADS,read_end=READ_END)

#Converts sorted.bam to sorted.sam
bfs.sorted_bam_to_sam(lib_names,destination=DESTINATION,\
                      modifier=MODIFIER,ref_genome=REF_GENOME,\
                      paired=PAIRED,read_end=READ_END,origin=ORIGIN)

#Remove duplicates via common mapped coordinates and UMIs.
if to_UMI==True:
    bfs.rmdup_via_UMI(lib_names,destination=DESTINATION,\
                      modifier=MODIFIER,ref_genome=REF_GENOME,\
                      paired=PAIRED,read_end=READ_END,origin=ORIGIN)



##Convert sorted_rmdup.sam (or sorted.sam but mainly rmdup version)
#to stranded bedgraph. Used for PE libs from m6ACE.
#Only plot first nucleotide.
bfs.five_prime_bg(lib_names,destination=DESTINATION,\
                  modifier=MODIFIER,ref_genome=REF_GENOME,\
                  paired=PAIRED,read_end=READ_END,origin=ORIGIN,UMIng=to_UMI)


#Normalized .bg file by RPM and gzip it.
bfs.normalize_bg(lib_names,destination=DESTINATION,\
                 modifier=MODIFIER,ref_genome=REF_GENOME,\
                 paired=PAIRED,read_end=READ_END,origin=ORIGIN)
##
##






