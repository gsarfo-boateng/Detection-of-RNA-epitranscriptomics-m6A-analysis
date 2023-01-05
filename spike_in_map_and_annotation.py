#!/usr/bin/python
#Pipeline for mapping and annotating m6ACE spike-in. Based on bowtie2.

import subprocess as sp
import read
import fasty
import bowtie2_to_browser as btb

#INPUT Library names here.

lib_names=[]
for i in range(903,911):
    lib_names+=['LSG'+str(i),]


#For pg435 mapping, paired=true, polyA=true, to_umi=true
PAIRED=True
POLYA=True
to_UMI=True
CLIP_LENGTH=20
to_trim=False###

#PARAMETERS
ORIGIN='store/fastqgz/'
DESTINATION='store/BT2/'
READ_END=['_R1','_R2']#If want to do SE on R2, just flip the order in this list, and adapters.
MODIFIER='_clipped'
FIRST=1#For trimming
LAST=75#For trimming
RA3='TGGAATTCTCGGGTGCCAAGG'
RA5='GATCGTCGGACTGTAGAACTCTGAAC'
QUALITY=[20,False]#[1] denotes whether nextseq (True) or not nextseq run (False).
ATAIL='A{100}'
TTAIL='T{100}'
READWITHUMI='_R1'
READWITHOUTUMI='_R2'
UMIPOSITION=[1,8]

#START PIPELINE HERE:
####FASTQ SECTION.
####fasty.unzip_fastq, fasty.fastq_clipping and fasty.fastq_UMI all performed in m6ACE_lib_map_and_annotation.py

#BOWTIE2 MAPPING SECTION.
REF_GENOME='pg435'
THREADS=30
MISMATCH=1
if to_UMI==True:
    MODIFIER=MODIFIER+'_identified'
if to_trim==True:
    MODIFIER=MODIFIER+'_trimmed'

#######Run bowtie2.
btb.run_bowtie2(lib_names,destination=DESTINATION,\
                modifier=MODIFIER,ref_genome=REF_GENOME,\
                mismatch=MISMATCH,\
                paired=PAIRED,threads=THREADS,read_end=READ_END,origin=ORIGIN)

##Convert bowtie2 output .sam file into sorted.bam file.
btb.sam_to_sorted_bam(lib_names,destination=DESTINATION,\
                      modifier=MODIFIER,\
                      ref_genome=REF_GENOME,\
                      paired=PAIRED,read_end=READ_END,origin=ORIGIN)


#Converts sorted.bam to sorted.sam
btb.sorted_bam_to_sam(lib_names,destination=DESTINATION,\
                      modifier=MODIFIER,ref_genome=REF_GENOME,\
                      paired=PAIRED,read_end=READ_END,origin=ORIGIN)

#Remove duplicates via common mapped coordinates and UMIs.
if to_UMI==True:
    btb.rmdup_via_UMI(lib_names,destination=DESTINATION,\
                      modifier=MODIFIER,ref_genome=REF_GENOME,\
                      paired=PAIRED,read_end=READ_END,origin=ORIGIN)

##Convert PE sorted.sam to stranded bedgraph.
#Only plot first nucleotide.
btb.five_prime_bg(lib_names,destination=DESTINATION,\
                  modifier=MODIFIER,ref_genome=REF_GENOME,\
                  paired=PAIRED,read_end=READ_END,origin=ORIGIN,UMIng=to_UMI)







