#motif_generation.py

#Universal module for generating motif sequence if given
#ref_genome
#chromosome
#strand
#start (0-based) or end (1-based).
#coors is dictionary in format coors[chr][strand][start or end]
#Currently only for single-base sites (e.g. from m6ACE start sites)

import complimentary

def fasta(ref_genome,upstream,seqlen,coors,filename):
    #This generates fasta file of all the desired sites being input in an indiscriminate manner.
    
    sequences=[]
    for chro in coors.keys():

        #Upload genomic seqeunce for each chromosome.         
        z=open('store/STARdir/genomes/genomes_parsed_by_chr/'+ref_genome+'/'+chro+'.fa','r')
        for line in z:
            genome=line.strip()
        z.close()

        #Go through list of coordinates for matching sites and output motif sequence
        for strand in coors[chro].keys():
            for coor in coors[chro][strand]:
                if strand=='+':
                    sequences.append( str.upper(genome[coor-upstream:coor+seqlen]) )
                elif strand=='-':
                    sequences.append( complimentary.seq( str.upper(genome[coor-seqlen:coor+upstream]) ) ) 
        del genome

    #Output fasta file of motif sequences
    a=open('store/motif/'+filename+'.fa','w')
    for count,seq in enumerate(sequences,1):
        a.write('>'+str(count)+'\n'+seq+'\n')


def append(ref_genome,upstream,seqlen,coors,filename):
    #This generates a dictionary that matches each site individually.
    
    contexts={}

    for chro in coors.keys():
        contexts[chro]={}

        #Upload genomic seqeunce for each chromosome.         
        z=open('store/STARdir/genomes/genomes_parsed_by_chr/'+ref_genome+'/'+chro+'.fa','r')
        for line in z:
            genome=line.strip()
        z.close()

        #Go through list of coordinates for matching sites and output motif sequence
        for strand in coors[chro].keys():
            contexts[chro][strand]={}
            for coor in coors[chro][strand]:
                if strand=='+':
                    sequence= str.upper(genome[coor-upstream:coor+seqlen]) 
                elif strand=='-':
                    sequence= complimentary.seq( str.upper(genome[coor-seqlen:coor+upstream]) )
                contexts[chro][strand][coor]=sequence
        del genome

    #Output context dictionary.
    return contexts

