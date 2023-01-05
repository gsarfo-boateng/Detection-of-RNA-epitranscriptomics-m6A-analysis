#!/usr/bin/python
#Format data from 5prime.bg for input in DESeq2.
#Same as v1 but creates 2 outputs, 1 for 'Aonly' and 1 for 'all'.

import read
from collections import defaultdict
import complimentary

#Parameters
min_rpm=0.05#This is cutoff where (read_count*1e6/total mapped reads in all libraries) must >= (min_rpm / number of libraries)
dest='store/STARdir/'
paired='_PE'
ref_genome='dm6'
upstream=0
seqlen=1

#Arrange as input, m6ace, input, m6ace, input, m6ace....
##
##sample='Empty-E2'
##libs=['LSG363+LSG407',\
##      'LSG368+LSG408',\
##      'LSG364+LSG409',\
##      'LSG369+LSG410',\
##      'LSG379+LSG411',\
##      'LSG374+LSG412',]

##sample='Mettl3-134B4'
##libs=['LSG365+LSG413',\
##      'LSG370+LSG414',\
##      'LSG366+LSG415',\
##      'LSG371+LSG416',\
##      'LSG380+LSG417',\
##      'LSG375+LSG418',]

##sample='FTO-168C4'
##libs=['LSG398+LSG419',\
##      'LSG404+LSG420',\
##      'LSG399+LSG421',\
##      'LSG405+LSG422',\
##      'LSG400+LSG423',\
##      'LSG406+LSG424',]

##sample='Mettl4-63A3'
##libs=['LSG449+LSG477',\
##      'LSG446+LSG474',\
##      'LSG450+LSG478',\
##      'LSG447+LSG475',\
##      'LSG451+LSG479',\
##      'LSG448+LSG476']

##sample='Alkbh5-7G4'
##libs=['LSG425+LSG468',\
##      'LSG428+LSG471',\
##      'LSG426+LSG469',\
##      'LSG429+LSG472',\
##      'LSG427+LSG470',\
##      'LSG430+LSG473',]
####
##sample='Pcif1-A5'
##libs=['LSG455+LSG483',\
##      'LSG454+LSG482',\
##      'LSG527+LSG545',\
##      'LSG521+LSG539',\
##      'LSG528+LSG546',\
##      'LSG522+LSG540',]

##sample='Scr-B9'
##libs=['LSG459+LSG485',\
##      'LSG458+LSG484',\
##      'LSG529+LSG547',\
##      'LSG523+LSG541',\
##      'LSG530+LSG548',\
##      'LSG524+LSG542',]

##sample='Mettl16-A2'
##libs=['LSG453+LSG481',\
##      'LSG452+LSG480',\
##      'LSG525+LSG543',\
##      'LSG519+LSG537',\
##      'LSG505',\
##      'LSG504']

##sample='Hek-cyto'
##libs=['LSG432+LSG491',\
##      'LSG431+LSG490',\
##      'LSG434+LSG493',\
##      'LSG433+LSG492',\
##      'LSG440+LSG499',\
##      'LSG439+LSG498',]

##sample='Hek-nucl'
##libs=['LSG436+LSG495',\
##      'LSG435+LSG494',\
##      'LSG438+LSG497',\
##      'LSG437+LSG496',\
##      'LSG442+LSG501',\
##      'LSG441+LSG500',]

##sample='siRNA-screen'
##libs=['LSG459','LSG458',\
##      'LSG453','LSG452',\
##      'LSG455','LSG454',\
##      'LSG457','LSG456',\
##      'LSG461','LSG460',\
##      'LSG463','LSG462',\
##      'LSG465','LSG464',\
##      'LSG467','LSG466',]

##sample='Pcif1-183B2'
##libs=['LSG534','LSG531',\
##      'LSG535','LSG532',\
##      'LSG536','LSG533',]

##sample='Mettl4-63A3-OE-WT-Mettl4'
##libs=['LSG555+LSG567',\
##      'LSG549+LSG561',\
##      'LSG556+LSG568',\
##      'LSG550+LSG562',\
##      'LSG557+LSG569',\
##      'LSG551+LSG563']

##sample='Mettl4-63A3-OE-CD-Mettl4'
##libs=['LSG558+LSG570',\
##      'LSG552+LSG564',\
##      'LSG559+LSG571',\
##      'LSG553+LSG565',\
##      'LSG560+LSG572',\
##      'LSG554+LSG566']

##sample='Empty-E2-2'
##libs=['LSG590+LSG600',\
##      'LSG585+LSG595+LSG637',\
##      'LSG610',\
##      'LSG605',\
##      'LSG620',\
##      'LSG615']

##sample='Mettl3-134B4-2'
##libs=['LSG594+LSG604',\
##      'LSG589+LSG599',\
##      'LSG614',\
##      'LSG609',\
##      'LSG624',\
##      'LSG619']

##sample='Empty-E2-OE-FL-Fto'
##libs=['LSG631',\
##      'LSG625',\
##      'LSG632',\
##      'LSG626',\
##      'LSG633',\
##      'LSG627']

##sample='Mettl5-71C2'
##libs=['LSG634',\
##      'LSG628',\
##      'LSG635',\
##      'LSG629',\
##      'LSG636',\
##      'LSG630']

##sample='Empty-E2-abcam'
##libs=['LSG648',\
##      'LSG642',\
##      'LSG649',\
##      'LSG643',\
##      'LSG650',\
##      'LSG644']

##sample='Empty-E2-sysy202111'
##libs=['LSG651',\
##      'LSG645',\
##      'LSG652',\
##      'LSG646',\
##      'LSG653',\
##      'LSG647']
####
##sample='Hemk2-165D2'
##libs=['LSG657',\
##      'LSG654',\
##      'LSG658',\
##      'LSG655',\
##      'LSG659',\
##      'LSG656']

##sample='HCT116'
##libs=['LSG729','LSG726',\
##      'LSG730','LSG727',\
##      'LSG731','LSG728',]

##sample='9DIV'
##libs=['LSG777','LSG773',\
##      'LSG778','LSG774']

##sample='19DIV'
##libs=['LSG779','LSG775',\
##      'LSG780','LSG776']

##sample='Empty-E2+IAV'
##libs=['LSG823','LSG829',\
##      'LSG824','LSG830',\
##      'LSG825','LSG831']
##
##sample='Pcif1-183B2+IAV'
##libs=['LSG826','LSG832',\
##      'LSG827','LSG833',\
##      'LSG828','LSG834']

sample='W1118-ovary'
libs=['LSG867','LSG861',\
      'LSG868','LSG862',\
      'LSG869','LSG863']

sample='Mettl3-KO-ovary'
libs=['LSG870','LSG864',\
      'LSG871','LSG865',\
      'LSG872','LSG866']
####
sample='WT-S2'
libs=['LSG879','LSG873',\
      'LSG880','LSG874',\
      'LSG881','LSG875']
####
sample='Mettl3-KO-S2'
libs=['LSG882','LSG876',\
      'LSG883','LSG877',\
      'LSG884','LSG878']

matrix={}
sites={}
#Actual run

#Collect data first.
total=0
for lib in libs:
    matrix[lib]=defaultdict(int)
    a=open(dest+lib+paired+'/'+lib+paired+'_5prime.bg','r')
    for i in range(3):
        next(a)
    for line in a:
        ar=read.tab(line)
        chro=str(ar[0])
        start=str(ar[1])
        end=str(ar[2])
        score=int(ar[3])
        total += abs(score)
        site=chro+':'+start+'-'+end
        if score>0:
            site=site+'+'
        else:
            site=site+'-'
            score=-score
        sites[site]=True
        matrix[lib][site]=score
    a.close()
    print (lib,'Original',len(sites.keys()))
min_read = min_rpm * float(total) / 1000000.0

#Output for DESeq2 '_all'
temps=sites
sites={}
for temp in temps.keys():
    sites[temp]=True
print ('all',len(sites.keys()))

finalcount=0
b=open('store/DESeq2/'+sample+'_for_DESeq2_all.txt','w')
for lib in libs:
    b.write('\t'+lib)
b.write('\n')
for site in sorted(sites.keys()):
    check=0
    for lib in libs:
        #This part is to remove anything with baseMean<min_rpm since those will end up with padj='NA' anyway.
        check+=matrix[lib][site]
    if check>=min_read:    
        b.write(site)
        for lib in libs:
            b.write('\t'+str(matrix[lib][site]))
        b.write('\n')
        finalcount+=1
print (sample,'Final',finalcount)


#Output for DESeq2 '_Aonly'
sites={}
chros=[]
y=open('store/STARdir/genomes/'+ref_genome+'/'+ref_genome+'_chromosomes.txt','r')
for line in y:
    chros+=[line.strip(),]
y.close()

#Use parse_chromosomes.py to parse fasta by chromosomes.

for chro in chros:
    many=0
    z=open('store/STARdir/genomes/genomes_parsed_by_chr/'+ref_genome+'/'+chro+'.fa','r')
    for line in z:
        genome=line.strip()
    z.close()
    passed=False
    for temp in sorted(temps.keys()):
        if chro == temp.split(':')[0]:
            passed=True
            many+=1
            strand=temp[-1]
            if strand == '+':
                start=int(temp[0:-1].split(':')[1].split('-')[0])
                sequence=str.upper(genome[start-upstream:start+seqlen])           
            elif strand == '-':
                end=int(temp[0:-1].split(':')[1].split('-')[1])
                sequence=str.upper(genome[end-seqlen:end+upstream])
                sequence=complimentary.seq(sequence)
            if sequence=='A':
                sites[temp]=True
            del temps[temp]
        elif ( chro != temp.split(':')[0] ) and ( passed==True ) :
            break

    print (chro+'\t'+str(many))
del temps
del genome
print ('Aonly',len(sites.keys()))

finalcount=0
b=open('store/DESeq2/'+sample+'_for_DESeq2_Aonly.txt','w')
for lib in libs:
    b.write('\t'+lib)
b.write('\n')
for site in sorted(sites.keys()):
    check=0
    for lib in libs:
        #This part is to remove anything with baseMean<min_rpm since those will end up with padj='NA' anyway.
        #DON'T FORGET THAT THIS CUTOFF ALSO INCLUDES INPUT READS (that's why have to set lower)
        check+=matrix[lib][site]
    if check>=min_read:    
        b.write(site)
        for lib in libs:
            b.write('\t'+str(matrix[lib][site]))
        b.write('\n')
        finalcount+=1
print (sample,'Final',finalcount)


