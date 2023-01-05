#!/usr/bin/python
#a)Normalize endogenous m6A signal to spike-in m6A signal
#b)Normalize endogenous input signal at m6A site and a fixed nt upstream of it, to
# same thing in spike-in
#THEN normalize (a) to (b)
#Same as V4 but adds +1 to all endogenous input score so as to avoid 'NA's.
#NOTE THAT though this allows me to find more sites, it also means that in cases of
#methylase-ko, a m6A site that's lost might be simply caused by loss of the whole mRNA itself.
#Same as V5 but uses rmdup _PE read counts mapped to pg435 instead of _SE.
#Same as V6 but +1 to endogenous input (and now endogenous m6ACE) is now optional
#; also only output lines that have no 0-read endogenous inputs in ALL replicates (the other option)

#REQUIRES:
#1) bed_from_DESeq_v2 to find sites; use sites found in wildtype and/or mutant.
#2) PE read counts mapped to pg435 spike-in

#NOTE: code is to be ran via 2 parts.


import read

dest='store/STARdir/'
paired='_PE'

part=2#1 is to generate individual output per lib. 2 is to make consolidated output for final_annotation.
z=0#input individual list index here. Only relevant for part 1.
avoidZero=1.0#Do 1.0 for TSS-m6Am, 0.0 for internal m6A/m6Am. Only relevant for part 2.

upstream=51 #Pair-end 75nt minus 8nt for UMI, minus minimum 20nt for mapping, plus 4nt since looking at read count up to 4nt upstream of m6A.
spike_site=21#Position of m6A in spike-in


#INPUT sample and LIB NAMES HERE!
#[Input,m6ACE] corresponding to sample type.

##samples=['Empty-E2','Mettl3-134B4']
######libs=[['LSG363','LSG368'],['LSG364','LSG369'],['LSG379','LSG374'],['LSG344','LSG349'],['LSG365','LSG370'],['LSG366','LSG371'],['LSG380','LSG375'],['LSG346','LSG351'],]#Quadplicate miseq lib
##libs=[['LSG363+LSG407','LSG368+LSG408'],['LSG364+LSG409','LSG369+LSG410'],['LSG379+LSG411','LSG374+LSG412'],['LSG365+LSG413','LSG370+LSG414'],['LSG366+LSG415','LSG371+LSG416'],['LSG380+LSG417','LSG375+LSG418'],]#Nextseq triplicate
####
##samples=['Empty-E2','Empty-E2']
##libs=[['LSG363+LSG407','LSG368+LSG408'],['LSG364+LSG409','LSG369+LSG410'],['LSG379+LSG411','LSG374+LSG412'],['LSG363+LSG407','LSG368+LSG408'],['LSG364+LSG409','LSG369+LSG410'],['LSG379+LSG411','LSG374+LSG412'],]#Nextseq triplicate

##samples=['Empty-E2','FTO-168C4']
##libs=[['LSG363','LSG368'],['LSG364','LSG369'],['LSG379','LSG374'],['LSG398','LSG404'],['LSG399','LSG405'],['LSG400','LSG406'],]#Triplicate miseq lib
##libs=[['LSG363+LSG407','LSG368+LSG408'],['LSG364+LSG409','LSG369+LSG410'],['LSG379+LSG411','LSG374+LSG412'],['LSG398+LSG419','LSG404+LSG420'],['LSG399+LSG421','LSG405+LSG422'],['LSG400+LSG423','LSG406+LSG424'],]#Nextseq triplicate

##samples=['Empty-E2','Alkbh5-7G4']
##libs=[['LSG363+LSG407','LSG368+LSG408'],['LSG364+LSG409','LSG369+LSG410'],['LSG379+LSG411','LSG374+LSG412'],\
##      ['LSG425+LSG468','LSG428+LSG471'],['LSG426+LSG469','LSG429+LSG472'],['LSG427+LSG470','LSG430+LSG473'],]#Triplicate nextseq lib

##samples=['Empty-E2','Mettl4-63A3']
##libs=[['LSG363+LSG407','LSG368+LSG408'],['LSG364+LSG409','LSG369+LSG410'],['LSG379+LSG411','LSG374+LSG412'],\
##      ['LSG449+LSG477','LSG446+LSG474'],['LSG450+LSG478','LSG447+LSG475'],['LSG451+LSG479','LSG448+LSG476'],]#Triplicate nextseq lib

##samples=['Empty-E2','Pcif1-A5']
##libs=[['LSG363+LSG407','LSG368+LSG408'],['LSG364+LSG409','LSG369+LSG410'],['LSG379+LSG411','LSG374+LSG412'],\
##      ['LSG455+LSG483','LSG454+LSG482'],['LSG486','LSG487'],['LSG488','LSG489'],]#Triplicate nextseq lib
##
##samples=['Hek-cyto','Hek-nucl']
##libs=[['LSG432+LSG491','LSG431+LSG490'],['LSG434+LSG493','LSG433+LSG492'],['LSG440+LSG499','LSG439+LSG498'],\
##      ['LSG436+LSG495','LSG435+LSG494'],['LSG438+LSG497','LSG437+LSG496'],['LSG442+LSG501','LSG441+LSG500'],]#Triplicate nextseq lib
##
##samples=['Scr-B9','Mettl16-A2']
##libs=[['LSG459+LSG485','LSG458+LSG484'],['LSG529+LSG547','LSG523+LSG541'],['LSG530+LSG548','LSG524+LSG542'],\
##      ['LSG453+LSG481','LSG452+LSG480'],['LSG525+LSG543','LSG519+LSG537'],['LSG505','LSG504'],]#Triplicate nextseq lib

##samples=['Scr-B9','Pcif1-A5']
##libs=[['LSG459+LSG485','LSG458+LSG484'],['LSG529+LSG547','LSG523+LSG541'],['LSG530+LSG548','LSG524+LSG542'],\
##      ['LSG455+LSG483','LSG454+LSG482'],['LSG527+LSG545','LSG521+LSG539'],['LSG528+LSG546','LSG522+LSG540'],]#Triplicate nextseq lib
####

##samples=['Empty-E2','Pcif1-183B2']
##libs=[['LSG363+LSG407','LSG368+LSG408'],['LSG364+LSG409','LSG369+LSG410'],['LSG379+LSG411','LSG374+LSG412'],\
##      ['LSG534','LSG531'],['LSG535','LSG532'],['LSG536','LSG533'],]#Triplicate nextseq lib

##samples=['Empty-E2','Mettl16-A2']
##libs=[['LSG363+LSG407','LSG368+LSG408'],['LSG364+LSG409','LSG369+LSG410'],['LSG379+LSG411','LSG374+LSG412'],\
##      ['LSG453+LSG481','LSG452+LSG480'],['LSG525+LSG543','LSG519+LSG537'],['LSG505','LSG504'],]#Triplicate nextseq lib

##samples=['siRNA-screen',]
##libs=[['LSG459','LSG458'],\
##      ['LSG453','LSG452'],\
##      ['LSG455','LSG454'],\
##      ['LSG457','LSG456'],\
##      ['LSG461','LSG460'],\
##      ['LSG463','LSG462'],\
##      ['LSG465','LSG464'],\
##      ['LSG467','LSG466']]

###Uses Empty-E2_condition1_m6A_peaks.bed (nextseq)
##samples=['crispr-screen',]
##libs=[['LSG363','LSG368'],\
##      ['LSG364','LSG369'],\
##      ['LSG379','LSG374'],\
##      ['LSG367','LSG372'],\
##      ['LSG378','LSG373'],\
##      ['LSG344','LSG349'],\
##      ['LSG348','LSG353'],\
##      ['LSG340','LSG341'],\
##      ['LSG342','LSG343']]

##samples=['Mettl4-63A3-OE-WT-Mettl4','Mettl4-63A3-OE-CD-Mettl4']
##libs=[['LSG555+LSG567','LSG549+LSG561'],\
##      ['LSG556+LSG568','LSG550+LSG562'],\
##      ['LSG557+LSG569','LSG551+LSG563'],\
##      ['LSG558+LSG570','LSG552+LSG564'],\
##      ['LSG559+LSG571','LSG553+LSG565'],\
##      ['LSG560+LSG572','LSG554+LSG566']]

###For 'Quantitative-Mix'
##samples=['Empty-E2-2','Mettl3-134B4-2']
##libs=[['LSG590+LSG600','LSG585+LSG595+LSG637'],\
##      ['LSG610','LSG605'],\
##      ['LSG620','LSG615'],\
##
##      ['LSG591+LSG601','LSG586+LSG596'],\
##      ['LSG611','LSG606'],\
##      ['LSG621','LSG616'],\
##
##      ['LSG592+LSG602','LSG587+LSG597'],\
##      ['LSG612','LSG607'],\
##      ['LSG622','LSG617'],\
##
##      ['LSG593+LSG603','LSG588+LSG598'],\
##      ['LSG613','LSG608'],\
##      ['LSG623','LSG618'],\
##      
##      ['LSG594+LSG604','LSG589+LSG599'],\
##      ['LSG614','LSG609'],\
##      ['LSG624','LSG619']]

###For just 2nd replicate of empty versus mettl3-ko (only need to do part 2)
##samples=['Empty-E2-2','Mettl3-134B4-2']
##libs=[['LSG590+LSG600','LSG585+LSG595+LSG637'],\
##      ['LSG610','LSG605'],\
##      ['LSG620','LSG615'],\
##
##      ['LSG594+LSG604','LSG589+LSG599'],\
##      ['LSG614','LSG609'],\
##      ['LSG624','LSG619']]

###Compare Empty-E2 triplicate 1 with Empty-E2 triplicate 2
##samples=['Empty-E2','Empty-E2-2']
##libs=[['LSG363+LSG407','LSG368+LSG408'],\
##      ['LSG364+LSG409','LSG369+LSG410'],\
##      ['LSG379+LSG411','LSG374+LSG412'],\
##      ['LSG590+LSG600','LSG585+LSG595+LSG637'],\
##      ['LSG610','LSG605'],\
##      ['LSG620','LSG615']]

###Compare Empty-E2 with FL-Fto overexpression in Empty-E2 cells.
##samples=['Empty-E2','Empty-E2-OE-FL-Fto']
##libs=[['LSG363+LSG407','LSG368+LSG408'],\
##      ['LSG364+LSG409','LSG369+LSG410'],\
##      ['LSG379+LSG411','LSG374+LSG412'],\
##      ['LSG631','LSG625'],\
##      ['LSG632','LSG626'],\
##      ['LSG633','LSG627']]
##
###Compare Empty-E2 with Mettl5-71C2.
##samples=['Empty-E2','Mettl5-71C2']
##libs=[['LSG363+LSG407','LSG368+LSG408'],\
##      ['LSG364+LSG409','LSG369+LSG410'],\
##      ['LSG379+LSG411','LSG374+LSG412'],\
##      ['LSG634','LSG628'],\
##      ['LSG635','LSG629'],\
##      ['LSG636','LSG630']]

###Compare Empty-E2 with Empty-E2-abcam
##samples=['Empty-E2','Empty-E2-abcam']
##libs=[['LSG363+LSG407','LSG368+LSG408'],\
##      ['LSG364+LSG409','LSG369+LSG410'],\
##      ['LSG379+LSG411','LSG374+LSG412'],\
##      ['LSG648','LSG642'],\
##      ['LSG649','LSG643'],\
##      ['LSG650','LSG644']]
##
###Compare Empty-E2 with Empty-E2-sysy202111
##samples=['Empty-E2','Empty-E2-sysy202111']
##libs=[['LSG363+LSG407','LSG368+LSG408'],\
##      ['LSG364+LSG409','LSG369+LSG410'],\
##      ['LSG379+LSG411','LSG374+LSG412'],\
##      ['LSG651','LSG645'],\
##      ['LSG652','LSG646'],\
##      ['LSG653','LSG647']]
####
###Compare Empty-E2 with Hemk2-165D2
##samples=['Empty-E2','Hemk2-165D2']
##libs=[['LSG363+LSG407','LSG368+LSG408'],\
##      ['LSG364+LSG409','LSG369+LSG410'],\
##      ['LSG379+LSG411','LSG374+LSG412'],\
##      ['LSG657','LSG654'],\
##      ['LSG658','LSG655'],\
##      ['LSG659','LSG656']]

###HCT116
##samples=['HCT116','HCT116']
##libs=[['LSG729','LSG726'],\
##      ['LSG730','LSG727'],\
##      ['LSG731','LSG728'],
##      ['LSG729','LSG726'],\
##      ['LSG730','LSG727'],\
##      ['LSG731','LSG728'],]

###For 2nd replicate of empty versus pcif1-ko
##samples=['Empty-E2-2','Pcif1-183B2']
##libs=[['LSG590+LSG600','LSG585+LSG595+LSG637'],['LSG610','LSG605'],['LSG620','LSG615'],\
##      ['LSG534','LSG531'],['LSG535','LSG532'],['LSG536','LSG533'],]

##samples=['9DIV','19DIV']
##libs=[['LSG777','LSG773'],['LSG778','LSG774'],\
##      ['LSG779','LSG775'],['LSG780','LSG776'],]

##samples=['Empty-E2+IAV','Pcif1-183B2+IAV']
##libs=[['LSG823','LSG829'],['LSG824','LSG830'],['LSG825','LSG831'],\
##      ['LSG826','LSG832'],['LSG827','LSG833'],['LSG828','LSG834']]

samples=['W1118-ovary','Mettl3-KO-ovary']
libs=[['LSG867','LSG861'],['LSG868','LSG862'],['LSG869','LSG863'],\
      ['LSG870','LSG864'],['LSG871','LSG865'],['LSG872','LSG866']]

##samples=['WT-S2','Mettl3-KO-S2']
##libs=[['LSG879','LSG873'],['LSG880','LSG874'],['LSG881','LSG875'],\
##      ['LSG882','LSG876'],['LSG883','LSG877'],['LSG884','LSG878']]



##Running of script from here onwards.
filename=''
for sample in samples:
    filename += sample+'_'


if part==1:
    libs=[libs[z],]
    for lib in libs:
        print (lib)
        #First, get list of significant m6A/m6Am sites, in BOTH sample types.
        sites={}
        sites['+']={}
        sites['-']={}
        for sample in samples:
            y=open('store/DESeq2/'+sample+'_good_m6A_peaks.bed','r')
            for i in range(3):
                next(y)
            for line in y:
                ar=read.tab(line)
                chro=ar[0]
                coor=int(ar[2])
                strand=ar[5]
                #2 columns are for endogenous input and m6A signal respectively.
                try:
                    sites[strand][chro][coor]=[0,0]
                except KeyError:
                    sites[strand][chro]={}
                    sites[strand][chro][coor]=[0,0]
            y.close()

        #Get spike-in input signal
        spike_in_input=0
        a=open('store/BT2/'+lib[0]+'_clipped_identified_PE_5prime.bg','r')
        for i in range(3):
            next(a)
        for line in a:
            ar=read.tab(line)
            if 1<=int(ar[2])<=spike_site:
                if int(ar[3])>0:
                    spike_in_input+=int(ar[3])
        a.close()

        #Get spike-in m6A signal
        spike_in_m6A=0
        a=open('store/BT2/'+lib[1]+'_clipped_identified_PE_5prime.bg','r')
        for i in range(3):
            next(a)
        for line in a:
            ar=read.tab(line)
            if (spike_site-4) <= int(ar[2]) <=spike_site:
                if int(ar[3])>0:
                    spike_in_m6A+=int(ar[3])
        a.close()

        #Get endogenous input signal
        a=open(dest+lib[0]+paired+'/'+lib[0]+paired+'_5prime.bg','r')
        for i in range(3):
            next(a)
        for line in a:
            ar=read.tab(line)
            chro=ar[0]
            coor=int(ar[2])
            score=int(ar[3])
            if score>0 and ( chro in sites['+'] ):
                for key in sites['+'][chro]:
                    if (key-upstream)<=coor<=key:
                        sites['+'][chro][key][0]+=score
            elif score<0 and ( chro in sites['-'] ):
                score=-score
                for key in sites['-'][chro]:
                    if key<=coor<=(key+upstream):
                        sites['-'][chro][key][0]+=score 
        a.close()

        
        #Get endogenous m6A signal
        a=open(dest+lib[1]+paired+'/'+lib[1]+paired+'_5prime.bg','r')
        for i in range(3):
            next(a)
        for line in a:
            ar=read.tab(line)
            chro=ar[0]
            coor=int(ar[2])
            score=int(ar[3])
            if score>0 and ( chro in sites['+'] ):
                for key in sites['+'][chro]:
                    if (key-4)<=coor<=key:
                        sites['+'][chro][key][1]+=score
            elif score<0 and ( chro in sites['-'] ):
                score=-score
                for key in sites['-'][chro]:
                    if key<=coor<=(key+4):
                        sites['-'][chro][key][1]+=score 
        a.close()

        #Write out output.
        b=open('store/quanti/'+filename+'_'+lib[0]+'_'+lib[1]+'_'+'output.txt','w')
        for strand in sites:
            for chro in sites[strand]:
                for coor in sites[strand][chro]:
                    b.write(chro+'\t'+str(coor-1)+'\t'+str(coor)+'\t'+strand+'\t'+\
                            str(spike_in_input)+'\t'+str(spike_in_m6A)+'\t'+\
                            str(sites[strand][chro][coor][0])+'\t'+\
                            str(sites[strand][chro][coor][1])+'\n')
        print (filename, lib,'done')


if part==2:
    ###To output quantifications into columns for viewing in excel.
    fulls={}
    count=0
    for lib in libs:
        a=open('store/quanti/'+filename+'_'+lib[0]+'_'+lib[1]+'_output.txt','r')
        for line in a:
            ar=read.tab(line)
            if float(ar[6]) != 0.0 and avoidZero==0.0:
                try:
                    fulls[ar[0],ar[1],ar[2],ar[3]][count]= ( float(ar[7]) / float(ar[6])  ) \
                                                           / ( float(ar[5]) / float(ar[4]) )
                except KeyError:
                    fulls[ar[0],ar[1],ar[2],ar[3]]=[0]*len(libs)           
                    fulls[ar[0],ar[1],ar[2],ar[3]][count]= ( float(ar[7]) / float(ar[6])  ) \
                                                           / ( float(ar[5]) / float(ar[4]) )
            elif avoidZero==1.0:
                try:
                    fulls[ar[0],ar[1],ar[2],ar[3]][count]= ( (float(ar[7]) + avoidZero) / (float(ar[6]) + avoidZero) ) \
                                                           / ( float(ar[5]) / float(ar[4]) )
                except KeyError:
                    fulls[ar[0],ar[1],ar[2],ar[3]]=[0]*len(libs)           
                    fulls[ar[0],ar[1],ar[2],ar[3]][count]= ( (float(ar[7]) + avoidZero) / (float(ar[6]) + avoidZero) ) \
                                                           / ( float(ar[5]) / float(ar[4]) )
        count+=1


    b=open('store/quanti/'+filename+'consolidate_output.txt','w')
    for key in fulls:
        to_write=1
        for score in fulls[key]:
            if str(score) == 0:
                to_write=0
                break
        if to_write==1:
            b.write(key[0]+'\t'+key[1]+'\t'+key[2]+'\t'+key[3]+'\t')
            for score in fulls[key]:
                b.write(str(score)+'\t')
            b.write('\n')

    #Write a log file to log which avoidZero option was used.
    c=open('store/quanti/'+filename+'parameters.txt','w')
    c.write('avoidZero\t'+str(avoidZero)+'\n')
    

            

