#complimentary.py
#Just to make complimentary sequences.

def seq(line):
    comp={}
    comp['A']='T'
    comp['G']='C'
    comp['C']='G'
    comp['T']='A'
    comp['U']='A'
    comp['N']='N'
    
    antisense=''
    sense=str.upper(line.strip())
    for i in range(len(sense)):
        antisense+=comp[sense[i]]
    antisense=antisense[::-1]
    return antisense
