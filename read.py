#read.py
#For reading text files delimited comma, tab or spaces.

def csv(line):
    return (line.strip()).split(',')

def tab(line):
    return (line.strip()).split('\t')

def space(line):
    from re import split as resplit
    regex='\s+'
    return resplit(regex,line.strip())

