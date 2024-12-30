# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 00:10:58 2022

@author: Sobhy
"""
'''
import pandas as pd
#========================== example 1 =========================================
infile = open(r'D:\Collage\Intrudction to Bio informatics\Sections\Section2\HAPPENN_dataset.fasta')
flag=0
tb=[]
for line in infile:
    if flag==0:
        s=line.split("|lcl|")
        print(s[3])
        flag=1
    else:
        print(line[:-1])
        flag=0
        if s[3]=='non-hemolytic' or s[3]=='non-hemolytic\n':
            tb.append([line[:-1],0])
        else:
            tb.append([line[:-1],1])
    print("--------------")

print(tb)
head=['Sequence','y']
df=pd.DataFrame(tb,columns=head)
df.to_csv(r"D:\Collage\Intrudction to Bio informatics\Sections\Section2\HAPPENN.csv")

'''
#========================== example 2 =========================================
import pandas as pd
infile = open(r"D:\Collage\Intrudction to Bio informatics\Sections\Section2\seq.fasta")
tb=[]
for line in infile:
    if line[0]==">":
        x=line[1:-1]
    else:
        seq=line[:-1]
        tb.append([x,seq])

print(tb)
head=['ID','Sequence']
df=pd.DataFrame(tb,columns=head)
df.to_csv(r"D:\Collage\Intrudction to Bio informatics\Sections\Section2\seq.csv")
'''
'''
#========================== example 2 =========================================
def GC_Content(seq):
    l=len(seq)
    num_G=seq.count("G")
    num_C=seq.count("C")
    total=num_C+num_G
    return total/l
