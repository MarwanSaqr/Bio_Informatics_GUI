# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 19:37:21 2022

@author: Mohamed Elhakim
"""
#import numpy as np
import bisect

def IndexSorted(seq,k):
    index = []
    for i in range(len(seq)-k+1):
        index.append((seq[i:i+k], i))
    index.sort() 
    return index

def query(t,p,index):
    keys = [r[0] for r in index]
    st = bisect.bisect_left(keys,p[:len(keys[0])])
    en = bisect.bisect(keys,p[:len(keys[0])])
    hits = index[st:en] 
    l=[h[1] for h in hits ]
    offsets=[]
    for i in l:
        if t[i:i+len(p)]==p:
            offsets.append(i)
    return offsets
        
file=open(r"D:\Collage\Intrudction to Bio informatics\Sections\sec5\dna1.fasta")
l=[i for i in file]
t=l[1][0:-1]
#index=IndexSorted(t,3)
t="GCTACGATCTAGAATCTA"
index=IndexSorted(t,2)
#p="GCGTCGCTGTGGAG"
p="TCTA"
print(query(t,p,index))
#print(t[7:11])


