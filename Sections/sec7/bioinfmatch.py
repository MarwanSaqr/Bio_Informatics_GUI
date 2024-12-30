# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 19:37:21 2022

@author: Mohamed Elhakim
"""
#import numpy as np
from itertools import permutations

def overlap(a, b, min_length=3):
    start=0
    while True:
        start=a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1
            
               
def native_overlap(reads, k):
    olap={}
    for a,b in permutations(reads, 2):
        olen=overlap(a, b, k)
        if olen > 0:
            olap[(a, b)]=olen
    return olap

print(native_overlap(["ACGGTA", "GGTACC","GCATACG"], 3))
