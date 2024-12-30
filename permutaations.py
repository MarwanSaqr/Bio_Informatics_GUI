from itertools import permutations
reads=["ACGGTA", "GGTACC","GCATACG"]
for a,b in permutations(reads, 2):
    print(f"a is equal {a} and b is equal {b}")