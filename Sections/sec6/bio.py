dec = {
    '$' : 0,
    'A' : 1,
    'C' : 2,
    'G' : 3,
    'T' : 4
}

T = 'ACGACTACGATAAC$'
table = []
i = 2**0
n = 0

# Initialize the table with the ranks of each character
row = [dec[T[j]] for j in range(len(T))]
table.append(row)

while True:
    l = []
    dec2 = {}
    
    # If i > 1, we need to process substrings of length i
    if i > 1:
        # Create list of unique substrings
        for j in range(len(T)):
            if not(table[n-1][j:j+i] in l):
                l.append(table[n-1][j:j+i])
        
        # Sort substrings
        l.sort()
        # print(f"Sorted substrings of length {i}: {l}")
        
        # Create new ranks for each substring
        for j in range(len(l)):
            dec2[tuple(l[j])] = j
    
    # Create the new row based on the previous ranks
    row = []
    for j in range(len(T)):
        if i == 1:
            row.append(dec[T[j]])
        else:
            row.append(dec2[tuple(table[n-1][j:j+i])])
    
    # Add the row to the table
    table.append(row)
    
    # Check if all ranks are unique in this row
    flag = 0
    for j in range(len(row)):
        c = row.count(j)
        if c > 1:
            flag = 1
            break
    
    # Print the current row for debugging
    print(f"Row after processing length {i}: {row}")
    print('#'*50)
    
    # If all ranks are unique, stop the loop
    if flag == 0:
        break
    
    # Prepare for the next iteration (doubling the length of the substrings)
    n += 1
    i = 2**n