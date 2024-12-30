def GC_Content(seq):
    l = len(seq)
    num_G = seq.count("G")
    num_C = seq.count("C")
    total = num_C + num_G
    return total / l if l > 0 else 0  # Avoid division by zero

def Complement(seq):
    dic = {"G": "C", "C": "G", "A": "T", "T": "A"}
    return ''.join(dic[base] for base in seq)

def Reverse(seq):
    return seq[::-1]

def Reverse_Complement(seq):
    return Complement(Reverse(seq))

def Translation_Table(seq):
    dic = {
        "TTT": "F", "CTT": "L", "ATT": "I", "GTT": "V",
        "TTC": "F", "CTC": "L", "ATC": "I", "GTC": "V",
        "TTA": "L", "CTA": "L", "ATA": "I", "GTA": "V",
        "TTG": "L", "CTG": "L", "ATG": "M", "GTG": "V",
        "TCT": "S", "CCT": "P", "ACT": "T", "GCT": "A",
        "TCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
        "TCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
        "TCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
        "TAT": "Y", "CAT": "H", "AAT": "N", "GAT": "D",
        "TAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
        "TAA": "*", "CAA": "Q", "AAA": "K", "GAA": "E",
        "TAG": "*", "CAG": "Q", "AAG": "K", "GAG": "E",
        "TGT": "C", "CGT": "R", "AGT": "S", "GGT": "G",
        "TGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
        "TGA": "*", "CGA": "R", "AGA": "R", "GGA": "G",
        "TGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"
    }
    s = ""
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if codon in dic:
            s += dic[codon]
        else:
            s += "X"  # Handle unknown codons
    return s

# Main script
try:
    with open(r"D:\Collage\Intrudction to Bio informatics\Sections\sec3\dna2.fna") as file:
        lines = file.readlines()
        s = ''.join(line.strip() for line in lines if not line.startswith('>'))  # Ignore header lines

    print("GC Content:", GC_Content(s))
    print("Reverse:", Reverse(s))
    print("Reverse Complement:", Reverse_Complement(s))
    print("Translation:", Translation_Table(s))

except FileNotFoundError:
    print("Error: File not found. Please check the file path.")
