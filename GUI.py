from tkinter import *
from tkinter import ttk, messagebox
from tkinter.filedialog import askopenfilename
import pandas as pd
import os
import numpy as np
from PIL import Image, ImageTk
import cv2
import bisect
from itertools import permutations


# ################### Creating the GUI ###################
root = Tk()
root.title("Bio Methods")
Font = ('Berlin Sans FB Demi', 20, 'bold')
Bg="lightblue"
Fg="black"
root.configure(bg="#fff")
Padx = 2
Pady = 5
ButtonWidth = 17
ButtonHeight = 1

# Load MP4 video
video_path = r"DNA_BACK.mp4" 
cap = cv2.VideoCapture(video_path)

# Function to play video
def play_video():
    ret, frame = cap.read()
    if ret:
        # Convert frame to RGB
        frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
        # Convert frame to ImageTk format
        frame_image = ImageTk.PhotoImage(image=Image.fromarray(frame))
        video_label.config(image=frame_image)
        video_label.image = frame_image
        root.after(10, play_video)  # Adjust to your video FPS
    else:
        cap.set(cv2.CAP_PROP_POS_FRAMES, 0)  # Restart video
        play_video()


video_label = Label(root)
video_label.pack(fill="both", expand=True)


def change_ico(x):
     x.iconbitmap("DNA.ico")
change_ico(root)

def center_screen(window):
    screen_width = window.winfo_screenwidth()
    screen_height = window.winfo_screenheight()
    window_width = screen_width - 470
    window_height = screen_height -250
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2
    window.geometry(f"{window_width}x{window_height}+{x}+{y}")



# Global variable to store the directory of the selected file
path_of_browserfile = ""

# Function to browse files
def browse_file():
    global path_of_browserfile
    file_name = askopenfilename(filetypes=[("Fasta Files", "*.fasta"),("Fna Files", "*.fna"), ("All Files", "*.*")])
    if file_name:
        file_path_entry.delete(0, END)
        file_path_entry.insert(0, file_name)
 
# Function to process and display the CSV data
def execute_and_show(kind_of_sequence):
    file_data = file_path_entry.get()
    extract_data = extract_data_entry.get()

    if not file_data or not extract_data:
        status_label.config(text="Please provide all inputs!", fg="red")
        return

    try:
        # Validate output file extension
        if not extract_data.endswith(".csv"):
            extract_data += ".csv"

        # Output file path
        output_path = os.path.join(path_of_browserfile, extract_data)

        # Process based on the selected kind of sequence
        if kind_of_sequence == 1:
            CSV1(file_data, output_path)

        elif kind_of_sequence ==2:
            CSV2(file_data, output_path)
      
        status_label.config(text="File processed successfully!", fg="green")

        # Load and display the CSV content
        df = pd.read_csv(output_path)
        for widget in display_frame.winfo_children():
            widget.destroy()

        # Create Treeview to display data
        columns = list(df.columns)
        tree = ttk.Treeview(display_frame, columns=columns, show='headings', height=10)

        for col in columns:
            tree.heading(col, text=col)
            tree.column(col, width=150, anchor='center')

        for _, row in df.iterrows():
            tree.insert("", "end", values=list(row))

        tree.pack(fill="both", expand=True)

    except Exception as e:
        status_label.config(text=f"Error: {str(e)}", fg="red")

# Function for CSV1 processing
def CSV1(file_data, output_path):
    try:
        with open(file_data, 'r') as infile:
            flag = 0
            tb = []
            for line in infile:
                if flag == 0:
                    s = line.split("|lcl|")
                    flag = 1
                else:
                    flag = 0
                    # Ensure the list tb gets correctly formatted rows
                    if len(s) > 3:  # Ensure index 3 exists in s
                        tb.append([line.strip(), 0 if s[3].strip() == 'non-hemolytic' else 1])
                    else:
                        raise ValueError("File format is invalid or unexpected.")
        
        # Ensure the DataFrame is correctly created
        head = ['Sequence', 'y']
        df = pd.DataFrame(tb, columns=head)
        
        # Save to CSV
        df.to_csv(output_path, index=False)
    except Exception as e:
            print(f"Error in sequence_DNA: {e}")
            raise
# Function for CSV2 processing
def CSV2(file_data, output_path):
    try:
        with open(file_data, 'r') as infile:
            tb = []
            for line in infile:
                if line[0]==">":
                    x=line[1:-1]
                else:
                    seq=line[:-1]
                    tb.append([x,seq])

        df = pd.DataFrame(tb, columns=["ID", "Sequence"])
        df.to_csv(output_path, index=False)
    except Exception as e:
        raise ValueError(f"Error in CSV2: {str(e)}")


def display_overlab():
    for widget in root.winfo_children():
        widget.destroy()

    #root.title("Bioinformatics Overlap GUI")
    # Create a frame for the buttons
    button_frame = Frame(root, bg="#fff", padx=10, pady=10)
    button_frame.pack(fill="x", pady=5)

    # Buttons inside the frame
    Label(button_frame, text="Reads (comma-separated):").grid(row=0, column=0, sticky="w")
    reads_entry = Entry(button_frame, width=50)
    reads_entry.grid(row=0, column=1, padx=5)

    Label(button_frame, text="Minimum Overlap Length:").grid(row=1, column=0, sticky="w")
    k_entry = Entry(button_frame, width=10)
    k_entry.grid(row=1, column=1, padx=5)

    Button(button_frame, text="Calculate Overlap", command=lambda: calculate_overlaps(), bg="lightblue").grid(row=2, columnspan=2, pady=10)

    result_text = Text(root, width=60, height=20)
    result_text.pack(fill="both", expand=True)

    def calculate_overlaps():
        result_text.config(state="normal")
        result_text.delete(1.0, END)

        reads = [read.strip() for read in reads_entry.get().split(',') if read.strip()]
        try:
            k = int(k_entry.get())
            if not reads or k < 1:
                raise ValueError("Invalid inputs.")

            overlaps = naive_overlap_map(reads, k)
            if overlaps:
                for (a, b), olen in overlaps.items():
                    result_text.insert(END, f"Overlap between '{a}' and '{b}': {olen}\n")
            else:
                result_text.insert(END, "No overlaps found.\n")
        except ValueError:
            result_text.insert(END, "Please enter valid reads and a minimum overlap length (integer).\n")

        result_text.config(state="disabled")

    # Back Button to return to the main menu
    Button(button_frame, text="Back", command=menu_window, bg="lightgray").grid(row=3, columnspan=2, pady=10)




# Functions to display the input frame
def display_input_frame(kind_of_sequence):
    for widget in root.winfo_children():
        widget.destroy()

    global file_path_entry, extract_data_entry, status_label, display_frame , extract_K_entry

    input_frame = Frame(root, padx=10, pady=10)
    input_frame.pack(fill="x", pady=5)
    if kind_of_sequence != 9:
        Label(input_frame, text="File Path: ").grid(row=0, column=0, sticky="w")
        file_path_entry = Entry(input_frame, width=50)
        file_path_entry.grid(row=0, column=1, padx=5)
        Button(input_frame, text="Browse", command=browse_file).grid(row=0, column=2, padx=5)
            
    
    #Label of output file or patern

    if kind_of_sequence ==1 or kind_of_sequence ==2: 
        Label(input_frame, text="Output File Name: ").grid(row=1, column=0, sticky="w")
        extract_data_entry = Entry(input_frame, width=50)
        extract_data_entry.grid(row=1, column=1, padx=5)
    elif kind_of_sequence ==7:
        Label(input_frame, text="Patern: ").grid(row=1, column=0, sticky="w")
        extract_data_entry = Entry(input_frame, width=50)
        extract_data_entry.grid(row=1, column=1, padx=5)
    elif kind_of_sequence ==8:
        Label(input_frame, text="Patern: ").grid(row=1, column=0, sticky="w")
        extract_data_entry = Entry(input_frame, width=50)
        extract_data_entry.grid(row=1, column=1, padx=5)

        Label(input_frame, text="K: ").grid(row=2, column=0, sticky="w")
        extract_K_entry = Entry(input_frame, width=50)
        extract_K_entry.grid(row=2, column=1, padx=5)
    elif kind_of_sequence ==9:
        Label(input_frame, text="Text: ").grid(row=1, column=0, sticky="w")
        extract_data_entry = Entry(input_frame, width=50)
        extract_data_entry.grid(row=1, column=1, padx=5)
    

    
    #run of process 
    if kind_of_sequence ==1 or kind_of_sequence ==2: 
        Button(input_frame, text="Process", command=lambda: execute_and_show(kind_of_sequence), bg="lightblue").grid(row=2, columnspan=3, pady=10)
    else:
        Button(input_frame, text="Process", command=lambda:display_without_CSV(kind_of_sequence), bg="lightblue").grid(row=3, columnspan=3, pady=10)

    
    # Back Button to return to the main menu
    Button(input_frame, text="Back", command=menu_window, bg="lightgray").grid(row=4, columnspan=3, pady=10)

    status_label = Label(input_frame, text="")
    status_label.grid(row=5, columnspan=3)

    display_frame = Frame(root, padx=10, pady=10)
    display_frame.pack(fill="both", expand=True)
#********************Section3**********************
def display_without_CSV(kind_of_sequence):
    
    file_data = "Sections\sec5\dna1.fasta"
    if kind_of_sequence != 9:
        if kind_of_sequence ==8:    
            k= int(extract_K_entry.get())
        file_data = file_path_entry.get()
    if kind_of_sequence >=7: 
        patern= extract_data_entry.get()
    



    if not file_data:
        status_label.config(text="Please provide all inputs!", fg="red")
        return

    try:
        
        # Read the file content
        with open(file_data, 'r') as infile:
            lines = infile.readlines()
        seq = ''.join(line.strip() for line in lines if not line.startswith('>'))

        # Process based on the selected kind of sequence
        global result
        result=[("Original Sequence", seq)]

        if kind_of_sequence == 3:
            gc_content = GC_Content(seq)
            result += [("GC Content", f"{gc_content:.2%}")]

        elif kind_of_sequence == 4:
            reversed_seq = Reverse(seq)
            result += [("Reversed Sequence", reversed_seq)]

        elif kind_of_sequence == 5:
            reverse_complement = Reverse_Complement(seq)
            result += [("Reverse Complement", reverse_complement)]

        elif kind_of_sequence == 6:
            translation = Translation_Table(seq)
            result += [("Translated Sequence", translation)]
        elif kind_of_sequence == 7:
            Match = match(seq,patern)
            result += [("Matching Number", Match)]
            badchar=Badchars(seq,patern)
            result += [("BadChar Number", badchar)]
        elif kind_of_sequence == 8:
            index=IndexSorted(seq,k)
            p=patern
            q=query(seq,p,index)
            result += [("Query Indexing", q)]
        elif kind_of_sequence == 9:
            result=[("Original Sequence", patern)]
            suffix_output = display_suffix(patern)
            for idx, row in enumerate(suffix_output):
                result.append((f"Suffix Array Row {idx+1}", row))
           
        

        status_label.config(text="File processed successfully!", fg="green")

        # Clear previous display and create Treeview for results
        for widget in display_frame.winfo_children():
            widget.destroy()

        tree = ttk.Treeview(display_frame, columns=["Key", "Value"], show='headings', height=10)
        tree.heading("Key", text="Key")
        tree.heading("Value", text="Value")
        tree.column("Key", width=200, anchor='center')
        tree.column("Value", width=500, anchor='center')

        for key, value in result:
            tree.insert("", "end", values=(key, value))

        tree.pack(fill="both", expand=True)

    except Exception as e:
        status_label.config(text=f"Error: {str(e)}", fg="red")

#GC content
def GC_Content(seq):
    l = len(seq)
    num_G = seq.count("G")
    num_C = seq.count("C")
    total = num_C + num_G
    return total / l if l > 0 else 0  # Avoid division by zero
#Complement
def Complement(seq):
    dic = {"G": "C", "C": "G", "A": "T", "T": "A"}
    return ''.join(dic[base] for base in seq)
#Reverse
def Reverse(seq):
    return seq[::-1]

#Reverse Complement
def Reverse_Complement(seq):
    return Complement(Reverse(seq))

#Translation Table
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
            s += "#"  # Handle unknown codons
    return s


#****************Section 4 ***************************************
#***************Bad Character Matching****************************

def match(seq,sub_seq):
    x=-1
    for i in range(len(seq)):
        if sub_seq==seq[i:i+len(sub_seq)]:
            x=i
            #break
    return x

def Badchars(seq,sub_seq):
    table=np.zeros([4,len(sub_seq)])     
    row=["A","C","G","T"]
    for i in range (4):
        num=-1
        for j in range (len(sub_seq)):
            if row[i]==sub_seq[j]:
                table[i,j]=-1
                num=-1
            else:
                num+=1
                table[i,j]=num
    x=-1
    i=0
    while(i<len(seq)-len(sub_seq)+1):
        if sub_seq==seq[i:i+len(sub_seq)]:
            x=i
            #break
        
        else:
            for j in range(len(sub_seq)-1,-1,-1):
                if seq[i+j] != sub_seq[j]:
                    k=row.index(seq[i+j])
                    i+=table[k,j]
                    break
        i=int(i+1)
    return x

#********************Section5***************************************
#*******************Indexing****************************************
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

#********************Sectio6***************************************
#*******************Suffix Array****************************************
def display_suffix(x):
    # Dictionary to represent rank encoding of characters
    dec = {
        '$': 0,
        'A': 1,
        'C': 2,
        'G': 3,
        'T': 4
    }

    T = x
    table = []  # Table will hold rows of rank information
    i = 2 ** 0  # Starting index size (1 character)
    n = 0  # Iteration counter
    output = []  # Output list to hold each individual row
    
    # Initialize the table with the ranks of each character
    row = [dec[T[j]] for j in range(len(T))]
    table.append(row)
    #output.append(f"Row {n}: {row}")  # Add first row to the output

    while True:
        l = []  # List to store unique substrings
        dec2 = {}  # Dictionary to store the new ranks
        
        # If i > 1, we need to process substrings of length i
        if i > 1:
            # Create list of unique substrings
            for j in range(len(T)):
                if table[n-1][j:j+i] not in l:
                    l.append(table[n-1][j:j+i])
            
            # Sort substrings
            l.sort()
            
            # Create new ranks for each substring
            for j in range(len(l)):
                dec2[tuple(l[j])] = j
        
        # Create the new row based on the previous ranks
        row = []
        for j in range(len(T)):
            if i == 1:
                row.append(dec[T[j]])  # Single characters
            else:
                row.append(dec2[tuple(table[n-1][j:j+i])])  # Rank based on previous table
        
        # Add the row to the table and the output list
        table.append(row)
        output.append(f"Row {n+1}: {row}")  # Add the row with its index

        # Check if all ranks are unique in this row
        flag = 0
        for j in range(len(row)):
            c = row.count(j)
            if c > 1:
                flag = 1
                break

        # If all ranks are unique, stop the loop
        if flag == 0:
            break
        
        # Prepare for the next iteration (doubling the length of the substrings)
        n += 1
        i = 2 ** n
    
    return output
#********************Section7***************************************
#*******************Overlap****************************************

def naive_overlap(a, b, min_length):
    start = 0  
    while True:
        start = a.find(b[:min_length], start)  # Find b[:min_length] in a
        if start == -1:  
            return 0
        # Check if the suffix of 'a' from 'start' matches the prefix of 'b'
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def naive_overlap_map(reads, k):
    overlaps = {}
    for a in reads:
        for b in reads:
            if a != b:  # Avoid self-overlaps
                olen = naive_overlap(a, b, k)
                if olen > 0:
                    overlaps[(a, b)] = olen
    return overlaps


#/*******************End of Functions************************************/
def play_video_exit():
    video_path = r"Karim_end.mp4"  
    cap = cv2.VideoCapture(video_path)
    
    # Create a pop-up window
    video_window = Toplevel(root)
    center_screen(video_window)
    video_window.title("Thank you for using our application")
    change_ico(video_window)
    video_label = Label(video_window)
    video_label.pack(fill="both", expand=True)

    def update_frame():
        ret, frame = cap.read()
        if ret:
            # Convert frame to RGB
            frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
            # Convert frame to ImageTk format
            frame_image = ImageTk.PhotoImage(image=Image.fromarray(frame))
            video_label.config(image=frame_image)
            video_label.image = frame_image
            video_window.after(15, update_frame)  # Adjust to your video FPS
        else:
            cap.release()
            video_window.destroy()  # Close the pop-up window after video ends
            messagebox.showinfo("I hope project is useful", "Let us see you in the next project")
            root.quit()  # Exit the application

    # Start the video
    update_frame()
def exit_after_video():
    root.after(5900, root.quit) 
def menu_window():
    for widget in root.winfo_children():
        widget.destroy()

    root.title("Bio Methods")
    # Create a frame for the buttons
    button_frame = Frame(root, bg="#fff", padx=10, pady=10)
    button_frame.pack(fill="x", pady=5)

    # Buttons inside the frame
    #Section 2
    Button(button_frame, text="To_CSV_1", command=lambda: display_input_frame(1), font=Font, bg=Bg,fg=Fg,
           width=ButtonWidth, height=ButtonHeight).grid(row=1, column=1, padx=Padx, pady=Pady)

    Button(button_frame, text="To_CSV_2", command=lambda: display_input_frame(2), font=Font, bg=Bg,fg=Fg,
           width=ButtonWidth, height=ButtonHeight).grid(row=1, column=3, padx=Padx, pady=Pady)
    #Section 3
    Button(button_frame, text="GC_Content", command=lambda: display_input_frame(3), font=Font, bg=Bg,fg=Fg,
           width=ButtonWidth, height=ButtonHeight).grid(row=2, column=3, padx=Padx, pady=Pady)

    Button(button_frame, text="Reverse", command=lambda: display_input_frame(4), font=Font, bg=Bg,fg=Fg,
           width=ButtonWidth, height=ButtonHeight).grid(row=2, column=1, padx=Padx, pady=Pady)

    Button(button_frame, text="Reverse_Complement", command=lambda: display_input_frame(5), font=Font, bg=Bg,fg=Fg,
           width=ButtonWidth, height=ButtonHeight).grid(row=3, column=1, padx=Padx, pady=Pady)

    Button(button_frame, text="Translation", command=lambda: display_input_frame(6), font=Font, bg=Bg,fg=Fg,
           width=ButtonWidth, height=ButtonHeight).grid(row=3, column=3, padx=Padx, pady=Pady)
    #Section 4
    Button(button_frame, text="Bad_Char", command=lambda: display_input_frame(7), font=Font, bg=Bg,fg=Fg,
           width=ButtonWidth, height=ButtonHeight).grid(row=4, column=1, padx=Padx, pady=Pady)
    #Section 5
    Button(button_frame, text="Indexing", command=lambda: display_input_frame(8), font=Font, bg=Bg,fg=Fg,
           width=ButtonWidth, height=ButtonHeight).grid(row=4, column=3, padx=Padx, pady=Pady)
    #Section 6
    Button(button_frame, text="Suffix", command=lambda: display_input_frame(9), font=Font, bg=Bg,fg=Fg,
           width=ButtonWidth, height=ButtonHeight).grid(row=5, column=1, padx=Padx, pady=Pady)
    
    #Section 7
    Button(button_frame, text="Overlap", command=lambda:display_overlab() , font=Font, bg=Bg,fg=Fg,
           width=ButtonWidth, height=ButtonHeight).grid(row=5, column=3, padx=Padx, pady=Pady)
    
    #Exit Button 
    #D:\Collage\Intrudction to Bio informatics\Project\karim_end.mp4
    Button(button_frame, text="Exit", command=play_video_exit, font=Font, bg='red',
           width=ButtonWidth, height=ButtonHeight).grid(row=8, columnspan=6, padx=Padx, pady=Pady)
    play_video()


# Center screen and start video
def start_frame():
    center_screen(root)
    start_button = Button(root, text="Start", command=menu_window, font=Font, bg="lightblue", fg="black")
    start_button.place(relx=0.5, rely=0.5, anchor="center")
    play_video()

start_frame()
root.mainloop()


