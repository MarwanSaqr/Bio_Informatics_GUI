import tkinter as tk

# Function to calculate and display the suffix array
def display_suffix(x, output_box, update_row):
    dec = {
        '$': 0,
        'A': 1,
        'C': 2,
        'G': 3,
        'T': 4
    }

    T = x
    table = []
    i = 2**0
    n = 0
    output = ""

    # Initialize the table with the ranks of each character
    row = [dec[T[j]] for j in range(len(T))]
    table.append(row)

    def display_next_row():
        nonlocal i, n, row, table, output
        
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
        
        output += f"Row after processing length {i}: \n{row}\n"
        output_box.delete(1.0, tk.END)  # Clear the text box
        output_box.insert(tk.END, output)  # Insert the updated output
        
        # If all ranks are unique, stop the loop
        if flag == 0:
            return
        
        # Prepare for the next iteration (doubling the length of the substrings)
        n += 1
        i = 2**n

        # Schedule the next iteration
        output_box.after(1000, display_next_row)  # 1000ms = 1 second

    # Start displaying rows by calling display_next_row initially
    display_next_row()

# Function to handle the button click event
def show_suffix():
    x = input_text.get()  # Get the input string from the entry widget
    # Start the display_suffix function to update output_box
    display_suffix(x, output_box, show_suffix)

# Create the main Tkinter window
root = tk.Tk()
root.title("Suffix Array Visualization")

# Create a label and input box for the user to input the string
input_label = tk.Label(root, text="Enter string:")
input_label.pack(pady=5)

input_text = tk.Entry(root, width=50)
input_text.pack(pady=5)

# Create a button to trigger the suffix array display
show_button = tk.Button(root, text="Show Suffix Array", command=show_suffix)
show_button.pack(pady=10)

# Create a text box to display the output
output_box = tk.Text(root, height=20, width=80)
output_box.pack(padx=10, pady=10)

# Start the Tkinter event loop
root.mainloop()
