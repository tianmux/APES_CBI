import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import pandas as pd

def open_file(fn):
    filename = filedialog.askopenfilename(initialdir="/home/txin/Documents/Work_Folder/Code/APES/APES_CBI_recon/view/canvas_test/", title="Select file",
                                          filetypes=(("CSV files", "*.csv"), ("all files", "*.*")))
    df = pd.read_csv(filename)
    return df

def plot_data(df, column):
    global canvas
    canvas.get_tk_widget().destroy()
    print(df)
    fig = Figure(figsize=(6,5), dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(df[column])
    #df[column].plot(kind='line', legend=True, ax=ax)
    ax.set_title('Plot of ' + column)
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack()

root = tk.Tk()

button_open = tk.Button(root, text="Open File", command=open_file)
button_open.pack()

canvas = None
file_name = "/home/txin/Documents/Work_Folder/Code/APES/APES_CBI_recon/view/canvas_test/data.csv"
df = None

df = open_file(file_name)
print(df.columns)
fig = Figure(figsize=(6,5), dpi=100)
ax = fig.add_subplot(111)
df['data1'].plot(kind='line', legend=True, ax=ax)
ax.set_title('Plot of ' + 'data1')
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack()

canvas.get_tk_widget().destroy()
fig = Figure(figsize=(6,5), dpi=100)
ax = fig.add_subplot(111)
df['data2'].plot(kind='line', legend=True, ax=ax)
ax.set_title('Plot of ' + 'data2')
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack()


# Create two buttons
button1 = tk.Button(root, text="Plot Column1", command=lambda: plot_data(df, 'data1'))
button1.pack()
button2 = tk.Button(root, text="Plot Column2", command=lambda: plot_data(df, 'data2'))
button2.pack()

root.mainloop()



