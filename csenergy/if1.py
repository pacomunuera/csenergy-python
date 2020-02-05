# -*- coding: utf-8 -*-
'''
Created on Wed Feb  5 12:08:03 2020

@author: fmunuera
'''

import tkinter as tk
from tkinter.filedialog import askopenfilename, asksaveasfile
from tkinter.messagebox import *
import tkinter.ttk as ttk

def simulation_new():
    pass

def simulation_open():
    pass

def simulation_save():
    name=asksaveasfile(mode='w',defaultextension='.json')
    text2save=str(text.get(0.0,END))
    name.write(text2save)
    name.close

def simulation_save_as():
    pass

def simulation_close():
    pass

def run_simulate():
    pass

def help_help():
    pass

def help_about():
    pass

def simulation_exit():
    root.destroy()


root = Tk()
root.geometry('800x600')
menubar=Menu(root)
#frame = Frame(root)

nb = ttk.Notebook(root)
f1 = tk.Frame(nb)
f2 = tk.Frame(nb)
f3 = tk.Frame(nb)

nb.add(f1, text='  Solar Field Layout ', padding=2)
nb.add(f2, text='   Site & Weather    ', padding=2)
nb.add(f3, text='         HTF         ', padding=2)
nb.select(f1)
nb.enable_traversal()
nb.pack()

#l1 = ttk.Label(f1, text='Subcampos').pack()
#t1 = ttk.Entry(f1, text='Numero de subcampos').pack()
#l2 = ttk.Label(f1, text='Lazos').pack()
#t2 = ttk.Entry(f1, text='Numero de Lazos').pack()
#text = ttk.Button(f1, text='boton', command= simulation_close).pack()
#lista = ttk.LabeledScale(f1, from_=0, to=10).pack()
#a = ttk.Combobox(f1).pack()
#b = ttk.Labelframe(f1).pack()
#c = ttk.Panedwindow(f1).pack()
#d = ttk.Scale(f1).pack()
#e = ttk.Scrollbar(f1).pack()
#f = ttk.Sizegrip(f1).pack()
#g = ttk.Spinbox(f1).pack()

config_table = ttk.Treeview(f1, columns=("NAME", "LOOPS", "SCAS", "HCES"))
config_table.pack()

config_table.column("NAME", width=150)
config_table.column("LOOPS", width=150)
config_table.column("SCAS", width=150)
config_table.column("HCES", width=150)
config_table.heading("NAME", text="Name")
config_table.heading("LOOPS", text="Loops")
config_table.heading("SCAS", text="SCAs")
config_table.heading("HCES", text="HCEs")


datacols = [('C1', 1, 4, 12),
           ('C2', 1, 4, 12),
           ('C3', 1, 4, 12),
           ('C4', 1, 4, 12)]

for item in datacols:
            config_table.insert('' ,'end', values=item)




simulation_menu=Menu(menubar,tearoff=0)
simulation_menu.add_command(label='New', command=simulation_new)
simulation_menu.add_command(label='Open', command=simulation_open)
simulation_menu.add_command(label='Save', command=simulation_save)
simulation_menu.add_command(label='Save as...', command=simulation_save_as)
simulation_menu.add_command(label='Close', command=simulation_close)
simulation_menu.add_separator()
simulation_menu.add_command(label='Exit', command=simulation_exit)
menubar.add_cascade(label='Simulation', menu=simulation_menu)

run_menu=Menu(menubar,tearoff=0)
run_menu.add_command(label='Run', command=run_simulate)
menubar.add_cascade(label='Run', menu=run_menu)

help_menu=Menu(menubar,tearoff=0)
help_menu.add_command(label='Help',command=help_help)
help_menu.add_command(label='About',command=help_about)
menubar.add_cascade(label='Help',menu=help_menu)
root.config(menu=menubar)
root.mainloop()