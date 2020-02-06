# -*- coding: utf-8 -*-
'''
Created on Wed Feb  5 12:08:03 2020

@author: fmunuera
'''

import tkinter as tk
from tkinter.filedialog import askopenfilename, asksaveasfile
#  from tkinter.messagebox import *
import tkinter.ttk as ttk
import recipe5807461 as rcp1
import json
from decimal import Decimal
from json import encoder

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

def solarfield_insert_row(config_table):
    config_table.insert_row(['','','','','',''])

def solarfield_del_rows(config_table):
    config_table.delete_all_selected_rows()

def run_simulate():
    pass

def help_help():
    pass

def help_about():
    pass

def simulation_exit():
    root.destroy()


#class solarfield_dialog(object):
#
#    def __init__(self, f1, config_table, title, labeltext = '' ):
#
#        path = askopenfilename(initialdir = ".saved_configurations/",
#                               title = "choose your file",
#                               filetypes = [("JSON files", "*.json")])
#
#        with open(path) as cfg_file:
#            cfg_settings = json.load(cfg_file)
#
#        datarow = []
#        for r in cfg_settings['solar_plant']['solarfields']:
#            datarow.append(tuple(r.values()))
#
#        config_table.table_data = datarow
#        config_table.pack(expand=True, fill= 'both')
#        f1.update()

def to_number(s):

    try:
    	i = int(s)
    	return(i)
    except ValueError:
        pass
    try:
        f = float(s)
        return(f)
    except ValueError:
        pass
    return s
    	
    
def solarfield_load_dialog(f1, config_table, title, labeltext = '' ):

    path = askopenfilename(initialdir = ".saved_configurations/",
                           title = "choose your file",
                           filetypes = [("JSON files", "*.json")])

    with open(path) as cfg_file:
        cfg_settings = json.load(cfg_file, 
                                 parse_float= float,
                                 parse_int= int)

    datarow = []
    for r in cfg_settings['solar_plant']['solarfields']:
        datarow.append(list(r.values()))

    config_table.table_data = datarow
    config_table.pack(expand=True, fill= 'both')
    f1.update()

def solarfield_save_dialog(f1, config_table, name, title, labeltext = '' ):

    #encoder.FLOAT_REPR = lambda o: format(o, '.2f')
    f = asksaveasfile(initialdir = ".saved_configurations/",
                           title = "choose your file name",
                           filetypes = [("JSON files", "*.json")])

    cfg_settings = dict({"solar_plant" : {}})
    cfg_settings["solar_plant"].update(dict({"name" : name}))

    datarow = list(config_table.table_data)
    print(datarow)
    dictkeys =["name", "loops", "scas", "hces", "massflow", "row_spacing"]

    solarfields = []

    for r in datarow:
        sf = {}
        index = 0
        for v in r:
            k = dictkeys[index]
            sf[k]= to_number(v)
            index += 1
        solarfields.append(sf)
        
    cfg_settings['solar_plant'].update({"solarfields" : solarfields})
        
    print(cfg_settings)
    f.write(json.dumps(cfg_settings))
    f.close()



root = tk.Tk()
root.geometry('800x600')
menubar = tk.Menu(root)
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

#  Solar Field Layout contruction

config_table = rcp1.Tk_Table(
                f1,
                ["NAME", "LOOPS", "SCAS", "HCES","MASSFLOW","ROWSPACING"],
                row_numbers=True,
                stripped_rows = ("white","#f2f2f2"),
                select_mode = "none",
                cell_anchor="center")


enname = ttk.Entry(f1, text= "Name").pack()

btloadcfg = ttk.Button(
        f1,
        text= "Load config",
        command= lambda : solarfield_load_dialog(
                f1, config_table, "Load config file", labeltext="Solar Field Config"))
btloadcfg.pack()
btsavecfg = ttk.Button(
        f1,
        text= "Save config",
        command= lambda : solarfield_save_dialog(
                f1, config_table, "nombre",  "Save config file", labeltext="Solar Field Config"))
btsavecfg.pack()





#datarows = [('C1', 1, 4, 12),
#            ('C2', 1, 4, 12),
#            ('C3', 1, 4, 12),
#            ('C4', 1, 4, 12)]
#
#config_table.table_data = datarows

btnewrow = ttk.Button(f1, text= "Insert",
                      command = lambda : solarfield_insert_row(config_table))

btdelrows = ttk.Button(f1, text="Delete",
                       command = lambda : solarfield_del_rows(config_table))

btnewrow.pack()
btdelrows.pack()



#entry_name = ttk.Entry(f1).pack()
#entry_loops= ttk.Entry(f1).pack()
#entre_scas = ttk.Entry(f1).pack()
#entry_hces = ttk.Entry(f1).pack()


#config_table = ttk.Treeview(f1,
#                            columns=("NAME", "LOOPS", "SCAS", "HCES"),
#                            selectmode = "extended", show ="headings")
#config_table.pack()
#config_table.heading("NAME", text="Name")
#config_table.heading("LOOPS", text="Loops")
#config_table.heading("SCAS", text="SCAs")
#config_table.heading("HCES", text="HCEs")
#
#config_table.column("NAME", width=150)
#config_table.column("LOOPS", width=150)
#config_table.column("SCAS", width=150)
#config_table.column("HCES", width=150)
#
#config_table.grid(row=0, column=0, columnspan=3, padx=10, pady=10,
#                  sticky = "nsew")
#for item in datacols:
#            config_table.insert('' ,'end', values=item)


simulation_menu= tk.Menu(menubar,tearoff=0)
simulation_menu.add_command(label='New', command=simulation_new)
simulation_menu.add_command(label='Open', command=simulation_open)
simulation_menu.add_command(label='Save', command=simulation_save)
simulation_menu.add_command(label='Save as...', command=simulation_save_as)
simulation_menu.add_command(label='Close', command=simulation_close)
simulation_menu.add_separator()
simulation_menu.add_command(label='Exit', command=simulation_exit)
menubar.add_cascade(label='Simulation', menu=simulation_menu)

run_menu= tk.Menu(menubar,tearoff=0)
run_menu.add_command(label='Run', command=run_simulate)
menubar.add_cascade(label='Run', menu=run_menu)

help_menu= tk.Menu(menubar,tearoff=0)
help_menu.add_command(label='Help',command=help_help)
help_menu.add_command(label='About',command=help_about)
menubar.add_cascade(label='Help',menu=help_menu)
root.config(menu=menubar)
root.mainloop()