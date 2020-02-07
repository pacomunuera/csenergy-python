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
import os.path

def simulation_new():
    pass

def simulation_open():
    pass

def simulation_save(text):
    name=asksaveasfile(mode='w',defaultextension='.json')
    text2save=str(text.get(0.0,tk.END))
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
        cfg = json.load(cfg_file,
                                 parse_float= float,
                                 parse_int= int)

    datarow = []
    for r in cfg['solar_plant']['solarfields']:
        datarow.append(list(r.values()))

    config_table.table_data = datarow

    config_table.grid(row = 1, columnspan =4)
    f1.update()

def solarfield_save_dialog(f1, config_table, name, title, labeltext = '' ):

    #encoder.FLOAT_REPR = lambda o: format(o, '.2f')
    f = asksaveasfile(initialdir = ".saved_configurations/",
                           title = "choose your file name",
                           filetypes = [("JSON files", "*.json")],
                           defaultextension = "json")

    cfg = dict({"solar_plant" : {}})
    cfg["solar_plant"].update(dict({"name" : name}))

    datarow = list(config_table.table_data)
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

    cfg['solar_plant'].update({"solarfields" : solarfields})
    cfg_settings['solar_plant'].update(dict(cfg))
    f.write(json.dumps(cfg))
    f.close()

def weather_load_dialog(f2, title, e, labeltext = '', path = None):

    try:
        if path is None:

            path = askopenfilename(initialdir = ".weather_files/",
                               title = "choose your file",
                               filetypes = (("TMY files","*.tm2"),
                                            ("TMY files","*.tm3"),
                                            ("csv files","*.csv"),
                                            ("all files","*.*")))

            if path is None:
                return
            else:
                filedir = os.path.dirname(path)+"/"
                filename = os.path.basename(path)
                cfg = dict({"weather": {}})
                cfg["weather"].update(
                        {"filepath": filedir,
                         "filename": filename})
                cfg_settings["weather"].update(dict(cfg))
                e.delete(0, tk.END)
                e.insert(0, path)

    except Exception:
        raise
        tk.txMessageBox.showerror('Error loading Weather Data File',
                               'Unable to open file: %r', path)
        f2.update()

def weather_save_dialog(f2, name, title, labeltext = '' ):

    #encoder.FLOAT_REPR = lambda o: format(o, '.2f')
    f = asksaveasfile(initialdir = ".saved_configurations/",
                           title = "choose your file name",
                           filetypes = [("JSON files", "*.json")],
                           defaultextension = "json")

    f.write(json.dumps(cfg_settings['weather']))
    f.close()


root = tk.Tk()
root.geometry('800x600')
menubar = tk.Menu(root)

cfg_settings = {"simulation" : {},
                "solar_plant": {},
                "site": {},
                "solar_system": {},
                "weather": {},
                "hce": {},
                "sca": {},
                "hot_fluid": {},
                "cold_fluid": {},
                "cycle" : {},
                "hce_model_settings": {},
                "hce_scattered_params": {},
                "sca_scattered_params": {}
                }

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


#  Solar Field Layout contruction Tab

config_table = rcp1.Tk_Table(
                f1,
                ["NAME", "LOOPS", "SCAS", "HCES","MASSFLOW","ROWSPACING"],
                row_numbers=True,
                stripped_rows = ("white","#f2f2f2"),
                select_mode = "none",
                cell_anchor="center",
                adjust_heading_to_content = True)

lbname = ttk.Label(f1, text= "Name").grid(row = 0, column = 0)
enname = ttk.Entry(f1, text= "Name").grid(row = 0, column = 1)

btloadcfgsolarfield = ttk.Button(
        f1,
        text= "Load config",
        command= lambda : solarfield_load_dialog(
                f1, config_table, "Load config file", labeltext="Solar Field Config"))
btloadcfgsolarfield.grid(row = 0, column = 2)
btsavecfgsolarfield = ttk.Button(
        f1,
        text= "Save config",
        command= lambda : solarfield_save_dialog(
                f1, config_table, "nombre",  "Save config file", labeltext="Solar Field Config"))
btsavecfgsolarfield.grid(row = 0, column = 3)

btnewrow = ttk.Button(f1, text= "Insert",
                      command = lambda : solarfield_insert_row(config_table))

btdelrows = ttk.Button(f1, text="Delete",
                       command = lambda : solarfield_del_rows(config_table))

btnewrow.grid(row = 2, column = 0)
btdelrows.grid(row = 2, column = 1)


#  Site & Weather contruction Tab

lbsitename = ttk.Label(f2, text= "Name").grid(row = 0, column = 0)
ensitename = ttk.Entry(f2, text= "Site").grid(row = 0, column = 1)
lbsitelat = ttk.Label(f2, text= "Latitude").grid(row = 0, column = 2)
ensitelat = ttk.Entry(f2, text= "Latitude").grid(row = 0, column = 3)
lbsitelong = ttk.Label(f2, text= "Longitude").grid(row = 0, column = 4)
ensitelong = ttk.Entry(f2, text= "Longitude").grid(row = 0, column = 5)
lbsitealt = ttk.Label(f2, text= "Altitude").grid(row = 0, column = 6)
ensitealt = ttk.Entry(f2, text= "Altitude").grid(row = 0, column = 7)

lbweatherfile = ttk.Label(f2, text= "Weather File")
lbweatherfile.grid(row = 1, column = 0)
enweatherfile = ttk.Entry(f2, width=75)
enweatherfile.insert(0, "Path to Weather File")
enweatherfile.grid(row = 1, column = 1, columnspan = 5)

btloadcfgweather = ttk.Button(
        f2,
        text= "Load Weather",
        command= lambda : weather_load_dialog(
                f2, "Load config file", enweatherfile, labeltext="Wather File"))
btloadcfgweather.grid(row = 1, column = 6)
btsavecfgweather = ttk.Button(
        f2,
        text= "Save config",
        command= lambda : weather_save_dialog(
                f2, "nombre",  "Save config file", labeltext="Weather"))
btsavecfgweather.grid(row = 1, column = 7)

#  Fluid contruction Tab

def fluid_load_dialog(f3, config_fluid_table, tmaxentry, tminentry, 
                      entnamefluid, title, labeltext = '' ):

    path = askopenfilename(initialdir = ".saved_configurations/",
                           title = "choose your file",
                           filetypes = [("JSON files", "*.json")])

    with open(path) as cfg_file:
        cfg = json.load(cfg_file, parse_float= float, parse_int= int)

    cp_coefs = [[]]*7
    rho_coefs = [[]]*7
    mu_coefs  = [[]]*7
    kt_coefs  = [[]]*7

    temp_cp = ["cp"]
    temp_rho = ["rho"]
    temp_mu = ["mu"]
    temp_kt = ["kt"]

    grades = ["Factor"]
    grades.extend([0, 1, 2, 3, 4, 5])

    temp_cp.extend(list(cfg['hot_fluid']['cp']))
    temp_rho.extend(list(cfg['hot_fluid']['rho']))
    temp_mu.extend(list(cfg['hot_fluid']['mu']))
    temp_kt.extend(list(cfg['hot_fluid']['kt']))

    for index in range(len(temp_cp)):
        cp_coefs[index] = temp_cp[index]
    for index in range(len(temp_rho)):
        rho_coefs[index] = temp_rho[index]
    for index in range(len(temp_mu)):
        mu_coefs[index] = temp_mu[index]
    for index in range(len(temp_kt)):
        kt_coefs[index] = temp_kt[index]

    datarow = []
    #datarow.append(grades)
    datarow.append(cp_coefs)
    datarow.append(rho_coefs)
    datarow.append(mu_coefs)
    datarow.append(kt_coefs)

    config_fluid_table.table_data = datarow
    config_fluid_table.grid(row = 1, column = 1)
    
    entmaxfluid.delete(0, tk.END)
    entmaxfluid.insert(0, cfg['hot_fluid']['tmax'])
    entminfluid.delete(0, tk.END)
    entminfluid.insert(0, cfg['hot_fluid']['tmin'])
    entnamefluid.delete(0, tk.END)
    entnamefluid.insert(0, cfg['hot_fluid']['name'])
    
    f3.update()

def fluid_save_dialog(f3, config_table, entmax, entmin,
                      name, title, labeltext = '' ):

    #encoder.FLOAT_REPR = lambda o: format(o, '.2f')
    f = asksaveasfile(initialdir = ".saved_configurations/",
                           title = "choose your file name",
                           filetypes = [("JSON files", "*.json")],
                           defaultextension = "json")

    cfg = dict({"hot_fluid" : {}})
    cfg["hot_fluid"].update(dict({"name" : name}))

    datarow = list(config_table.table_data)

    for r in datarow:
        param_name = r[0]
        param_values = list(map(to_number, r[1:]))
        cfg["hot_fluid"].update(dict({param_name : param_values}))

    cfg['hot_fluid'].update({"tmax" : float(entmax.get())})
    cfg['hot_fluid'].update({"tmin" : float(entmin.get())})
    f.write(json.dumps(cfg))
    f.close()


config_fluid_table = rcp1.Tk_Table(
                f3,
                ["Parameter", " Grade 0", "Grade 1", "Grade 2","Grade 3",
                 "Grade 4", "Grade 5"],
                row_numbers=True,
                stripped_rows = ("white","#f2f2f2"),
                select_mode = "none",
                cell_anchor="center",
                adjust_heading_to_content = True)

lbnamefluid = ttk.Label(f3, text= "Name")
lbnamefluid.grid(row = 0, column = 0)
ennamefluid = ttk.Entry(f3, text= "Name")
ennamefluid.grid(row = 0, column = 1)

entmaxfluid = ttk.Entry(f3, text= "Tmax")
entminfluid = ttk.Entry(f3, text= "Tmin")
entmaxfluid.grid(row = 0, column = 3)
entminfluid.grid(row = 0, column = 4)


btloadcfgfluid = ttk.Button( 
    f3, text= "Load config", 
    command= lambda : fluid_load_dialog(
            f3, config_fluid_table, entmaxfluid, entminfluid, ennamefluid,
            "Load config file", labeltext="HTF Config"))

btloadcfgfluid.grid(row = 0, column = 2)

btsavecfgfluid = ttk.Button(
        f3,
        text= "Save config",
        command= lambda : fluid_save_dialog(
                f3, config_fluid_table, entmaxfluid, entminfluid, "Name",  "Save config file", labeltext="HTF Config"))
btsavecfgfluid.grid(row = 0, column = 5)


#entry_name = ttk.Entry(f1).pack()
#entry_loops= ttk.Entry(f1).pack()
#entre_scas = ttk.Entry(f1).pack()
#entry_hces = ttk.Entry(f1).pack()


#config_table = ttk.Treeview(f1,
#                            columns=("NAME", "LOOPS", "SCAS", "HCES"),
#                            selectmode = "extended", show ="headings")

#config_table.configure_column(width = 100, anchor = 'center')
#config_table.heading("LOOPS", text="Loops")
#config_table.heading("SCAS", text="SCAs")
#config_table.heading("HCES", text="HCEs")
#
#config_table.column("NAME", width=150)
#config_table.column("LOOPS", width=150)
#config_table.column("SCAS", width=150)
#config_table.column("HCES", width=150)

#config_table.grid(row=0, column=0, columnspan=3, padx=10, pady=10,
#                  sticky = "nsew")



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