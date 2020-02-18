# -*- coding: utf-8 -*-
'''
Created on Wed Feb  5 12:08:03 2020

@author: fmunuera
'''

import sys
sys.path.append('./libs')

import os.path

import tkinter as tk
from tkinter.filedialog import askopenfilename, asksaveasfile
import tkinter.ttk as ttk

# recipe-580746-1.py from
# http://code.activestate.com/recipes/
# 580746-t kinter-treeview-like-a-table-or-multicolumn-listb/
import recipe5807461 as table

import json
from decimal import Decimal
from json import encoder


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

def solarfield_insert_row(solarfield_table):
    solarfield_table.insert_row(['','','','','',''])

def solarfield_del_rows(solarfield_table):
    solarfield_table.delete_all_selected_rows()

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

_DIR = {'saved_configurations' : './saved_configurations/',
        'weather_files' : './weather_files/',
        'fluid_files' : './fluid_files',
        'hce_files' : './hce_files',
        'sca_files' : './sca_files'}

root = tk.Tk()
root.attributes('-fullscreen', True)
root.title("wm min/max")
w, h = root.winfo_screenwidth(), root.winfo_screenheight()
root.geometry("%dx%d+0+0" % (w, h))
menubar = tk.Menu(root)

cfg_settings = {"simulation" : {},
                "solar_plant": {},
                "site": {},
                "solar_system": {},
                "weather": {},
                "hce": {},
                "sca": {},
                "hot_hce": {},
                "cold_hce": {},
                "cycle" : {},
                "hce_model_settings": {},
                "hce_scattered_params": {},
                "sca_scattered_params": {}
                }

nb = ttk.Notebook(root)
f1 = tk.Frame(nb, width=100, height=100)
f2 = tk.Frame(nb)
f3 = tk.Frame(nb)
f4 = tk.Frame(nb)
f5 = tk.Frame(nb)

nb.add(f1, text='  Solar Field Layout ', padding=2)
nb.add(f2, text='   Site & Weather    ', padding=2)
nb.add(f3, text='         HTF         ', padding=2)
nb.add(f4, text='HCE: Heat Collector Element', padding=2)
nb.add(f5, text='SCA: Solar Collector Assembly', padding=2)
nb.select(f1)
nb.enable_traversal()
nb.pack()


#  Solar Field Layout contruction Tab

def solarfield_load_dialog(f1, solarfield_table, ennamesolarfield,
                           title, labeltext = '' ):

    path = askopenfilename(initialdir = _DIR['saved_configurations'],
                           title = "choose your file",
                           filetypes = [("JSON files", "*.json")])

    with open(path) as cfg_file:
        cfg = json.load(cfg_file,
                                 parse_float= float,
                                 parse_int= int)

    datarow = []
    for r in cfg['solar_plant']['solarfields']:
        datarow.append(list(r.values()))

    solarfield_table.table_data = datarow
    solarfield_table.grid(row = 1, column = 0, columnspan =4)
    ennamesolarfield.delete(0, tk.END)
    ennamesolarfield.insert(0, cfg['solar_plant']['name'])
    f1.update()

def solarfield_save_dialog(f1, solarfield_table, ennamesolarfield, title, labeltext = '' ):

    #encoder.FLOAT_REPR = lambda o: format(o, '.2f')
    f = asksaveasfile(initialdir = _DIR['saved_configurations'],
                           title = "choose your file name",
                           filetypes = [("JSON files", "*.json")],
                           defaultextension = "json")

    cfg = dict({"solar_plant" : {}})
    cfg["solar_plant"].update(dict({"name" : ennamesolarfield.get()}))

    datarow = list(solarfield_table.table_data)
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

solarfield_table = table.Tk_Table(
                f1,
                ["NAME", "LOOPS", "SCAS", "HCES","MASSFLOW","ROWSPACING"],
                row_numbers=True,
                stripped_rows = ("white","#f2f2f2"),
                select_mode = "none",
                cell_anchor="center",
                adjust_heading_to_content = True)

lbnamesolarfield = ttk.Label(f1, text= "Name")
lbnamesolarfield.grid(row = 0, column = 0)
ennamesolarfield = ttk.Entry(f1, text= ":::::::::")
ennamesolarfield.grid(row = 0, column = 1)

btloadcfgsolarfield = ttk.Button(
        f1,
        text= "Load config",
        command= lambda : solarfield_load_dialog(
                f1, solarfield_table, ennamesolarfield, "Load config file", labeltext="Solar Field Config"))
btloadcfgsolarfield.grid(row = 0, column = 2)

btsavecfgsolarfield = ttk.Button(
        f1,
        text= "Save config",
        command= lambda : solarfield_save_dialog(
            f1, solarfield_table, ennamesolarfield, "Save config file",
            labeltext="Solar Field Config"))
btsavecfgsolarfield.grid(row = 0, column = 3)

btnewrow = ttk.Button(f1, text= "Insert",
                      command = lambda : solarfield_insert_row(solarfield_table))
btnewrow.grid(row = 2, column = 2)
btdelrows = ttk.Button(f1, text="Delete",
                       command = lambda : solarfield_del_rows(solarfield_table))
btdelrows.grid(row = 2, column = 3)


#  Site & Weather contruction Tab

def weather_load_dialog(f2, title, e, labeltext = '', path = None):

    try:
        if path is None:

            path = askopenfilename(initialdir = _DIR['weather_files'],
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
    f = asksaveasfile(initialdir = _DIR['weather_files'],
                           title = "choose your file name",
                           filetypes = [("JSON files", "*.json")],
                           defaultextension = "json")

    f.write(json.dumps(cfg_settings['weather']))
    f.close()

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
enweatherfile = ttk.Entry(f2)
enweatherfile.insert(0, "Path to Weather File")
enweatherfile.grid(row = 1, column = 1, columnspan = 4)
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
def fluid_load_dialog(f3, fluid_table, tmaxentry, tminentry,
                      entnamefluid, title, labeltext = '' ):

    path = askopenfilename(initialdir = _DIR['fluid_files'],
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

    fluid_table.table_data = datarow
    fluid_table.grid(row = 1, column = 1, columnspan = 7)

    entmaxfluid.delete(0, tk.END)
    entmaxfluid.insert(0, cfg['hot_fluid']['tmax'])
    entminfluid.delete(0, tk.END)
    entminfluid.insert(0, cfg['hot_fluid']['tmin'])
    entnamefluid.delete(0, tk.END)
    entnamefluid.insert(0, cfg['hot_fluid']['name'])

    f3.update()

def fluid_save_dialog(f3, fluid_table, entmax, entmin,
                      ennamefluid, title, labeltext = '' ):

    #encoder.FLOAT_REPR = lambda o: format(o, '.2f')
    f = asksaveasfile(initialdir = _DIR['fluid_files'],
                           title = "choose your file name",
                           filetypes = [("JSON files", "*.json")],
                           defaultextension = "json")

    cfg = dict({"hot_fluid" : {}})
    cfg["hot_fluid"].update(dict({"name" : ennamefluid.get()}))

    datarow = list(fluid_table.table_data)

    for r in datarow:
        param_name = r[0]
        param_values = list(map(to_number, r[1:]))
        cfg["hot_fluid"].update(dict({param_name : param_values}))

    cfg['hot_fluid'].update({"tmax" : float(entmax.get())})
    cfg['hot_fluid'].update({"tmin" : float(entmin.get())})
    f.write(json.dumps(cfg))
    f.close()

fluid_table = table.Tk_Table(
                f3,
                ["Parameter [x]", " A x^0", "B x^1", "C x^2","D x^3",
                 "E x^4", "F x^5"],
                row_numbers=True,
                stripped_rows = ("white","#f2f2f2"),
                select_mode = "none",
                cell_anchor="center",
                adjust_heading_to_content = True)

lbnamefluid = ttk.Label(f3, text= "Name")
lbnamefluid.grid(row = 0, column = 0)
ennamefluid = ttk.Entry(f3, text= "Name")
ennamefluid.grid(row = 0, column = 1)
btloadcfgfluid = ttk.Button(
    f3, text= "Load config",
    command= lambda : fluid_load_dialog(
            f3, fluid_table, entmaxfluid, entminfluid, ennamefluid,
            "Load config file", labeltext="HTF Config"))
btloadcfgfluid.grid(row = 0, column = 2)

lbtminfluid = ttk.Label(f3, text = "Tmin["+chr(176)+"C]")
lbtminfluid.grid(row=0,column = 3)
entminfluid = ttk.Entry(f3, text= "Tmin")
entminfluid.grid(row = 0, column = 4)
lbtmaxfluid = ttk.Label(f3, text = "Tmax ["+chr(176)+"C]")
lbtmaxfluid.grid(row=0,column = 5)
entmaxfluid = ttk.Entry(f3, text= "Tmax")
entmaxfluid.grid(row = 0, column = 6)

btsavecfgfluid = ttk.Button(
        f3,
        text= "Save config",
        command= lambda : fluid_save_dialog(
                f3, fluid_table, entmaxfluid, entminfluid, ennamefluid,
                "Save config file", labeltext="HTF Config"))
btsavecfgfluid.grid(row = 0, column = 7)


# HCE Construction tab

def hce_load_dialog(f4, hce_table, title, labeltext = '' ):

    path = askopenfilename(initialdir = _DIR['hce_files'],
                               title = "choose your file",
                               filetypes = [("JSON files", "*.json")])
    with open(path) as cfg_file:
        cfg = json.load(cfg_file, parse_float= float, parse_int= int)

    data = list(cfg)
    datarow = []
    for r in data:
        datarow.append(list(r.values()))

    hce_table.table_data = datarow


def hce_save_dialog(f4, hce_table, title, labeltext = '' ):

    #encoder.FLOAT_REPR = lambda o: format(o, '.2f')
    file = asksaveasfile(initialdir = _DIR['hce_files'],
                           title = "choose your file name",
                           filetypes = [("JSON files", "*.json")],
                           defaultextension = "json")
    cfg = []
    data = hce_table.table_data
    param_name = ["Name", "Description", "Condition", "Broken", "Bellows",
                 "Transmissivity", "Absorption", "Unaccounted",
                 "A0", "A1", "A2", "A3", "A4", "A5", "A6", "Factor"]
    for r in data:
        param_values = list(map(to_number, r[0:]))
        cfg.append(dict(zip(param_name, param_values)))

    file.write(json.dumps(cfg))
    file.close()

hce_table = table.Tk_Table(
                f4,
                ["Name", "Description", "Condition", "Broken", "Bellows",
                 "Transmissivity", "Absorption", "Unaccounted",
                 "A0", "A1", "A2", "A3", "A4", "A5", "A6", "Factor"],
                row_numbers=True,
                stripped_rows = ("white","#f2f2f2"),
                select_mode = "none",
                cell_anchor="center",
                adjust_heading_to_content = True)

hce_table.grid(row = 0, columnspan = 2)
btloadcfghce = ttk.Button(
    f4, text= "Load config",
    command= lambda : hce_load_dialog(
            f4, hce_table,
            "Load config file", labeltext="HCE Config"))
btloadcfghce.grid(row = 1, column = 0)
btsavecfghce = ttk.Button(
        f4,
        text= "Save config",
        command= lambda : hce_save_dialog(
                f4, hce_table,
                "Save config file", labeltext="HCE Config"))
btsavecfghce.grid(row = 1, column = 1)
f4.update()



# SCA Construction tab

def sca_load_dialog(f5, sca_table, title, labeltext = '' ):

    path = askopenfilename(initialdir = _DIR['sca_files'],
                           title = "choose your file",
                           filetypes = [("JSON files", "*.json")])

    with open(path) as cfg_file:
        cfg = json.load(cfg_file, parse_float= float, parse_int= int)

    data = list(cfg)
    datarow = []
    for r in data:
        datarow.append(list(r.values()))

    sca_table.table_data = datarow
    sca_table.grid(row = 1, column = 0)
    f5.update()

def sca_save_dialog(f5, sca_table, title, labeltext = '' ):

    #encoder.FLOAT_REPR = lambda o: format(o, '.2f')
    file = asksaveasfile(initialdir = _DIR['sca_files'],
                           title = "choose your file name",
                           filetypes = [("JSON files", "*.json")],
                           defaultextension = "json")
    cfg = []
    data = sca_table.table_data
    param_name = ["Name", "Description", "Condition", "Broken", "Bellows",
                 "Transmissivity", "Absorption", "Unaccounted",
                 "A0", "A1", "A2", "A3", "A4", "A5", "A6", "Factor"]
    for r in data:
        param_values = list(map(to_number, r[0:]))
        cfg.append(dict(zip(param_name, param_values)))

    file.write(json.dumps(cfg))
    file.close()

sca_table = table.Tk_Table(
                f5, ["Name","SCA Length","Aperture","Aperture Area","Focal Len",
                    "IAM Coefficient F0","IAM Coefficient F1", "AM Coefficient F2",
                    "Track  Twist", "Geom.Accuracy", "Reflectance", "Cleanliness",
                    "Dust", "Factor", "Availability"],
                row_numbers=True,
                stripped_rows = ("white","#f2f2f2"),
                select_mode = "none",
                cell_anchor="center",
                adjust_heading_to_content = True)

#lbnamesca = ttk.Label(f5, text= "Name")
#lbnamesca.grid(row = 0, column = 0)
#ennamesca = ttk.Entry(f5, text= "Name")
#ennamesca.grid(row = 0, column = 1)
#
#entmaxsca = ttk.Entry(f5, text= "Tmax")
#entminsca = ttk.Entry(f5, text= "Tmin")
#entmaxsca.grid(row = 0, column = 3)
#entminsca.grid(row = 0, column = 4)

btloadcfgsca = ttk.Button(
    f5, text= "Load config",
    command= lambda : sca_load_dialog(
            f5, sca_table,
            "Load config file", labeltext="SCA Config"))
btloadcfgsca.grid(row = 0, column = 2)

btsavecfgsca = ttk.Button(
        f5,
        text= "Save config",
        command= lambda : sca_save_dialog(
                f5, sca_table,
                "Save config file", labeltext="HTF Config"))
btsavecfgsca.grid(row = 0, column = 5)





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