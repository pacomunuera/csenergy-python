# -*- coding: utf-8 -*-
'''
Created on Wed Feb  5 12:08:03 2020

@author: fmunuera
'''

import sys
sys.path.append('./libs')
import os.path

import pvlib as pvlib
from pvlib import iotools

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

import pandas as pd

class Interface(object):

    _COOLPROP_FLUIDS = ['Water', 'INCOMP::TVP1', 'INCOMP::S800']

    _DIR = {'saved_configurations' : './saved_configurations/',
        'site_files' : './site_files/',
        'fluid_files' : './fluid_files',
        'hce_files' : './hce_files',
        'sca_files' : './sca_files',
        'fielddata_files': './fielddata_files',
        'weather_files': './weather_files'}

    cfg_settings = {"simulation" : {},
                "solar_plant": {},
                "site": {},
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

    def __init__(self):

        # Main window
        self.root = tk.Tk()
        self.root.attributes('-fullscreen', True)
        self.root.title("wm min/max")
        w, h = self.root.winfo_screenwidth(), self.root.winfo_screenheight()
        self.root.geometry("%dx%d+0+0" % (w, h))

        # Menu
        self.menubar = tk.Menu(self.root)

        self.simulation_menu = tk.Menu(self.menubar,tearoff=0)
        self.simulation_menu.add_command(label='New', command=self.simulation_new)
        self.simulation_menu.add_command(label='Open', command=self.simulation_open)
        self.simulation_menu.add_command(label='Save', command=self.simulation_save)
        self.simulation_menu.add_command(label='Save as...', command=self.simulation_save_as)
        self.simulation_menu.add_command(label='Close', command=self.simulation_close)
        self.simulation_menu.add_separator()
        self.simulation_menu.add_command(label='Exit', command=self.simulation_exit)
        self.menubar.add_cascade(label='Simulation', menu=self.simulation_menu)

        self.run_menu= tk.Menu(self.menubar,tearoff=0)
        self.run_menu.add_command(label='Run', command=self.run_simulate)
        self.menubar.add_cascade(label='Run', menu=self.run_menu)

        self.help_menu= tk.Menu(self.menubar,tearoff=0)
        self.help_menu.add_command(label='Help',command=self.help_help)
        self.help_menu.add_command(label='About',command=self.help_about)
        self.menubar.add_cascade(label='Help',menu=self.help_menu)
        self.root.config(menu=self.menubar)


        # Notebook (tabs)
        self.nb = ttk.Notebook(self.root)

        self.fr_simulation = tk.Frame(self.nb)
        self.fr_solarfield = tk.Frame(self.nb)
        #self.fr_data = tk.Frame(self.nb)
        self.fr_fluid = tk.Frame(self.nb)
        self.fr_hce = tk.Frame(self.nb)
        self.fr_sca = tk.Frame(self.nb)


        self.buildNotebook()
        self.buildSolarFieldFrame()
        self.buildSimulationFrame()
        self.buildFluidFrame()
        #self.buildDataFrame()
        #self.buildHCEFrame()
        # self.buildSCAFrame()

    def simulation_new(self):
        pass

    def simulation_open(self):

        path = askopenfilename(initialdir = self._DIR['saved_configurations'],
                               title = "choose your file",
                               filetypes = [("JSON files", "*.json")])

        with open(path) as cfg_file:
            cfg = json.load(cfg_file,
                                     parse_float= float,
                                     parse_int= int)

        self.varsimID.set(cfg['simulation']['ID'])
        self.varsimdatatype.set(cfg['simulation']['datatype'])
        self.varsimulation.set(cfg['simulation']['benchmark'])
        self.varfastmode.set(cfg['simulation']['fastmode'])
        self.vardatafilename .set(cfg['simulation']['filename'])
        self.vardatafilepath.set(cfg['simulation']['filepath'])

        datarow = []
        for r in cfg['solarfield']['subfields']:
            datarow.append(list(r.values()))

        self.solarfield_table.table_data = datarow
        self.enname.delete(0, tk.END)
        self.varsolarfieldname.set(cfg['solarfield']['name'])
        self.vartin.set(cfg['solarfield']['rated_tin'])
        self.vartout.set(cfg['solarfield']['rated_tout'])
        self.varpin.set(cfg['solarfield']['rated_pin'])
        self.varpout.set(cfg['solarfield']['rated_pout'])
        self.vartmin.set(cfg['solarfield']['tmin'])
        self.vartmax.set(cfg['solarfield']['tmax'])
        self.varratedmassflow.set(cfg['solarfield']['rated_massflow'])
        self.varrecirculation.set(cfg['solarfield']['recirculation_massflow'])
        self.varscas.set(cfg['solarfield']['loop']['scas'])
        self.varhces.set(cfg['solarfield']['loop']['hces'])

        #self.enratedtin.insert(0, cfg['solarfield']['rated_tin'])
        # self.enratedtout.delete(0, tk.END)
        # self.enratedtout.insert(0, cfg['solarfield']['rated_tout'])

        self.fr_solarfield.update()

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

    def __insert_rows__(self, table):

        lst = []
        for c in table._multicolumn_listbox._columns:
            lst.append('')

        rows = table.number_of_rows
        table.insert_row(lst, index = rows)

    def __del_rows__(self, table):
        table.delete_all_selected_rows()

    def run_simulate():
        pass

    def help_help():
        pass

    def help_about():
        pass

    def simulation_exit(self):

        self.root.destroy()


    def to_number(self, s):

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


    def buildNotebook(self):

        self.nb.add(self.fr_simulation, text='Simulation Configuration', padding=2)
        self.nb.add(self.fr_solarfield, text='  Solar Field Layout ', padding=2)
        #self.nb.add(self.fr_data, text='   Site & Weather    ', padding=2)
        self.nb.add(self.fr_fluid, text='         HTF         ', padding=2)
        self.nb.add(self.fr_hce, text='HCE: Heat Collector Element', padding=2)
        self.nb.add(self.fr_sca, text='SCA: Solar Collector Assembly', padding=2)
        self.nb.select(self.fr_simulation)
        self.nb.enable_traversal()
        self.nb.pack()

    #  Simulaton Configuration tab

    def simulationLoadDialog(self, title, labeltext=''):

        path = askopenfilename(initialdir = self._DIR['saved_configurations'],
                               title = "choose your file",
                               filetypes = [("JSON files", "*.json")])

        with open(path) as cfg_file:
            cfg = json.load(cfg_file,
                                     parse_float= float,
                                     parse_int= int)
        # datarow = []
        # for r in cfg['simulation']:
        #     datarow.append(list(r.values()))

        self.varsimID.set(cfg['simulation']['ID'])
        self.varsimdatatype.set(cfg['simulation']['datatype'])
        self.varsimulation.set(cfg['simulation']['simulation'])
        self.varbenchmark.set(cfg['simulation']['benchmark'])
        self.varfastmode.set(cfg['simulation']['fastmode'])

    def checkoptions(self):

        var = self.varsimdatatype.get()

        if var == 1: #  Weather File
            self.varbenchmark.set(False)
            self.cbbenchmark['state'] ='disabled'
            self.cbloadsitedata['state'] = 'normal'
        elif var == 2: #  Field Data File
            self.varloadsitedata.set(False)
            self.cbbenchmark['state']='normal'
            self.cbloadsitedata['state'] = 'disabled'
        else:
            pass


    def buildSimulationFrame(self):

        self.varsimID = tk.StringVar(self.fr_simulation)
        self.varsimdatatype = tk.IntVar(self.fr_simulation)
        self.varsimulation = tk.BooleanVar(self.fr_simulation)
        self.varbenchmark = tk.BooleanVar(self.fr_simulation)
        self.varfastmode = tk.BooleanVar(self.fr_simulation)
        self.varloadsitedata = tk.BooleanVar(self.fr_simulation)
        self.varloadsitedata.set(False)

        self.varsitename = tk.StringVar(self.fr_simulation)
        self.varsitelat = tk.DoubleVar(self.fr_simulation)
        self.varsitelong = tk.DoubleVar(self.fr_simulation)
        self.varsitealt = tk.DoubleVar(self.fr_simulation)

        self.vardatafileurl = tk.StringVar(self.fr_simulation)
        self.vardatafilename = tk.StringVar(self.fr_simulation)
        self.vardatafilepath = tk.StringVar(self.fr_simulation)

        self.lbsimID = ttk.Label(
            self.fr_simulation,
            text= "Simulation ID").grid(
                row = 0, column = 0, sticky='W', padx=2, pady=5)
        self.ensimID = ttk.Entry(
            self.fr_simulation,
            textvariable = self.varsimID ).grid(
                row = 0, column = 1, sticky='W', padx=2, pady=5)

        self.lbdatatype = ttk.Label(
            self.fr_simulation,
            text= "Choose a Data Source Type:").grid(
                row = 1, column = 0, columnspan = 2,  sticky='W', padx=2, pady=5)

        self.rbweather = tk.Radiobutton(
            self.fr_simulation,
            padx = 5,
            text = "Weather File",
            variable=self.varsimdatatype, value=1,
            command = lambda: self.checkoptions()).grid(
                row = 2, column = 0, sticky='W', padx=2, pady=5)

        #  RadioButton for Field Data File
        self.rbfielddata = tk.Radiobutton(
            self.fr_simulation,
            padx = 5,
            text = "Field Data File",
            variable=self.varsimdatatype, value=2,
            command = lambda: self.checkoptions()).grid(
                row = 2, column = 1, sticky='W', padx=2, pady=5)

        self.btselectdatasource = ttk.Button(
            self.fr_simulation, text='Select File',
            command= lambda : self.dataLoadDialog(
                "Select File",
                labeltext = "Select File"))
        self.btselectdatasource.grid(
                row = 3, column = 0, sticky='W', padx=2, pady=5)

                #  Data source path
        self.vardatafileurl.set('Data source file path...')
        self.lbdatasourcepath = ttk.Label(
            self.fr_simulation, textvariable = self.vardatafileurl).grid(
                row = 3, column= 1, columnspan=4, sticky='W', padx=2, pady=5)



        #  Checkbox for Simulation
        self.lbsimulation = ttk.Label(
            self.fr_simulation, text= "Run test type...").grid(
                row = 4, column = 0, sticky='W', padx=2, pady=5)
        self.cbsimulation = ttk.Checkbutton(
            self.fr_simulation,
            text="Simulation",
            variable=self.varsimulation)
        self.cbsimulation.grid(
            row = 4, column = 1, sticky='W', padx=2, pady=5)



        #  Checkbox for Benchmark
        # self.lbbenckmark = ttk.Label(
        #     self.fr_simulation,
        #     text= "Benchmark").grid(
        #         row = 4, column = 2, sticky='W', padx=5, pady=5)
        self.cbbenchmark = ttk.Checkbutton(
            self.fr_simulation,
            text="Benchmark",
            variable=self.varbenchmark)
        self.cbbenchmark.grid(
            row = 4, column = 2, sticky='W', padx=2, pady=5)

        #  Checkbox for fastmode
        # self.lbfastmode = ttk.Label(
        #     self.fr_simulation,
        #     text= "Fast mode (Only prototype loop)").grid(
        #         row = 3, column = 0, sticky='W', padx=5, pady=5)
        self.cbfastmode = ttk.Checkbutton(
            self.fr_simulation,
            text="Fast mode ON",
            variable=self.varfastmode).grid(
                row = 4, column = 3, sticky='W', padx=80, pady=5)

        #  Checkbox to load site data from weather file
        # self.lbloadsitedata = ttk.Label(
        #     self.fr_simulation,
        #     text= "Copy site data from weather file").grid(
        #         row = 4, column = 4, sticky='W', padx=2, pady=5)

        self.cbloadsitedata = ttk.Checkbutton(
            self.fr_simulation,
            text="Load site data from weather file",
            variable=self.varloadsitedata)
        self.cbloadsitedata.grid(
            row = 4, column = 4, sticky='W', padx=2, pady=5)



        self.lbsitename = ttk.Label(
              self.fr_simulation, text= "Name").grid(row = 5, column = 0, sticky='W', padx=2, pady=5)
        self.ensitename = ttk.Entry(
              self.fr_simulation, textvariable= self.varsitename).grid(row = 5, column = 1, sticky='W', padx=2, pady=5)
        self.lbsitelat = ttk.Label(
              self.fr_simulation, text= "Latitude").grid(row = 5, column = 2, sticky='W', padx=2, pady=5)
        self.ensitelat = ttk.Entry(
              self.fr_simulation, textvariable = self.varsitelat).grid(row = 5, column = 3, sticky='W', padx=2, pady=5)
        self.lbsitelong = ttk.Label(
              self.fr_simulation, text= "Longitude").grid(row = 5, column = 4, sticky='W', padx=2, pady=5)
        self.ensitelong = ttk.Entry(
              self.fr_simulation, textvariable= self.varsitelong).grid(row = 5, column = 5, sticky='W', padx=2, pady=5)
        self.lbsitealt = ttk.Label(
              self.fr_simulation, text= "Altitude").grid(row = 5, column = 6, sticky='W', padx=2, pady=5)
        self.ensitealt = ttk.Entry(
              self.fr_simulation, textvariable = self.varsitealt).grid(row = 5, column = 7, sticky='W', padx=2, pady=5)

        self.varsimdatatype.set(1)  #  1 for Weather File, 2 for Field Data File
        self.checkoptions()


        #  Button Load Configuration File
        # self.btloadcfgsimulation = ttk.Button(
        #     self.fr_simulation,
        #     text= "Load config",
        #     command= lambda : self.simulationLoadDialog(
        #         "Load config file",
        #         labeltext = "Load Simulation Config"))
        # self.btloadcfgsimulation.grid(row = 4, column = 0, sticky='W', padx=5, pady=5)

        # self.btsavecfgsimulation = ttk.Button(
        #     self.fr_simulation,
        #     text= "Save config",
        #     command= lambda : self.simulation_save_dialog(
        #         "Save config file",
        #         labeltext = "Save Simulation Config"))
        # self.btsavecfgsimulation.grid(row = 4, column = 1, sticky='W', padx=5, pady=5)


    #  Solar Field Layout contruction Tab

    def solarfieldLoadDialog(self, title, labeltext = '' ):

        path = askopenfilename(initialdir = self._DIR['saved_configurations'],
                               title = "choose your file",
                               filetypes = [("JSON files", "*.json")])

        with open(path) as cfg_file:
            cfg = json.load(cfg_file,
                                     parse_float= float,
                                     parse_int= int)

        datarow = []
        for r in cfg['solarfield']['subfields']:
            datarow.append(list(r.values()))

        self.solarfield_table = table.Tk_Table(
                        self.fr_solarfield,
                        ["NAME", "LOOPS"],
                        row_numbers=True,
                        stripped_rows=("white", "#f2f2f2"),
                        select_mode="none",
                        cell_anchor="center",
                        adjust_heading_to_content=True)



        self.solarfield_table.table_data = datarow
        self.solarfield_table.grid(row = 8, column = 0, columnspan =4)

        self.varsolarfieldname.set(cfg['solarfield']['name'])
        self.vartin.set(cfg['solarfield']['rated_tin'])
        self.vartout.set(cfg['solarfield']['rated_tout'])
        self.varpin.set(cfg['solarfield']['rated_pin'])
        self.varpout.set(cfg['solarfield']['rated_pout'])
        self.vartmin.set(cfg['solarfield']['tmin'])
        self.vartmax.set(cfg['solarfield']['tmax'])
        self.varratedmassflow.set(cfg['solarfield']['rated_massflow'])
        self.varrecirculation.set(cfg['solarfield']['recirculation_massflow'])
        self.varscas.set(cfg['solarfield']['loop']['scas'])
        self.varhces.set(cfg['solarfield']['loop']['hces'])

        #self.enratedtin.insert(0, cfg['solarfield']['rated_tin'])
        # self.enratedtout.delete(0, tk.END)
        # self.enratedtout.insert(0, cfg['solarfield']['rated_tout'])

        self.fr_solarfield.update()

    def solarfield_save_dialog(self, title, labeltext = '' ):

        #encoder.FLOAT_REPR = lambda o: format(o, '.2f')
        f = asksaveasfile(initialdir = self._DIR['saved_configurations'],
                               title = "choose your file name",
                               filetypes = [("JSON files", "*.json")],
                               defaultextension = "json")

        cfg = dict({"solarfield" : {}})
        cfg["solarfield"].update(dict({"name" : self.enname.get()}))
        cfg["solarfield"].update(dict({"rated_tin" : self.enratedtin.get()}))
        cfg["solarfield"].update(dict({"rated_tout" : self.enratedtout.get()}))
        cfg["solarfield"].update(dict({"rated_pin" : self.enratedpin.get()}))
        cfg["solarfield"].update(dict({"rated_pout" : self.enratedpout.get()}))
        cfg["solarfield"].update(dict({"rated_masslflow" : self.enratedmassflow.get()}))
        cfg["solarfield"].update(dict({"recirculation_massflow" : self.enrecirculationmassflow.get()}))
        cfg["solarfield"].update(dict({"tmin" : self.entmin.get()}))
        cfg["solarfield"].update(dict({"tmin" : self.entmax.get()}))
        cfg["solarfield"].update(dict({"loop": dict({"scas": self.enscas.get(),
                                                     "hces": self.enhces.get()})}))


        datarow = list(self.solarfield_table.table_data)
        dictkeys =["name", "loops"]

        subfields = []

        for r in datarow:
            sf = {}
            index = 0
            for v in r:
                k = dictkeys[index]
                sf[k]= self.to_number(v)
                index += 1
            subfields.append(sf)

        cfg['solarfield'].update({"subfields" : subfields})
        # cfg_settings['solarfield'].update(dict(cfg))
        f.write(json.dumps(cfg))
        f.close()


    def buildSolarFieldFrame(self):

        self.varsolarfieldname = tk.StringVar()
        self.lbname = ttk.Label(self.fr_solarfield, text ="Solar Field Name" )
        self.lbname.grid(row = 0, column = 0, sticky='W', padx=5, pady=5)
        self.enname = ttk.Entry(self.fr_solarfield, textvariable = self.varsolarfieldname)
        self.enname.grid(row = 0, column = 1, sticky='W', padx=5, pady=5)

        self.vartin = tk.StringVar()
        self.lbratedtin = ttk.Label(self.fr_solarfield, text= "Rated Tin [K]")
        self.lbratedtin.grid(row = 1, column = 0, sticky='W', padx=5, pady=5)
        self.enratedtin = ttk.Entry(self.fr_solarfield, textvariable = self.vartin)
        self.enratedtin.grid(row = 1, column = 1, sticky='W', padx=5, pady=5)

        self.vartout = tk.StringVar()
        self.lbratedtout = ttk.Label(self.fr_solarfield, text= "Rated Tout [K]")
        self.lbratedtout.grid(row = 1, column = 2, sticky='W', padx=5, pady=5)
        self.enratedtout = ttk.Entry(self.fr_solarfield, textvariable= self.vartout)
        self.enratedtout.grid(row = 1, column = 3, sticky='W', padx=5, pady=5)

        self.varpin = tk.StringVar()
        self.lbratedpin = ttk.Label(self.fr_solarfield, text= "Rated Pin [Pa]")
        self.lbratedpin.grid(row = 2, column = 0, sticky='W', padx=5, pady=5)
        self.enratedpin = ttk.Entry(self.fr_solarfield, textvariable= self.varpin)
        self.enratedpin.grid(row = 2, column = 1, sticky='W', padx=5, pady=5)

        self.varpout = tk.StringVar()
        self.lbratedpout = ttk.Label(self.fr_solarfield, text= "Rated Pout [Pa]")
        self.lbratedpout.grid(row = 2, column = 2, sticky='W', padx=5, pady=5)
        self.enratedpout = ttk.Entry(self.fr_solarfield, textvariable = self.varpout)
        self.enratedpout.grid(row = 2, column = 3, sticky='W', padx=5, pady=5)

        self.vartmin = tk.StringVar()
        self.lbtmin = ttk.Label(self.fr_solarfield, text= "Tmin [K]")
        self.lbtmin.grid(row = 3, column = 0, sticky='W', padx=5, pady=5)
        self.entmin = ttk.Entry(self.fr_solarfield, textvariable = self.vartmin)
        self.entmin.grid(row = 3, column = 1, sticky='W', padx=5, pady=5)

        self.vartmax = tk.StringVar()
        self.lbtmax = ttk.Label(self.fr_solarfield, text= "Tmax [K]")
        self.lbtmax.grid(row = 3, column = 2, sticky='W', padx=5, pady=5)
        self.entmax = ttk.Entry(self.fr_solarfield, textvariable = self.vartmax)
        self.entmax.grid(row = 3, column = 3, sticky='W', padx=5, pady=5)

        self.varratedmassflow = tk.StringVar()
        self.lbratedmassflow = ttk.Label(self.fr_solarfield, text= "Rated massflow (per loop) [Kg/s]")
        self.lbratedmassflow.grid(row = 4, column = 0, sticky='W', padx=5, pady=5)
        self.enratedmassflow = ttk.Entry(self.fr_solarfield, textvariable = self.varratedmassflow)
        self.enratedmassflow.grid(row = 4, column = 1, sticky='W', padx=5, pady=5)

        self.varrecirculation = tk.StringVar()
        self.lbrecirculationmassflow = ttk.Label(self.fr_solarfield, text= "Rated massflow (per loop) [Kg/s]")
        self.lbrecirculationmassflow.grid(row = 5, column = 0, sticky='W', padx=5, pady=5)
        self.enrecirculationmassflow = ttk.Entry(self.fr_solarfield, textvariable = self.varrecirculation)
        self.enrecirculationmassflow.grid(row = 5, column = 1, sticky='W', padx=5, pady=5)

        self.varscas = tk.StringVar()
        self.lbscas = ttk.Label(self.fr_solarfield, text= "SCAs per Loop")
        self.lbscas.grid(row = 6, column = 0, sticky='W', padx=5, pady=5)
        self.enscas = ttk.Entry(self.fr_solarfield, textvariable = self.varscas)
        self.enscas.grid(row = 6, column = 1, sticky='W', padx=5, pady=5)

        self.varhces = tk.StringVar()
        self.lbhces = ttk.Label(self.fr_solarfield, text= "HCEs per SCA")
        self.lbhces.grid(row = 6, column = 2, sticky='W', padx=5, pady=5)
        self.enhces = ttk.Entry(self.fr_solarfield, textvariable = self.varhces)
        self.enhces.grid(row = 6, column = 3, sticky='W', padx=5, pady=5)

        self.solarfield_table = table.Tk_Table(
            self.fr_solarfield,
            ["SUBFIELD NAME", "NUMBER OF LOOPS"],
            row_numbers=True,
            stripped_rows=("white", "#f2f2f2"),
            select_mode="none",
            cell_anchor="center",
            adjust_heading_to_content=True)
        self.solarfield_table.grid(row = 7, column = 0, columnspan =2)

        # self.btloadcfgsolarfield = ttk.Button(
        #     self.fr_solarfield,
        #     text= "Load config",
        #     command= lambda : self.solarfieldLoadDialog(
        #         "Load config file",
        #         labeltext = "Solar Field Config"))
        # self.btloadcfgsolarfield.grid(row = 7, column = 0, sticky='W', padx=5, pady=5)

        # self.btsavecfgsolarfield = ttk.Button(
        #     self.fr_solarfield,
        #     text= "Save config",
        #     command= lambda : self.solarfield_save_dialog(
        #         "Save config file",
        #         labeltext = "Solar Field Config"))
        # self.btsavecfgsolarfield.grid(row = 7, column = 1, sticky='W', padx=5, pady=5)

        self.btnewrow = ttk.Button(
            self.fr_solarfield,
            text= "Insert",
            command = lambda : self.__insert_rows__(self.solarfield_table))
        self.btnewrow.grid(row = 9, column = 0, sticky='W', padx=5, pady=5)
        self.btdelrows = ttk.Button(
            self.fr_solarfield,
            text="Delete",
            command = lambda : self.__del_rows__(self.solarfield_table))
        self.btdelrows.grid(row = 9, column = 1, sticky='W', padx=5, pady=5)


    #  Site & Weather contruction Tab

    def dataLoadDialog(self, title, labeltext = ''):

        if self.varsimdatatype.get() == 1:
            path = askopenfilename(initialdir = self._DIR['weather_files'],
                                   title = "choose your file",
                                   filetypes = (("TMY files","*.tm2"),
                                        ("TMY files","*.tm3"),
                                        ("csv files","*.csv"),
                                        ("all files","*.*")))
            self.vardatafilepath.set(os.path.dirname(path)+"/")
            self.vardatafilename.set(os.path.basename(path))
            self.vardatafileurl.set(os.path.dirname(path)+"/" +
                                    os.path.basename(path))
        elif self.varsimdatatype.get() == 2:
            path = askopenfilename(initialdir = self._DIR['fielddata_files'],
                       title = "choose your file",
                       filetypes = (("csv files","*.csv"),
                                    ("all files","*.*")))
            self.vardatafilepath.set(os.path.dirname(path)+"/")
            self.vardatafilename.set(os.path.basename(path))
            self.vardatafileurl.set(os.path.dirname(path)+"/" +
                                    os.path.basename(path))
        else:

            tk.messagebox.showwarning(
                title='Warning',
                message='Check Data Source Selection')


        if self.varloadsitedata.get():

            strfilename, strext = os.path.splitext(path)

            if  strext == ".csv":
                weatherdata = pvlib.iotools.tmy.read_tmy3(path)
                file = path
            elif (strext == ".tm2" or strext == ".tmy"):
                weatherdata = pvlib.iotools.tmy.read_tmy2(path)
                file = path
            elif strext == ".xls":
                pass
            else:
                print("unknow extension ", strext)
                return

            self.varsitename.set(weatherdata[1]['City'])
            self.varsitelat.set(weatherdata[1]['latitude'])
            self.varsitelong.set(weatherdata[1]['longitude'])
            self.varsitealt.set(weatherdata[1]['altitude'])


    def dataSaveDialog(self, title, labeltext = '' ):

        #encoder.FLOAT_REPR = lambda o: format(o, '.2f')
        f = asksaveasfile(initialdir = self._DIR['site_files'],
                               title = "choose your file name",
                               filetypes = [("JSON files", "*.json")],
                               defaultextension = "json")

        f.write(json.dumps(cfg_settings['weather']))
        f.close()

    def buildDataFrame(self):

        # self.varsitename = tk.StringVar(self.fr_data)
        # self.varsitelat = tk.DoubleVar(self.fr_data)
        # self.varsitelong = tk.DoubleVar(self.fr_data)
        # self.varsitealt = tk.DoubleVar(self.fr_data)

        # self.lbsitename = ttk.Label(
        #     self.fr_data, text= "Name").grid(row = 0, column = 0, sticky='W', padx=5, pady=5)
        # self.ensitename = ttk.Entry(
        #     self.fr_data, textvariable= self.varsitename).grid(row = 0, column = 1, sticky='W', padx=5, pady=5)
        # self.lbsitelat = ttk.Label(
        #     self.fr_data, text= "Latitude").grid(row = 0, column = 2, sticky='W', padx=5, pady=5)
        # self.ensitelat = ttk.Entry(
        #     self.fr_data, textvariable = self.varsitelat).grid(row = 0, column = 3, sticky='W', padx=5, pady=5)
        # self.lbsitelong = ttk.Label(
        #     self.fr_data, text= "Longitude").grid(row = 0, column = 4, sticky='W', padx=5, pady=5)
        # self.ensitelong = ttk.Entry(
        #     self.fr_data, textvariable= self.varsitelong).grid(row = 0, column = 5, sticky='W', padx=5, pady=5)
        # self.lbsitealt = ttk.Label(
        #     self.fr_data, text= "Altitude").grid(row = 0, column = 6, sticky='W', padx=5, pady=5)
        # self.ensitealt = ttk.Entry(
        #     self.fr_data, textvariable = self.varsitealt).grid(row = 0, column = 7, sticky='W', padx=5, pady=5)

        self.lbweatherfile = ttk.Label(self.fr_data, text= "Weather File")
        self.lbweatherfile.grid(row = 1, column = 0)
        self.enweatherfile = ttk.Entry(self.fr_data)
        self.enweatherfile.insert(0, "Path to Weather File")
        self.enweatherfile.grid(row = 1, column = 1, columnspan = 4)
        self.btloadcfgweather = ttk.Button(
            self.fr_data,
            text= "Load Weather",
            command= lambda : self.dataLoadDialog(
                "Load config file",
                labeltext="Weather File"))
        self.btloadcfgweather.grid(row = 1, column = 6)
        self.btsavecfgweather = ttk.Button(
                text= "Save config",
                command= lambda : self.dataSaveDialog(
                        self.fr_data,
                        "nombre",
                        "Save config file",
                        labeltext="Weather"))
        #self.btsavecfgweather.grid(row = 1, column = 7)



    #  Fluid contruction Tab

    def select_fluid_from_csv(self):

        # abre una ventana de seleccón de csv.

        self.loadWindow= tk.Tk()
        self.loadWindow.attributes('-fullscreen', False)
        self.loadWindow.title("Selection Window")
        # w, h = self.root.winfo_screenwidth(), self.root.winfo_screenheight()
        # self.root.geometry("%dx%d+0+0" % (w, h))

        # Menu
        self.loadWindow.menubar = tk.Menu(self.loadWindow)
        self.loadWindow.load_menu = tk.Menu(self.loadWindow.menubar,tearoff=0)
        self.loadWindow.load_menu.add_command(label='Open', command=self.open_csv())
        self.loadWindow.load_menu.add_separator()
        self.loadWindow.load_menu.add_command(label='Exit', command=None)
        self.loadWindow.menubar.add_cascade(label='Select File', menu=self.loadWindow.load_menu)

        self.loadWindow.config(menu=self.loadWindow.menubar)

        # carga el csv en una nueva ventana, en formato tabla
        # el usuario selección un registro (un fluido)
        # los datos del registo se pasan en json al programa principal


        pass


    def open_csv(self, path=None):

        df = pd.DataFrame()


        try:
            if path is None:
                root = Tk()
                root.withdraw()
                path = askopenfilename(initialdir = ".fielddata_files/",
                                   title = "choose your file",
                                   filetypes = (("csv files","*.csv"),
                                                ("all files","*.*")))
                root.update()
                root.destroy()

                if path is None:
                    return
                else:
                    strfilename, strext = os.path.splitext(path)

                    if  strext == ".csv":
                        print("csv........")

                        df = pd.read_csv(path, sep=';',
                                                     decimal= ',')
                        file = path
                    else:
                        print("unknow extension ", strext)
                        return
            else:
                strfilename, strext = os.path.splitext(path)

                if  strext == ".csv":
                    print("csv...")
                    df = pd.read_csv(path, sep=';',
                                                 decimal= ',')
                    file = path
                else:
                    print("unknow extension ", strext)
                    return

        except Exception:
            raise
            # txMessageBox.showerror('Error loading  File',
            #                        'Unable to open file: %r', self.file)


    def fluid_load_dialog(self):

        path = askopenfilename(initialdir = self._DIR['fluid_files'],
                               title = "choose your file",
                               filetypes = [("JSON files", "*.json")])

        with open(path) as cfg_file:
            cfg = json.load(cfg_file, parse_float= float, parse_int= int)

        cp_coefs = [[]]*7
        rho_coefs = [[]]*7
        mu_coefs  = [[]]*7
        kt_coefs  = [[]]*7
        h_coefs = [[]]*7
        t_coefs = [[]]*7

        temp_cp = ["cp"]
        temp_rho = ["rho"]
        temp_mu = ["mu"]
        temp_kt = ["kt"]
        temp_h = ["h"]
        temp_t = ["t"]

        grades = ["Factor"]
        grades.extend([0, 1, 2, 3, 4, 5])

        temp_cp.extend(list(cfg['hot_fluid']['cp']))
        temp_rho.extend(list(cfg['hot_fluid']['rho']))
        temp_mu.extend(list(cfg['hot_fluid']['mu']))
        temp_kt.extend(list(cfg['hot_fluid']['kt']))
        temp_h.extend(list(cfg['hot_fluid']['h']))
        temp_t.extend(list(cfg['hot_fluid']['t']))

        for index in range(len(temp_cp)):
            cp_coefs[index] = temp_cp[index]
        cp_coefs[1:] = [float(s) for s in cp_coefs[1:]]
        for index in range(len(temp_rho)):
            rho_coefs[index] = temp_rho[index]
        rho_coefs[1:] = [float(s) for s in rho_coefs[1:]]
        for index in range(len(temp_mu)):
            mu_coefs[index] = temp_mu[index]
        mu_coefs[1:] = [float(s) for s in mu_coefs[1:]]
        for index in range(len(temp_kt)):
            kt_coefs[index] = temp_kt[index]
        kt_coefs[1:] = [float(s) for s in kt_coefs[1:]]
        for index in range(len(temp_h)):
            h_coefs[index] = temp_h[index]
        h_coefs[1:] = [float(s) for s in h_coefs[1:]]
        for index in range(len(temp_t)):
            t_coefs[index] = temp_t[index]
        t_coefs[1:] = [float(s) for s in t_coefs[1:]]

        datarow = []
        datarow.append(cp_coefs)
        datarow.append(rho_coefs)
        datarow.append(mu_coefs)
        datarow.append(kt_coefs)
        datarow.append(h_coefs)
        datarow.append(t_coefs)

        self.fluid_table.table_data = datarow

        self.varminfluidtemp.set(cfg['hot_fluid']['tmin'])
        self.varmaxfluidtemp.set(cfg['hot_fluid']['tmax'])
        self.varnamefluid.set(cfg['hot_fluid']['name'])

        self.fr_fluid.update()

    def fluid_save_dialog(self, title, labeltext = '' ):

        #encoder.FLOAT_REPR = lambda o: format(o, '.2f')
        f = asksaveasfile(
            initialdir = self._DIR['fluid_files'],
            title = "choose your file name",
            filetypes = [("JSON files", "*.json")],
            defaultextension = "json")

        cfg = dict({"hot_fluid" : {}})
        cfg["hot_fluid"].update(dict({"name" : self.ennamefluid.get()}))

        datarow = list(self.fluid_table.table_data)

        for r in datarow:
            param_name = r[0]
            param_values = list(map(self.to_number, r[1:]))
            cfg["hot_fluid"].update(dict({param_name : param_values}))

        cfg['hot_fluid'].update({"tmax" : float(self.entmax.get())})
        cfg['hot_fluid'].update({"tmin" : float(self.entmin.get())})
        f.write(json.dumps(cfg))
        f.close()

    def checkfluid(self):

        var = self.varfluidtable.get()

        if var == 1: #  Fluid from table
            self.varcoolpropID.set('')
            #self.encoolpropID['state'] ='disabled'
            self.cmbcoolpropID['state'] = 'disabled'
            self.ennamefluid['state'] = 'normal'
            self.btloadcfgfluid['state'] = 'normal'
            self.entminfluid['state'] = 'normal'
            self.entmaxfluid['state'] = 'normal'
        elif var == 2: #  Field Data File
            self.varmaxfluidtemp.set('')
            self.varminfluidtemp.set('')
            self.varnamefluid.set('')
            # self.encoolpropID['state'] ='normal'
            self.cmbcoolpropID['state'] = 'readonly'
            self.ennamefluid['state'] = 'disabled'
            self.btloadcfgfluid['state'] = 'disabled'
            self.entminfluid['state'] = 'disabled'
            self.entmaxfluid['state'] = 'disabled'
            self.fluid_table.table_data = []
        else:
            pass

    def buildFluidFrame(self):

        self.varfluidtable = tk.IntVar(self.fr_fluid)
        self.varfluidtable.set(1)  # 1 for Table, 2 for CoolProp Library
        self.varcoolpropID = tk.StringVar(self.fr_fluid)

        #  RadioButton for fluid from table
        self.rbfluidtable = tk.Radiobutton(
            self.fr_fluid,
            padx = 2,
            text = "Fluid Data from table",
            variable=self.varfluidtable, value=1,
            command = lambda: self.checkfluid()).grid(
                row = 0, column = 0, sticky='W', padx=2, pady=5)

        #  RadioButton for fluid from library
        self.rbfluidlib = tk.Radiobutton(
            self.fr_fluid,
            padx = 2,
            text = "Fluid Data from CoolProp",
            variable=self.varfluidtable, value=2,
            command = lambda: self.checkfluid()).grid(
                row = 0, column = 3, sticky='W', padx=2, pady=5)

        self.lbcoolpropID = tk.Label(self.fr_fluid, text= "CoolProp ID (INCOMP::xxxx)")
        self.lbcoolpropID.grid(row = 0, column = 4, sticky='W', padx=2, pady=5)

        self.cmbcoolpropID = ttk.Combobox(self.fr_fluid)
        self.cmbcoolpropID["values"] = self._COOLPROP_FLUIDS
        self.cmbcoolpropID['state']='readonly'
        self.cmbcoolpropID.grid(row = 0, column = 5, sticky='W', padx=2, pady=5)

        # self.encoolpropID = ttk.Entry(self.fr_fluid, textvariable=self.varcoolpropID)
        # self.encoolpropID.grid(row = 0, column = 3)

        self.varnamefluid = tk.StringVar(self.fr_fluid)
        self.lbnamefluid = ttk.Label(self.fr_fluid, text= "Name")
        self.lbnamefluid.grid(row = 1, column = 0, sticky='W', padx=2, pady=5)
        self.ennamefluid = ttk.Entry(self.fr_fluid, textvariable=self.varnamefluid)
        self.ennamefluid.grid(row = 1, column = 1, sticky='W', padx=2, pady=5)

        self.btloadcfgfluid = ttk.Button(
            self.fr_fluid, text= "Load config",
            command= lambda : self.fluid_load_dialog())
        self.btloadcfgfluid.grid(row = 1, column = 2, sticky='W', padx=2, pady=5)

        self.varminfluidtemp = tk.DoubleVar(self.fr_fluid)
        self.lbtminfluid = ttk.Label(self.fr_fluid, text = "Tmin[K]")
        self.lbtminfluid.grid(row=2,column = 0, sticky='W', padx=2, pady=5)
        self.entminfluid = ttk.Entry(self.fr_fluid, textvariable=self.varminfluidtemp)
        self.entminfluid.grid(row = 2, column = 1, sticky='W', padx=2, pady=5)

        self.varmaxfluidtemp = tk.DoubleVar(self.fr_fluid)
        self.lbtmaxfluid = ttk.Label(self.fr_fluid, text = "Tmax [K]")
        self.lbtmaxfluid.grid(row=3,column = 0, sticky='W', padx=2, pady=5)
        self.entmaxfluid = ttk.Entry(self.fr_fluid, textvariable= self.varmaxfluidtemp)
        self.entmaxfluid.grid(row = 3, column = 1, sticky='W', padx=2, pady=5)

        self.fluid_table = table.Tk_Table(
                self.fr_fluid,
                ["Parameter [x]",
                 "       A x^0       ", "       B x^1       ",
                 "       C x^2       ", "       D x^3       ",
                 "       E x^4       ", "       F x^5       "],
                row_numbers=True,
                stripped_rows = ("white","#f2f2f2"),
                select_mode = "none",
                cell_anchor="center",
                adjust_heading_to_content = True)
        self.fluid_table.grid(row = 4, column = 0, columnspan = 6, sticky='W', padx=2, pady=5)

        self.btsavecfgfluid = ttk.Button(
                self.fr_fluid,
                text= "Save config",
                command= lambda : self.fluid_save_dialog(
                    "Save config file",
                    labeltext="HTF Config"))
        self. btsavecfgfluid.grid(row = 5, column = 0, sticky='W', padx=2, pady=5)

        self.checkfluid()

"""

    # HCE Construction tab

    def build_hce_table(self):

        self.hce_table = table.Tk_Table(
            self.fr_hce,
            param_name,
            row_numbers=True,
            stripped_rows = ("white","#fr_datafr_datafr_data"),
            select_mode = "none",
            cell_anchor="center",
            adjust_heading_to_content = True)

        self.hce_table.grid(row = 0, columnspan = 2)



    def hce_load_dialog(fr_hce, hce_table, title, labeltext = '' ):

        path = askopenfilename(initialdir = _DIR['hce_files'],
                                   title = "choose your file",
                                   filetypes = [("JSON files", "*.json")])
        with open(path) as cfg_file:
            cfg = json.load(cfg_file, parse_float= float, parse_int= int)

        data = list(cfg)
        param_names = data[0].values()

        datarow = []
        for r in data:
            datarow.append(list(r.values()))

        self.hce_table.table_data = datarow

    # param_name = ["Name", "Description", "Condition", "Broken", "Bellows",
    #              "Transmissivity", "Absorption", "Unaccounted",
    #              "A0", "A1", "A2", "A3", "A4", "A5", "A6", "Factor"]
    def hce_save_dialog(fr_hce, hce_table, title, labeltext = '' ):

        #encoder.FLOAT_REPR = lambda o: format(o, '.2f')
        file = asksaveasfile(initialdir = _DIR['hce_files'],
                               title = "choose your file name",
                               filetypes = [("JSON files", "*.json")],
                               defaultextension = "json")
        cfg = []
        data = hce_table.table_data
        # param_name = ["Name", "Description", "Condition", "Broken", "Bellows",
        #              "Transmissivity", "Absorption", "Unaccounted",
        #              "A0", "A1", "A2", "A3", "A4", "A5", "A6", "Factor"]
        for r in data:
            param_values = list(map(to_number, r[0:]))
            cfg.append(dict(zip(param_name, param_values)))

        file.write(json.dumps(cfg))
        file.close()

    # hce_table = table.Tk_Table(
    #                 fr_hce,
    #                 ["Name", "Description", "Condition", "Broken", "Bellows",
    #                  "Transmissivity", "Absorption", "Unaccounted",
    #                  "A0", "A1", "A2", "A3", "A4", "A5", "A6", "Factor"],
    #                 row_numbers=True,
    #                 stripped_rows = ("white","#fr_datafr_datafr_data"),
    #                 select_mode = "none",
    #                 cell_anchor="center",
    #                 adjust_heading_to_content = True)

    def buildHCEFrame(self):



        self.btloadcfghce = ttk.Button(
            self.fr_hce, text= "Load config",
            command= lambda : self.hce_load_dialog(
                self.fr_hce,
                self.hce_table,
                "Load config file",
                labeltext="HCE Config"))
        self.btloadcfghce.grid(row = 1, column = 0)
        self.btsavecfghce = ttk.Button(
            self.fr_hce,
            text= "Save config",
            command= lambda : hce_save_dialog(
                self.fr_hce,
                self.hce_table,
                "Save config file",
                labeltext="HCE Config"))
        self.btsavecfghce.grid(row = 1, column = 1)
        self.fr_hce.update()



    # SCA Construction tab

    def sca_load_dialog(fr_sca, sca_table, title, labeltext = '' ):

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
        fr_sca.update()

    def sca_save_dialog(fr_sca, sca_table, title, labeltext = '' ):

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
                    fr_sca, ["Name","SCA Length","Aperture","Aperture Area","Focal Len",
                        "IAM Coefficient F0","IAM Coefficient fr_solarfield", "AM Coefficient fr_data",
                        "Track  Twist", "Geom.Accuracy", "Reflectance", "Cleanliness",
                        "Dust", "Factor", "Availability"],
                    row_numbers=True,
                    stripped_rows = ("white","#fr_datafr_datafr_data"),
                    select_mode = "none",
                    cell_anchor="center",
                    adjust_heading_to_content = True)

    #lbnamesca = ttk.Label(fr_sca, text= "Name")
    #lbnamesca.grid(row = 0, column = 0)
    #ennamesca = ttk.Entry(fr_sca, text= "Name")
    #ennamesca.grid(row = 0, column = 1)
    #
    #entmaxsca = ttk.Entry(fr_sca, text= "Tmax")
    #entminsca = ttk.Entry(fr_sca, text= "Tmin")
    #entmaxsca.grid(row = 0, column = 3)
    #entminsca.grid(row = 0, column = 4)

    btloadcfgsca = ttk.Button(
        fr_sca, text= "Load config",
        command= lambda : sca_load_dialog(
                fr_sca, sca_table,
                "Load config file", labeltext="SCA Config"))
    btloadcfgsca.grid(row = 0, column = 2)

    btsavecfgsca = ttk.Button(
            fr_sca,
            text= "Save config",
            command= lambda : sca_save_dialog(
                    fr_sca, sca_table,
                    "Save config file", labeltext="HTF Config"))
    btsavecfgsca.grid(row = 0, column = 5)


    def buildSCAFrame(self):
"""






if __name__ == '__main__':

    try:
        from Tkinter import Tk
        import tkMessageBox as messagebox
    except ImportError:
        from tkinter import Tk
        from tkinter import messagebox

    interface = Interface()

    interface.root.mainloop()