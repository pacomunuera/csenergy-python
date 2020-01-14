# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 11:53:18 2020

@author: fmunuera
"""
from tkinter import * 
from tkinter.filedialog import askopenfilename
import pandas as pd
from datetime import datetime
import numpy as np
import pvlib as pvlib






#print("Definir el  número de subcampos", end='\n')
#input_solarfields = int(input())
#print("Definir el número de lazos en cada subcampo", end='\n')
#input_loops = int(input())
#print("Definir el número de SCA en cada lazo", end='\n')
#input_sca = int(input())
#print("Definir el número de HCE en cada SCA", end='\n')
#input_HCE = int(input())

# Window and form
class main(Frame):
    def __init__(self, master):
        Frame.__init__(self, master = None)
        self.master.title("CSENERGY")
        self.master.geometry("800x600")
        
        bt_site = Button(self.master, text="Site", width=20, 
                         command=self.site_dialog)
        bt_site.grid(row=0, column=0)
        lb_site = Label(self.master, text="Site selection")
        lb_site.grid(row=1, column=0)
        
        bt_plant = Button(self.master, text="Plant", width=20, 
                          command=self.plant_dialog)
        bt_plant.grid(row=0, column=1)
        lb_plant = Label(self.master, text="Plant configuration")
        lb_plant.grid(row=1, column=1)
        
        bt_weather = Button(self.master, text="Weather", width=20, 
                            command=self.weather_dialog)
        bt_weather.grid(row=0, column=2)
        lb_weather = Label(self.master, text="Weather data file")
        lb_weather.grid(row=1, column=2)
        
        bt_simulation = Button(self.master, text="Simulation", width=20, 
                               command=self.simulation_dialog)
        bt_simulation.grid(row=0, column=3)
        lb_simulation = Label(self.master, text="Simulation configuration")
        lb_simulation.grid(row=1, column=3)
        
        bt_exit = Button(self.master, text="Exit", width=20, command=self.site_dialog)
        bt_exit.grid(row=0, column=4)
        lb_exit = Label(self.master, text="Exit")
        lb_exit.grid(row=1, column=4)
        
        bt_operation = Button(self.master, text="Operation Settings", 
                                    width=20, command=self.operation_dialog)
        bt_operation.grid(row=0, column=4)
        lb_operation = Label(self.master, text="Operation Settings")
        lb_operation.grid(row=1, column=4)
#        Button(root, text="cambiar valor", command=self.site_dialog).pack()
#        self.modalwindowtext = StringVar()
#        self.modalwindowtext.set("Entrada de datos")
#        Label(self.master, textvariable=self.modalwindowtext).pack()
 
    def site_dialog(self):
        d = SiteDialog(root, "Site", "Select Site")
        root.wait_window(d.top)
        #self.modalwindowtext.set(d.ejemplo)
        
    def fluid_dialog(self):
        d = FluidDialog(root, "Fluid", "Select Fluid")
        root.wait_window(d.top)
    
    def plant_dialog(self):
        d = PlantDialog(root, "Plant", "Plant Settings")
        root.wait_window(d.top)    
        
    def simulation_dialog(self):
        d = SimulationDialog(root, "Simulation", "Simulation Settings")
        root.wait_window(d.top)
        
    def weather_dialog(self):
        d = WeatherDialog(root, "Weather", "Select file")
        root.wait_window(d.top)
        
    def operation_dialog(self):
        d = OperationalDialog(root, "Operation Mode", "Operation Settings")
        root.wait_window(d.top)
        
class SiteDialog(object):
    def __init__(self, parent, valor, title, labeltext = '' ):
        self.modalwindowtext = valor
 
        self.top = Toplevel(parent)
        self.top.transient(parent)
        self.top.grab_set()
#        if len(title) > 0: self.top.title(title)
#        if len(labeltext) == 0: labeltext = 'Site'
        Label(self.top, text=labeltext).pack()
        self.top.bind("<Return>", self.ok)
        self.e = Entry(self.top, text=valor.get())
        self.e.bind("<Return>", self.ok)
        self.e.bind("<Escape>", self.cancel)
        self.e.pack(padx=15)
        self.e.focus_set()
        b = Button(self.top, text="OK", command=self.ok)
        b.pack(pady=5)
 
    def ok(self, event=None):
        print("You have selected...", self.e.get())
        self.modalwindowtext.set(self.e.get())
        self.top.destroy()
 
    def cancel(self, event=None):
        self.top.destroy()
        
        
class FluidDialog(object):
    def __init__(self, parent, valor, title, labeltext = '' ):
        self.modalwindowtext = valor
 
        self.top = Toplevel(parent)
        self.top.transient(parent)
        self.top.grab_set()
#        if len(title) > 0: self.top.title(title)
#        if len(labeltext) == 0: labeltext = 'Site'
        Label(self.top, text=labeltext).pack()
        self.top.bind("<Return>", self.ok)
        self.e = Entry(self.top, text=valor.get())
        self.e.bind("<Return>", self.ok)
        self.e.bind("<Escape>", self.cancel)
        self.e.pack(padx=15)
        self.e.focus_set()
        b = Button(self.top, text="OK", command=self.ok)
        b.pack(pady=5)
 
    def ok(self, event=None):
        print("You have selected...", self.e.get())
        self.modalwindowtext.set(self.e.get())
        self.top.destroy()
 
    def cancel(self, event=None):
        self.top.destroy()
        
class PlantDialog(object):
    def __init__(self, parent, valor, title, labeltext = '' ):
        self.modalwindowtext = valor
 
        self.top = Toplevel(parent)
        self.top.transient(parent)
        self.top.grab_set()
#        if len(title) > 0: self.top.title(title)
#        if len(labeltext) == 0: labeltext = 'Site'
        Label(self.top, text=labeltext).pack()
        self.top.bind("<Return>", self.ok)
        self.e = Entry(self.top, text=valor.get())
        self.e.bind("<Return>", self.ok)
        self.e.bind("<Escape>", self.cancel)
        self.e.pack(padx=15)
        self.e.focus_set()
        b = Button(self.top, text="OK", command=self.ok)
        b.pack(pady=5)
 
    def ok(self, event=None):
        print("You have selected...", self.e.get())
        self.modalwindowtext.set(self.e.get())
        self.top.destroy()
 
    def cancel(self, event=None):
        self.top.destroy()      

class SimulationDialog(object):
    def __init__(self, parent, valor, title, labeltext = '' ):
        self.modalwindowtext = valor
 
        self.top = Toplevel(parent)
        self.top.transient(parent)
        self.top.grab_set()
        Label(self.top, text=labeltext).pack()
        self.top.bind("<Return>", self.ok)
        self.e = Entry(self.top, text=valor.get())
        self.e.bind("<Return>", self.ok)
        self.e.bind("<Escape>", self.cancel)
        self.e.pack(padx=15)
        self.e.focus_set()
        b = Button(self.top, text="OK", command=self.ok)
        b.pack(pady=5)
 
    def ok(self, event=None):
        print("You have selected...", self.e.get())
        self.modalwindowtext.set(self.e.get())
        self.top.destroy()
 
    def cancel(self, event=None):
        self.top.destroy()       

class WeatherDialog(object):
    def __init__(self, parent, title, labeltext = '' ):
        
        '''resample(self, rule, how=None, axis=0, fill_method=None, closed=None, 
    label=None, convention='start', kind=None, loffset=None, limit=None,
    base=0, on=None, level=None)
'''
        

#        self.filename = askopenfilename(initialdir = ".",
#                               title = "choose your file",
#                               filetypes = (("TMY files","*.tm2"),
#                                            ("TMY files","*.tm3"),
#                                            ("all files","*.*")))
#     
#        if self.filename:
#            print(f"Loading data from: {self.filename}")
#        else:
#            print("No hay nombre de fichero")
            
        weatherdata = pvlib.iotools.tmy.read_tmy2(filename = None)
        


        # #date_rng = pd.date_range(start='1/1/2014',end='31/12/2014',freq='H')
        # weatherdata = pd.read_csv(file, sep=';', decimal=',', index_col=0)
        # weatherdata.index = pd.to_datetime(weatherdata.index)
        # weatherdata = weatherdata.apply(pd.to_numeric, errors='coerce')
        # print(weatherdata)
        # robj = weatherdata.resample('10T').mean()
        # print(robj)
 
        self.top = Toplevel(parent)
        self.top.transient(parent)
        self.top.grab_set()
        Label(self.top, text=title).pack()
        self.top.bind("<Return>", self.ok)
        self.e = Entry(self.top, text=title)
        self.e.bind("<Return>", self.ok)
        self.e.bind("<Escape>", self.cancel)
        self.e.pack(padx=15)
        self.e.focus_set()
        b = Button(self.top, text="OK", command=self.ok)
        b.pack(pady=5)
 
    def ok(self, event=None):
        print("You have selected...", self.e.get())
        self.top.destroy()
 
    def cancel(self, event=None):
        self.top.destroy()
        

class OperationalDialog(object):
    def __init__(self, parent, valor, title, labeltext = '' ):
        self.modalwindowtext = valor
 
        self.top = Toplevel(parent)
        self.top.transient(parent)
        self.top.grab_set()
#        if len(title) > 0: self.top.title(title)
#        if len(labeltext) == 0: labeltext = 'Site'
        Label(self.top, text=labeltext).pack()
        self.top.bind("<Return>", self.ok)
        self.e = Entry(self.top, text=valor.get())
        self.e.bind("<Return>", self.ok)
        self.e.bind("<Escape>", self.cancel)
        self.e.pack(padx=15)
        self.e.focus_set()
        b = Button(self.top, text="OK", command=self.ok)
        b.pack(pady=5)
 
    def ok(self, event=None):
        print("You have selected...", self.e.get())
        self.modalwindowtext.set(self.e.get())
        self.top.destroy()
 
    def cancel(self, event=None):
        self.top.destroy()        
        
        





#wd = Tk()
#wd.title("Data input")
#wd.resizable( True, True)
##wd.iconbitmap("appicon.ico")
#wd.geometry()
#
#
#fm = Frame(wd)
#fm.config(width="650", height="350")
#fm.pack(side="right", anchor="s", fill="both", expand="true")
#
#
##tb = Entry(fm).place()
##tb.pack()
#
#lb = Label(fm, text="Datos").place(x=100, y=200)
#
#btn = Button(fm, text="Click Me")
#btn.grid(column=1, row=0)
#
##tb = Entry(fm)
##tb.place(x=100,y=100)
#
#wd.mainloop()
#    
#wd.destroy