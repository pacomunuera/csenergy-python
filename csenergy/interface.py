# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 11:53:18 2020

@author: fmunuera
"""
from tkinter import * 

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
        Frame.__init__(self, master=None)
        self.master.title("CSENERGY")
        self.master.geometry("800x600")
        
        bt_site = Button(root, text="Site", width=20, 
                         command=self.site_dialog)
        bt_site.grid(row=0, column=0)
        self.modalwindowtext = StringVar()
        self.modalwindowtext.set("Site selection")
        lb_site = Label(self.master, textvariable=self.modalwindowtext)
        lb_site.grid(row=1, column=0)
        
        bt_plant = Button(root, text="Plant", width=20, 
                          command=self.plant_dialog)
        bt_plant.grid(row=0, column=1)
        self.modalwindowtext = StringVar()
        self.modalwindowtext.set("Plant configuration")
        lb_plant = Label(self.master, textvariable=self.modalwindowtext)
        lb_plant.grid(row=1, column=1)
        
        bt_weather = Button(root, text="Weather", width=20, 
                            command=self.weather_dialog)
        bt_weather.grid(row=0, column=2)
        self.modalwindowtext = StringVar()
        self.modalwindowtext.set("Weather data file")
        lb_weather = Label(self.master, textvariable=self.modalwindowtext)
        lb_weather.grid(row=1, column=2)
        
        bt_simulation = Button(root, text="Simulation", width=20, 
                               command=self.simulation_dialog)
        bt_simulation.grid(row=0, column=3)
        self.modalwindowtext = StringVar()
        self.modalwindowtext.set("Simulation configuration")
        lb_simulation = Label(self.master, textvariable=self.modalwindowtext)
        lb_simulation.grid(row=1, column=3)
        
        bt_exit = Button(root, text="Exit", width=20, command=self.site_dialog)
        bt_exit.grid(row=0, column=4)
        self.modalwindowtext = StringVar()
        self.modalwindowtext.set("Exit")
        lb_exit = Label(self.master, textvariable=self.modalwindowtext)
        lb_exit.grid(row=1, column=4)
        
        bt_operationaldata = Button(root, text="Operational Data", 
                                    width=20, command=self.operationaldata_dialog)
        bt_operationaldata.grid(row=0, column=4)
#        self.modalwindowtext = StringVar()
#        self.modalwindowtext.set("Operational Data")
        lb_operationaldata = Label(self.master, textvariable="Operational Data")
        lb_operationaldata.grid(row=1, column=4)
#        Button(root, text="cambiar valor", command=self.site_dialog).pack()
#        self.modalwindowtext = StringVar()
#        self.modalwindowtext.set("Entrada de datos")
#        Label(self.master, textvariable=self.modalwindowtext).pack()
 
    def site_dialog(self):
        d = SiteDialog(root, self.modalwindowtext, "Site", "Select Site")
        root.wait_window(d.top)
        #self.modalwindowtext.set(d.ejemplo)
        
    def fluid_dialog(self):
        d = FluidDialog(root, self.modalwindowtext, "Fluid", "Select Fluid")
        root.wait_window(d.top)
    
    def plant_dialog(self):
        d = PlantDialog(root, self.modalwindowtext, "Plant", "Plant Settings")
        root.wait_window(d.top)    
        
    def simulation_dialog(self):
        d = SimulationDialog(root, self.modalwindowtext, "Simulation", "Simulation Settings")
        root.wait_window(d.top)
        
    def weather_dialog(self):
        d = WeatherDialog(root, self.modalwindowtext, "Weather", "Select file")
        root.wait_window(d.top)
        
    def operationaldata_dialog(self):
        d = OperationalDataDialog(root, self.modalwindowtext, "Weather", "Select file")
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

class WeatherDialog(object):
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
        

class OperationalDataDialog(object):
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
        
        

root = Tk()
a = main(root)
root.mainloop()



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