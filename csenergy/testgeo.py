import csenergy as cs
import pandas as pd
from tkinter import *
import json

with open("./saved_configurations/geo.json") as simulation_file:
    simulation_settings = json.load(simulation_file)

simulation = cs.Simulation(simulation_settings['simulation'])

if simulation.type == "type0":
    datasource = cs.Weather(simulation_settings['weather_data_file'])
    print("Simulation based on weather data")
elif simulation.type == "type1":
    datasource = cs.FieldData(simulation_settings['field_data_file'])
    print("Benchmarking based on actual data")

site = cs.Site(simulation_settings['site'])


solarplant = cs.SolarPlant(simulation_settings['solar_plant'],
                           simulation_settings['sca'],
                           simulation_settings['hce'],
                           simulation_settings['hce_model_settings'])

prototypeloop = cs.PrototypeLoop(simulation_settings['solar_plant'],
                                 simulation_settings['sca'],
                                 simulation_settings['hce'],
                                 simulation_settings['hce_model_settings']
                                 )

simulation.hotfluid = hotfluid
simulation.coldfluid = coldfluid
simulation.site = site
simulation.solarplant = solarplant
simulation.datasource = datasource
simulation.prototypeloop = prototypeloop

simulation.testgeo()






