# -*- coding: utf-8 -*-
import csenergy as cs
import pandas as pd
from tkinter import *
import json
import copy
from datetime import datetime

with open("./saved_configurations/simulation4x36.json") as simulation_file:
    simulation_settings = json.load(simulation_file)

simulation = cs.Simulation(simulation_settings['simulation'])


if simulation.type == "type0":
    datasource = cs.Weather(simulation_settings['weather_data_file'])
    print("Simulation based on weather data")
elif simulation.type == "type1":
    datasource = cs.FieldData(simulation_settings['field_data_file'])
    print("Benchmarking based on actual data")


if not hasattr(datasource, 'site'):
    site = cs.Site(simulation_settings['site'])
else:
    site = cs.Site(datasource.site_to_dict())

coolPropFluids = ['Water', 'INCOMP::TVP1', 'INCOMP::S800']

if simulation_settings['hot_fluid']['CoolPropID'] in coolPropFluids:
    hotfluid = cs.Fluid_CoolProp(simulation_settings['hot_fluid'])
    print("Fluid data from CoolProp: ", hotfluid.name)
else:
    hotfluid = cs.Fluid_Tabular(simulation_settings['hot_fluid'])
    print("Fluid data from table: ", hotfluid.name)

if simulation_settings['cold_fluid']['CoolPropID'] in coolPropFluids:
    coldfluid = cs.Fluid_CoolProp(simulation_settings['cold_fluid'])
    print("Cold fluid data from CoolProp: ", coldfluid.name)
else:
    coldfluid = cs.Fluid_Tabular(simulation_settings['cold_fluid'])
    print("Cold fluid data from table: ", coldfluid.name)




solarfield = cs.SolarField(simulation_settings['solarfield'],
                           simulation_settings['sca'],
                           simulation_settings['hce'],
                           simulation_settings['hce_model_settings'])



model = cs.ModelBarbero4grade(simulation_settings['hce_model_settings'])

prototype_loop = cs.PrototypeLoop(solarfield)

simulation.hotfluid = hotfluid
simulation.coldfluid = coldfluid
simulation.site = site
simulation.solarfield = solarfield
simulation.model = model
simulation.datasource = datasource
simulation.prototype_loop = prototype_loop


flag_00 = datetime.now()

for row in simulation.datasource.dataframe.iterrows():

    solarpos = simulation.get_solarposition(row)
    aoi =

    (self, aoi, solarpos, row):
    values = {'tin': 563,
              'pin': 1900000,
              'massflow': 4}

    prototype_loop.initialize(values)

    print(prototype_loop.tin, prototype_loop.pin, prototype_loop.massflow)

    prototype_loop.calc_loop_pr_for_tout(row, solarpos, hotfluid, model)

    simulation.datasource.dataframe.at[row[0], 'PTLoop_mf'] = prototype_loop.massflow
    simulation.datasource.dataframe.at[row[0], 'PTLoop_tin'] = prototype_loop.tin
    simulation.datasource.dataframe.at[row[0], 'PTLoop_tout'] = prototype_loop.tout


flag_01 = datetime.now()
delta_01 = flag_01 - flag_00
print("Total runtime: ", delta_01.total_seconds())
