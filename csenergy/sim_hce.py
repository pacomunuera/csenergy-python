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

if simulation.datatype == "weather":
    datasource = cs.Weather(simulation_settings['weather_data_file'])
    print("Simulation based on weather data")
elif simulation.datatype == "field data":
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

    aoi = prototype_loop.scas[0].get_aoi(solarpos)

    values = {'tin': 563,
              'pin': 1900000,
              'massflow': 4}

    prototype_loop.initialize(values)

    #print(prototype_loop.tin, prototype_loop.pin, prototype_loop.massflow)

    prototype_loop.calc_loop_pr_for_tout(row, solarpos, hotfluid, model)

    pr_geo = prototype_loop.scas[0].hces[0].get_pr_geo(aoi, solarpos, row)
    pr_opt_peak =  prototype_loop.scas[0].hces[0].get_pr_opt_peak(aoi, solarpos, row)
    pr_shadows = prototype_loop.scas[0].hces[0].get_pr_shadows(aoi, solarpos, row)
    solarfraction = prototype_loop.scas[0].get_solar_fraction(aoi, solarpos, row)
    iam = prototype_loop.scas[0].get_IAM(aoi)

    simulation.datasource.dataframe.at[row[0], 'pr_geo'] = pr_geo
    simulation.datasource.dataframe.at[row[0], 'pr_opt_peak'] = pr_opt_peak
    simulation.datasource.dataframe.at[row[0], 'pr_shadows'] = pr_shadows
    simulation.datasource.dataframe.at[row[0], 'solar_fraction'] = solarfraction
    simulation.datasource.dataframe.at[row[0], 'elevation'] = solarpos['elevation'][0]
    simulation.datasource.dataframe.at[row[0], 'azimuth'] = solarpos['azimuth'][0]
    simulation.datasource.dataframe.at[row[0], 'iam'] = iam
    simulation.datasource.dataframe.at[row[0], 'aoi'] = aoi




flag_01 = datetime.now()
delta_01 = flag_01 - flag_00
print("Total runtime: ", delta_01.total_seconds())
