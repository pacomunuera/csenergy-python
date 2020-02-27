import csenergy as cs
import pandas as pd
from tkinter import *
import json

with open("./saved_configurations/simulation1x1.json") as simulation_file:
    simulation_settings = json.load(simulation_file)

simulation = cs.Simulation(simulation_settings['simulation'])

if simulation.type == "type0":
    datasource = cs.Weather(simulation_settings['weather_data_file'])
    print("Simulation based on weather data")
elif simulation.type == "type1":
    datasource = cs.FieldData(simulation_settings['field_data_file'])
    print("Benchmarking based on actual data")

if datasource.site is None:
    site = cs.Site(simulation_settings['site'])
else:
    site = cs.Site(datasource.site_to_dict())

coolPropFluids = ['Water', 'INCOMP::TVP1', 'INCOMP::S800']

if simulation_settings['hot_fluid']['CoolPropID'] in coolPropFluids:
    hotfluid = cs.Fluid_CoolProp(simulation_settings['hot_fluid'])
    print("Fluid data from CoolProp: ", hotfluid.name)
else:
    hotfluid = cs.Fluid_Tabular(simulation_settings['hot_fluid'])
    print("Fluid data from table: ", hotfluid.name )

if simulation_settings['cold_fluid']['CoolPropID'] in coolPropFluids:
    coldfluid = cs.Fluid_CoolProp(simulation_settings['cold_fluid'])
    print("Cold fluid data from CoolProp: ", coldfluid.name)
else:
    coldfluid = cs.Fluid_Tabular(simulation_settings['cold_fluid'])
    print("Cold fluid data from table: ", coldfluid.name)

solarplant = cs.SolarPlant(simulation_settings['solar_plant'],
                           simulation_settings['sca'],
                           simulation_settings['hce'],
                           simulation_settings['hce_model_settings'])

hcemask = cs.HCEScatterMask(simulation_settings['solar_plant'],
                            simulation_settings['hce_scattered_params'])
scamask = cs.SCAScatterMask(simulation_settings['solar_plant'],
                            simulation_settings['sca_scattered_params'])

# TO-DO: MÁSCARAS PARA INTRODUCIR DISPERSIÓN EN LOS PARÁMETROS.
#hcemask.applyMask(solarplant)
#scamask.applyMask(solarplant)

powercycle = cs.PowerCycle(simulation_settings['powercycle'])

heatexchanger = cs.HeatExchanger(simulation_settings['heatexchanger'],
                                 hotfluid, coldfluid)

generator = cs.Generator(simulation_settings['generator'])

model = cs.ModelBarbero4grade(simulation_settings['hce_model_settings'])

powersystem = cs.PowerSystem(simulation_settings['powersystem'])

#simulation.precalc(powersystem, solarplant, hotfluid, simulation,
#                   simulation_settings['hce'])

prototypeloop = cs.PrototypeLoop(simulation_settings['solar_plant'],
                                 simulation_settings['sca'],
                                 simulation_settings['hce'],
                                 simulation_settings['hce_model_settings']
                                 )

simulation.hotfluid = hotfluid
simulation.coldfluid = coldfluid
simulation.site = site
simulation.solarplant = solarplant
simulation.model = model
simulation.datasource = datasource
simulation.prototypeloop = prototypeloop
simulation.powercycle = powercycle
simulation.heatexchanger = heatexchanger
simulation.generator = generator
simulation.powersystem = powersystem

simulation.runSimulation()






