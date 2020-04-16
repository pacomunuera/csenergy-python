import csenergy as cs
import pandas as pd
from tkinter import *
import json
import copy
from datetime import datetime

with open("./saved_configurations/simulation4x36_4.json") as simulation_file:
    simulation_settings = json.load(simulation_file)

simulation = cs.Simulation(simulation_settings['simulation'])

if simulation.datatype == 1:
    datasource = cs.Weather(simulation_settings['simulation'])
    print("Simulation based on weather data")
elif simulation.datatype == 2:
    datasource = cs.FieldData(simulation_settings['simulation'])
    print("Benchmarking based on actual data")


if not hasattr(datasource, 'site'):
    site = cs.Site(simulation_settings['site'])
else:
    site = cs.Site(datasource.site_to_dict())

coolPropFluids = ['Water', 'INCOMP::TVP1', 'INCOMP::S800']

if simulation_settings['HTF']['source'] == "CoolProp":
    if simulation_settings['HTF']['CoolPropID'] not in coolPropFluids:
        print("Not CoolPropID valid")
    else:
        htf = cs.Fluid_CoolProp(simulation_settings['HTF'])
        print("Fluid data from CoolProp: ", htf.name)
else:
    htf = cs.Fluid_Tabular(simulation_settings['HTF'])
    print("Fluid data from table: ", htf.name)

# if simulation_settings['cold_fluid']['CoolPropID'] in coolPropFluids:
#     coldfluid = cs.Fluid_CoolProp(simulation_settings['cold_fluid'])
#     print("Cold fluid data from CoolProp: ", coldfluid.name)
# else:
#     coldfluid = cs.Fluid_Tabular(simulation_settings['cold_fluid'])
#     print("Cold fluid data from table: ", coldfluid.name)

solarfield = cs.SolarField(simulation_settings['solarfield'],
                           simulation_settings['SCA'],
                           simulation_settings['HCE'])

# hcemask = cs.HCEScatterMask(simulation_settings['solarfield'],
#                             simulation_settings['hce_scattered_params'])
# scamask = cs.SCAScatterMask(simulation_settings['solarfield'],
#                             simulation_settings['sca_scattered_params'])

# TO-DO: MÁSCARAS PARA INTRODUCIR DISPERSIÓN EN LOS PARÁMETROS.
#hcemask.applyMask(solarfield)
#scamask.applyMask(solarfield)

# powercycle = cs.PowerCycle(simulation_settings['powercycle'])

# heatexchanger = cs.HeatExchanger(simulation_settings['heatexchanger'],
#                                  htf, coldfluid)

# generator = cs.Generator(simulation_settings['generator'])

model = cs.ModelBarbero4grade()

# powersystem = cs.PowerSystem(simulation_settings['powersystem'])

#simulation.precalc(powersystem, solarfield, htf, simulation,
#                   simulation_settings['hce'])

prototype_loop = cs.PrototypeLoop(solarfield)

simulation.htf = htf
# simulation.coldfluid = coldfluid
simulation.site = site
simulation.solarfield = solarfield
simulation.model = model
simulation.datasource = datasource
simulation.prototype_loop = prototype_loop
# simulation.powercycle = powercycle
# simulation.heatexchanger = heatexchanger
# simulation.generator = generator
# simulation.powersystem = powersystem

flag_00 = datetime.now()
simulation.runSimulation()
flag_01 = datetime.now()
delta_01 = flag_01 - flag_00
print("Total runtime: ", delta_01.total_seconds())






