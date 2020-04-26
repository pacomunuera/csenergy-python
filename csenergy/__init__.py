import csenergy as cs
import pandas as pd
from tkinter import *
import json
import copy
import sys
from datetime import datetime

with open("./saved_configurations/simulation_B_F_UVAC.json") as simulation_file:
    simulation_settings = json.load(simulation_file)

simulation = cs.Simulation(simulation_settings['simulation'])

if simulation.datatype == 1:
    datasource = cs.Weather(simulation_settings['simulation'])
elif simulation.datatype == 2:
    datasource = cs.FieldData(simulation_settings['simulation'],
                              simulation_settings['tags'])

if not hasattr(datasource, 'site'):
    site = cs.Site(simulation_settings['site'])
else:
    site = cs.Site(datasource.site_to_dict())

coolPropFluids = ['Water', 'INCOMP::TVP1', 'INCOMP::S800']

if simulation_settings['HTF']['source'] == "CoolProp":
    if simulation_settings['HTF']['CoolPropID'] not in coolPropFluids:
        print("Not CoolPropID valid")
        sys.exit()
    else:
        htf = cs.FluidCoolProp(simulation_settings['HTF'])

else:
    htf = cs.FluidTabular(simulation_settings['HTF'])

solarfield = cs.SolarField(simulation_settings['subfields'],
                           simulation_settings['loop'],
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

model = cs.ModelBarbero4thOrder()

# powersystem = cs.PowerSystem(simulation_settings['powersystem'])

#simulation.precalc(powersystem, solarfield, htf, simulation,
#                   simulation_settings['hce'])

base_loop = cs.BaseLoop(simulation_settings['loop'],
                        simulation_settings['SCA'],
                        simulation_settings['HCE'])

simulation.htf = htf
# simulation.coldfluid = coldfluid
simulation.site = site
simulation.solarfield = solarfield
simulation.model = model
simulation.datasource = datasource
simulation.base_loop = base_loop
# simulation.powercycle = powercycle
# simulation.heatexchanger = heatexchanger
# simulation.generator = generator
# simulation.powersystem = powersystem

flag_00 = datetime.now()
print("Running simulation for source data file: {0}".format(
    simulation_settings['simulation']['filename']))
print("Simulation: {0} ; Benchmark: {1} ; FastMode: {2}".format(
    simulation_settings['simulation']['simulation'],
    simulation_settings['simulation']['benchmark'],
    simulation_settings['simulation']['fastmode']))

print("Site: {0} @ Lat: {1:.2f}º, Long: {2:.2f}º, Alt: {3} m".format(
    site.name, site.latitude, site.longitude, site.altitude))

print("Loops:", solarfield.total_loops,
      'SCA/loop:', simulation_settings['loop']['scas'],
      'HCE/SCA:', simulation_settings['loop']['hces'])
print("SCA model:", simulation_settings['SCA']['Name'])
print("HCE model:", simulation_settings['HCE']['Name'])
print("HTF:", simulation_settings['HTF']['name'])
print("---------------------------------------------------")
simulation.runSimulation()
flag_01 = datetime.now()
delta_01 = flag_01 - flag_00
print("Total runtime: ", delta_01.total_seconds())






