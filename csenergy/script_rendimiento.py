# -*- coding: utf-8 -*-
import csenergy as cs
import pandas as pd
from tkinter import *
import json
import copy
from datetime import datetime
import matplotlib.pyplot as plt

with open("./saved_configurations/simulation_S_F_UVAC.json") as simulation_file:
    simulation_settings = json.load(simulation_file)

simulation = cs.Simulation(simulation_settings['simulation'])


data_file = {'filename': 'serie_qabs_corto.csv', 'filepath': 'data_files/'}
datasource = cs.TableData(data_file)

site = cs.Site(simulation_settings['site'])


coolPropFluids = ['Water', 'INCOMP::TVP1', 'INCOMP::S800']

hotfluid = cs.FluidTabular(simulation_settings['HTF'])
print("Fluid data from table: ", hotfluid.name)

solarfield = cs.SolarField(simulation_settings['solarfield'],
                           simulation_settings['SCA'],
                           simulation_settings['HCE'])


model = cs.ModelBarberoSimplified()

prototype_loop = cs.PrototypeLoop(solarfield)

simulation.hotfluid = hotfluid
simulation.site = site
simulation.solarfield = solarfield
simulation.model = model
simulation.datasource = datasource
simulation.prototype_loop = prototype_loop


flag_00 = datetime.now()

for row in simulation.datasource.dataframe.iterrows():

    solarpos = simulation.get_solarposition(row)

    aoi = prototype_loop.scas[0].get_aoi(solarpos)

    values = {'tin': 573,
              'pin': 1900000,
              'massflow': 4}

    prototype_loop.initialize(values)

    #print(prototype_loop.tin, prototype_loop.pin, prototype_loop.massflow)
    for s in prototype_loop.scas:
            aoi = s.get_aoi(solarpos)
            for h in s.hces:
                qabs = row[1]['qabs']
                model.calc_pr(h, hotfluid, qabs, row)

    tout = prototype_loop.scas[-1].hces[-1].tout
    pout = prototype_loop.scas[-1].hces[-1].pout
    pr = prototype_loop.get_loop_avg_pr()

    simulation.datasource.dataframe.at[row[0], 'qabs'] = qabs
    simulation.datasource.dataframe.at[row[0], 'pr'] = pr
    simulation.datasource.dataframe.at[row[0], 'tout'] = tout
    simulation.datasource.dataframe.at[row[0], 'pout'] = pout


dfsalida = simulation.datasource.dataframe[
    simulation.datasource.dataframe['pr']>0.00]

# dfsalida = pd.DataFrame([[simulation.datasource.dataframe['pr'] > 0.94]])

print(dfsalida)

dfsalida[['pr']].plot(figsize=(20,10), linewidth=5, fontsize=20)
plt.xlabel('qabs', fontsize=20)

# simulation.datasource.dataframe[['pr']].plot(figsize=(20,10), linewidth=5, fontsize=20)
# plt.xlabel('qabs', fontsize=20)


# simulation.datasource.dataframe[[
#      simulation.datasource.dataframe['pr']> 0.94]].plot(figsize=(20,10), linewidth=5, fontsize=20)
# plt.xlabel('qabs', fontsize=20)

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
    # pd.set_option('display.max_colwidth', -1)


flag_01 = datetime.now()
delta_01 = flag_01 - flag_00
print("Total runtime: ", delta_01.total_seconds())
