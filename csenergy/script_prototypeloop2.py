# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 11:29:44 2020

@author: paco
"""

# -*- coding: utf-8 -*-
import csenergy as cs
import pandas as pd
from tkinter import *
import json
import copy
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import rc




with open("./saved_configurations/simulation_loop.json") as simulation_file:
    simulation_settings = json.load(simulation_file)

simulation = cs.Simulation(simulation_settings['simulation'])


data_file = {'filename': 'serie_aperture.csv', 'filepath': './data_files/'}
datasource = cs.TableData(data_file)

site = cs.Site(simulation_settings['site'])


coolPropFluids = ['Water', 'INCOMP::TVP1', 'INCOMP::S800']

hotfluid = cs.FluidTabular(simulation_settings['HTF'])
print("Fluid data from table: ", hotfluid.name)

model = cs.ModelBarbero4thOrder()

prototype_loop = cs.PrototypeLoop(simulation_settings['loop'])

simulation.hotfluid = hotfluid
simulation.site = site
simulation.model = model
simulation.datasource = datasource
simulation.prototype_loop = prototype_loop


flag_00 = datetime.now()

print('row spacing', prototype_loop.row_spacing)

for row in simulation.datasource.dataframe.iterrows():

    solarpos = simulation.get_solarposition(row)
    aoi = prototype_loop.scas[0].get_aoi(solarpos)
    values = {'tin': 573,
              'pin': 1900000,
              'massflow': 4}

    prototype_loop.initialize(values)

    for s in prototype_loop.scas:
            aoi = s.get_aoi(solarpos)
            for h in s.hces:
                h.sca.parameters['Aperture'] = row[1]['aperture']
                qabs = h.get_qabs(aoi, solarpos, row)
                model.calc_pr(h, hotfluid, qabs, row)

    tout = prototype_loop.scas[-1].hces[-1].tout
    pout = prototype_loop.scas[-1].hces[-1].pout
    pr = prototype_loop.get_loop_avg_pr()
    print('pr', pr, 'aperture', row[1]['aperture'])

    simulation.datasource.dataframe.at[row[0], 'aperture'] = row[1]['aperture']
    simulation.datasource.dataframe.at[row[0], 'pr'] = pr
    simulation.datasource.dataframe.at[row[0], 'tout'] = tout
    simulation.datasource.dataframe.at[row[0], 'pout'] = pout



dfsalida = simulation.datasource.dataframe[
    simulation.datasource.dataframe['pr']>0.00]

print(dfsalida)

dfsalida[['pr']].plot(figsize=(20,10), linewidth=5, fontsize=20)

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
## rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

plt.ylabel('Performance, pr', fontsize=20)
plt.xlabel('$q_{abs} [kW]$', fontsize=20)

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

flag_01 = datetime.now()
delta_01 = flag_01 - flag_00
print("Total runtime: ", delta_01.total_seconds())
