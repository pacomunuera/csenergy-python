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
# from matplotlib import rc




with open("./saved_configurations/test_4.json") as simulation_file:
    simulation_settings = json.load(simulation_file)


modelos = ['Barbero4thOrder', 'Barbero1stOrder', 'BarberoSimplified']

dict_resultados = {}

for m in modelos:

    qabs_index = []

    simulation_settings['model']['name'] = m
    simulation = cs.LoopSimulation(simulation_settings)
    dict_resultados[m] = []
    for row in simulation.datasource.dataframe.iterrows():

        solarpos = simulation.site.get_solarposition(row)
        aoi = simulation.base_loop.scas[0].get_aoi(solarpos)
        values = {'tin': 573,
                  'pin': 1900000,
                  'massflow': 6}

        simulation.base_loop.initialize('values', values)

        for s in simulation.base_loop.scas:
            aoi = s.get_aoi(solarpos)
            for h in s.hces:
                h.set_pr_opt(solarpos)
                tro = h.tin + row[1]['qabs'] * \
                    h.pr / h.get_urec(h.tin, h.pin, simulation.htf)
                h.qabs = row[1]['qabs']
                h.set_tin()
                h.set_pin()
                simulation.model.calc_pr(h, simulation.htf, row, row[1]['qabs'])
                print(h.get_urec(h.tin, h.pin, simulation.htf))
                print(h.get_hint(h.tin, h.pin, simulation.htf))
                print(h.get_eext(tro,row[1]['Wspd']))

        tout = simulation.base_loop.scas[-1].hces[-1].tout
        pout = simulation.base_loop.scas[-1].hces[-1].pout
        simulation.base_loop.set_loop_values_from_HCEs()
        pr = simulation.base_loop.pr
        dict_resultados[m].append(pr)
        qabs_index.append(row[1]['qabs'])

        simulation.datasource.dataframe.at[row[0], 'pr'] = pr
        simulation.datasource.dataframe.at[row[0], 'tout'] = tout
        simulation.datasource.dataframe.at[row[0], 'pout'] = pout

dfsalida = pd.DataFrame(dict_resultados, index=qabs_index)

print(dfsalida)

dfsalida.plot(figsize=(20,10), linewidth=5, fontsize=20)

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})

rc('text', usetex=True)

plt.ylabel('Performance, pr', fontsize=20)
plt.xlabel('HCE', fontsize=20)

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

