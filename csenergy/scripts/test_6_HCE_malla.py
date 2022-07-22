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




with open("./saved_configurations/test_6_malla.json") as simulation_file:
    simulation_settings = json.load(simulation_file)


with open("./hce_files/HCE_library_malla.json") as hces_file:
    hces_configurations = json.load(hces_file)

dict_resultados = {}


modelos = ['Barbero4thOrder', 'Barbero1stOrder', 'BarberoSimplified']

for m in modelos:

    simulation_settings['model']['name'] = m

    for hce_conf in hces_configurations:

        simulation_settings['HCE'].update(hce_conf)
        simulation_settings['loop']['hces'] = round(
            simulation_settings['SCA']['SCA Length'] /
            simulation_settings['HCE']['Length'])

        dni_index = []

        simulation = cs.LoopSimulation(simulation_settings)
        dict_resultados[hce_conf['Length']] = []
        for row in simulation.datasource.dataframe.iterrows():

            solarpos = simulation.site.get_solarposition(row)
            aoi = simulation.base_loop.scas[0].get_aoi(solarpos)
            values = {'tin': 573,
                      'pin': 1900000,
                      'massflow': 6}

            simulation.base_loop.initialize('values', values)

            simulation.base_loop.calc_loop_pr_for_massflow(
                row,
                solarpos,
                simulation.htf,
                simulation.model)

            tout = simulation.base_loop.scas[-1].hces[-1].tout
            pout = simulation.base_loop.scas[-1].hces[-1].pout
            simulation.base_loop.set_loop_values_from_HCEs()
            pr = simulation.base_loop.pr
            dict_resultados[hce_conf['Length']].append(pr)
            dni_index.append(row[1]['DNI'])

            simulation.datasource.dataframe.at[row[0], 'pr'] = pr
            simulation.datasource.dataframe.at[row[0], 'tout'] = tout
            simulation.datasource.dataframe.at[row[0], 'pout'] = pout

    dfsalida = pd.DataFrame(dict_resultados, index=dni_index)

    print(dfsalida)
    dfsalida.to_csv(m, sep=';', decimal=',')

# dfsalida.plot(figsize=(20,10), linewidth=5, fontsize=20)

# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# # for Palatino and other serif fonts use:
# # rc('font',**{'family':'serif','serif':['Palatino']})

# rc('text', usetex=True)

# plt.ylabel('Performance, pr', fontsize=20)
# plt.xlabel('HCE', fontsize=20)

# pd.set_option('display.max_rows', None)
# pd.set_option('display.max_columns', None)
# pd.set_option('display.width', None)
