# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 11:29:44 2020

@author: paco
"""

# -*- coding: utf-8 -*-
import csenergy as cs
import pandas as pd
import json
from datetime import datetime

FLAG_00 = datetime.now()
with open("./saved_configurations/test_1.json") as simulation_file:
    simulation_settings = json.load(simulation_file)


with open("./hce_files/HCE_library.json") as hces_file:
    hces_configurations = json.load(hces_file)

dict_resultados = {}

for hce_conf in hces_configurations:

    simulation_settings['HCE'].update(hce_conf)
    dni_index = []

    simulation = cs.LoopSimulation(simulation_settings)
    dict_resultados[hce_conf['Name']] = []
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
        dict_resultados[hce_conf['Name']].append(pr)
        dni_index.append(row[1]['DNI'])

        simulation.datasource.dataframe.at[row[0], 'pr'] = pr
        simulation.datasource.dataframe.at[row[0], 'tout'] = tout
        simulation.datasource.dataframe.at[row[0], 'pout'] = pout

dfsalida = pd.DataFrame(dict_resultados, index=dni_index)

dfsalida.to_csv('rendimiento_dni.csv', sep=';', decimal=',')
FLAG_01 = datetime.now()
DELTA_01 = FLAG_01 - FLAG_00
print("Total runtime: ", DELTA_01.total_seconds())


