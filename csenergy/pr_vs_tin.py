# -*- coding: utf-8 -*-
"""
@author: paco munuera
2020
"""

import csenergy as cs
import pandas as pd
import json
import matplotlib.pyplot as plt

with open("./saved_configurations/test_2.json") as simulation_file:
    simulation_settings = json.load(simulation_file)

with open("./hce_files/HCE_library.json") as hces_file:
    hces_configurations = json.load(hces_file)

dict_resultados = {}

for hce_conf in hces_configurations:

    simulation_settings['HCE'].update(hce_conf)
    tin_index = []
    simulation = cs.LoopSimulation(simulation_settings)
    dict_resultados[hce_conf['Name']] = []

    for row in simulation.datasource.dataframe.iterrows():

        solarpos = simulation.site.get_solarposition(row)
        aoi = simulation.base_loop.scas[0].get_aoi(solarpos)

        for tin in range(473, 693, 10):
            values = {'tin': tin,
                      'pin': 1900000,
                      'massflow': 6}
            simulation.base_loop.initialize('values', values)
            simulation.base_loop.calc_loop_pr_for_massflow(
                row,
                solarpos,
                simulation.htf,
                simulation.model)
            simulation.base_loop.set_loop_values_from_HCEs()
            pr = simulation.base_loop.pr
            dict_resultados[hce_conf['Name']].append(pr)
            tin_index.append(tin-273)
            simulation.datasource.dataframe.at[row[0], 'pr'] = pr

dfsalida = pd.DataFrame(dict_resultados, index=tin_index)

print(dfsalida)

grafica = dfsalida.plot()
plt.ylabel('Rendimiento', fontsize=20)
plt.xlabel('Temperatura de entrada (°C)', fontsize=20)



