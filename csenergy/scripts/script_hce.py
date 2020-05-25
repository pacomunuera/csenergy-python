# -*- coding: utf-8 -*-
import csenergy as cs
import pandas as pd
from tkinter import *
import json
import copy
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import rc


with open("./saved_configurations/simulation_S_F_UVAC.json") as simulation_file:
    simulation_settings = json.load(simulation_file)


hce_settings = simulation_settings['HCE']
sca_settings = simulation_settings['SCA']
loop_settings = simulation_settings['loop']
model_settings = simulation_settings['model']

site = cs.Site(simulation_settings['site'])

loop = cs.BaseLoop(loop_settings, sca_settings, hce_settings)
sca = cs.SCA(loop, 0, sca_settings)
hce = cs.HCE(sca, 0, hce_settings)

htf = cs.FluidTabular(simulation_settings['HTF'])

model = cs.ModelBarbero4thOrder(model_settings)

data_file = {'filename': 'serie_qabs_corto.csv', 'filepath': './data_files/'}

datasource = cs.TableData(data_file)

site = cs.Site(simulation_settings['site'])

row = next(datasource.dataframe.iterrows())[1]

row2 = tuple([pd.to_datetime(row[0]), row])

print(row2[0], '\n')
print(row2[1])

flag_00 = datetime.now()

results = pd.DataFrame(columns=['qabs', 'pr', 'tout', 'pout'])

for t in range(350, 673, 5):

    hce.tin = t
    hce.pin = 1900000
    hce.sca.loop.massflow = 4.0

    # solarpos = site.get_solarposition(row)
    # aoi = sca.get_aoi(solarpos)

    # qabs = row[1]['qabs']
    qabs = 15000
    aoi = 15

    model.calc_pr(hce, htf, qabs, row2)

    print(hce.pr, hce.tout, qabs)

    results.at[t, 'qabs'] = qabs
    results.at[t, 'pr'] = hce.pr
    results.at[t, 'tout'] = hce.tout
    results.at[t, 'pout'] = hce.pout


results[['pr']].plot(figsize=(20,10), linewidth=5, fontsize=20)

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
## rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

plt.ylabel('Performance, pr', fontsize=20)
plt.xlabel('$t{in} [K]$', fontsize=20)

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

flag_01 = datetime.now()
delta_01 = flag_01 - flag_00
print("Total runtime: ", delta_01.total_seconds())
