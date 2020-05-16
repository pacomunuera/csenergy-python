# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 11:29:44 2020

@author: paco
"""
# -*- coding: utf-8 -*-
import csenergy as cs
import numpy as np
import pandas as pd
from tkinter import *
import json
import copy
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import rc

def save_results(df):


    try:
        initialdir = "./simulations_outputs/"
        prefix = datetime.today().strftime("%Y%m%d %H%M%S")
        filename = " CALCULO PERDIDAS"
        sufix = ".csv"

        path = initialdir + prefix + filename + sufix

        df.to_csv(path, sep=';', decimal = ',')

    except Exception:
        raise
        print('Error saving results, unable to save file: %r', path)


with open("./saved_configurations/perdidas.json") as simulation_file:
    settings = json.load(simulation_file)

datasource = cs.TableData(settings['simulation'])
htf  = cs.FluidTabular(settings['HTF'])


flag_0 = datetime.now()

for row in datasource.dataframe.iterrows():

    SF_MF = row[1]['NO.mf'] + row[1]['NE.mf'] +row[1]['SO.mf'] + row[1]['SE.mf']
    SF_P_IN = np.min([100000*row[1]['NO.pin'],
                     100000*row[1]['NE.pin'],
                     100000*row[1]['SO.pin'],
                     100000*row[1]['SE.pin']])
    SF_P_OUT = np.min([100000*row[1]['NO.pout'],
                     100000*row[1]['NE.pout'],
                     100000*row[1]['SO.pout'],
                     100000*row[1]['SE.pout']])

    NO_H_IN = row[1]['NO.mf'] * htf.get_deltaH(
        row[1]['NO.tin']+273.15, 100000*row[1]['NO.pin'])
    NE_H_IN = row[1]['NE.mf'] * htf.get_deltaH(
        row[1]['NE.tin']+273.15, 100000*row[1]['NE.pin'])
    SO_H_IN = row[1]['SO.mf'] * htf.get_deltaH(
        row[1]['SO.tin']+273.15, 100000*row[1]['SO.pin'])
    SE_H_IN = row[1]['SE.mf'] * htf.get_deltaH(
        row[1]['SE.tin']+273.15, 100000*row[1]['SE.pin'])

    SF_H_IN = (NO_H_IN + NE_H_IN + SO_H_IN + SE_H_IN) / SF_MF
    SF_T_IN = htf.get_T(SF_H_IN, SF_P_IN)


    NO_H_OUT = row[1]['NO.mf'] * htf.get_deltaH(
        row[1]['NO.tout']+273.15, 100000*row[1]['NO.pout'])
    NE_H_OUT = row[1]['NE.mf'] * htf.get_deltaH(
        row[1]['NE.tout']+273.15, 100000*row[1]['NE.pout'])
    SO_H_OUT = row[1]['SO.mf'] * htf.get_deltaH(
        row[1]['SO.tout']+273.15, 100000*row[1]['SO.pout'])
    SE_H_OUT = row[1]['SE.mf'] * htf.get_deltaH(
        row[1]['SE.tout']+273.15, 100000*row[1]['SE.pout'])

    SF_H_OUT = (NO_H_OUT + NE_H_OUT + SO_H_OUT + SE_H_OUT) / SF_MF
    SF_T_OUT = htf.get_T(SF_H_OUT, SF_P_OUT)

    SF_POWER = SF_MF * (SF_H_OUT - SF_H_IN)

    IDP_H_IN = htf.get_deltaH(row[1]['IDP.tin']+273.15, SF_P_OUT)
    IDP_H_OUT = htf.get_deltaH(row[1]['IDP.tout']+273.15, SF_P_IN)

    COLD_PIPE_LOSSES = (SF_H_IN - IDP_H_OUT) * SF_MF
    HOT_PIPE_LOSSES = (IDP_H_IN - SF_H_OUT) * SF_MF

    IDP_POWER = SF_MF * (IDP_H_IN - IDP_H_OUT)

    datasource.dataframe.at[row[0], 'SF_POWER'] = SF_POWER / 1000000
    datasource.dataframe.at[row[0], 'SF_MF'] = SF_MF
    datasource.dataframe.at[row[0], 'SF_T_IN'] = SF_T_IN
    datasource.dataframe.at[row[0], 'SF_T_OUT'] = SF_T_OUT
    datasource.dataframe.at[row[0], 'COLD_PIPE_LOSSES'] =  \
        COLD_PIPE_LOSSES / 1000000
    datasource.dataframe.at[row[0], 'IDP_POWER'] =  \
        IDP_POWER / 1000000
    datasource.dataframe.at[row[0], 'IDP_T_IN'] = row[1]['IDP.tin'] +273.15
    datasource.dataframe.at[row[0], 'IDP_T_OUT'] = row[1]['IDP.tout'] +273.15
    datasource.dataframe.at[row[0], 'HOT_PIPE_LOSSES'] =  \
        HOT_PIPE_LOSSES / 1000000

    print('SF_POWER', SF_POWER/1000000,
          'COLD_PIPE_LOSSES', COLD_PIPE_LOSSES/1000000,
          'IDP_POWER', IDP_POWER/1000000,
          'HOT_PIPE_LOSSES', HOT_PIPE_LOSSES / 1000000)



print(datasource.dataframe)

flag_1 = datetime.now()
delta_t = flag_1 - flag_0
print("Total runtime: ", delta_t.total_seconds())

save_results(datasource.dataframe)


