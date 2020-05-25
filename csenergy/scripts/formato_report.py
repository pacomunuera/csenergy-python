# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 12:06:36 2020

@author: paco
"""

import numpy as np
import scipy as sc
from scipy import constants
import math as mt
from CoolProp.CoolProp import PropsSI
from datetime import datetime
import os.path

import pandas as pd                    ## Este proporciona una estructura similiar a los data.frame
import statsmodels.api as sm           ## Este proporciona funciones para la estimación de muchos modelos estadísticos
import statsmodels.formula.api as smf

from sklearn import datasets, linear_model
from matplotlib import pyplot             # Permite la generación de gráficos
from mpl_toolkits.mplot3d import Axes3D   # Permite agregar eje tridimensionales
import random


#df = pd.read_fwf('../fielddata_files/ReportM4-1.txt')

#df = pd.read_csv('../fielddata_files/ReportM4-1.txt', skiprows=0, sep=' ', decimal=b',', index_col=0)
#df.index = pd.to_datetime(df.index)
#df = df.resample('1H').mean()


tags = []
datos = []
grupoA = []
grupoB =[]
points = 0

with open('../fielddata_files/ReportM4-1.txt', 'r') as f:
    for l in f:
        if 'Point Name:' in l:
            points += 1
            tag = l.split()
            for key in tag:
                if '@' in key:
                    tags.append(key)
            if len(grupoA)>0:
                datos.append([grupoA.copy(),grupoB.copy()])
            del grupoA[:]
            del grupoB[:]
        else:
            if (not l.isspace() > 0 and len(tags) > 0):
                lst = l.split()
                del lst[3:]
                grupoA.append(lst[0]+" "+lst[1])
                grupoB.append(lst[2])

    datos.append([grupoA.copy(),grupoB.copy()])

print("Points:", points)
print("Tags:", len(tags))
print("Datos:", len(datos))
#print("GrupoA:", datos[70][0])

dic = dict(zip(tags, datos))

date_rng = pd.date_range(start='17/06/2019',end='18/06/2019',freq='H')
base = {'date': date_rng}
df_H = pd.DataFrame(base, index=base['date'])
df_H.set_index('date', inplace=True)
claves = list(dic.keys())
cont = 0
for k in claves:
    #if k == '1100-AXT-1871R.UNIT1@ASTE1A':
    cont += 1
    d = {'date': dic[k][0], k: dic[k][1]}
    df = pd.DataFrame(d,
                      index=d['date'], dtype=float)
    df.index = pd.to_datetime(df.index)
    df = df.resample('1H').mean()
    
    #df.to_csv('../fielddata_files/report'+str(cont)+'.csv', header=True, decimal=',', sep=';', float_format='%.3f')
    
    #df_H.join(df, on='date')
    df_merged = df_H.merge(df, left_index=True, right_index=True)
    df_H = df_merged.copy()
    df.drop(index = df.index)
    print("paso", cont, " de ", points, " TAG:", k)

#df.index = pd.to_datetime(df.index)
#
#for i in range(2):
#    for j in range(2):
#        print("datos[{0}][{1}]={2}".format(i, j, datos[i][j]))

#df = pd.DataFrame(dict)
#print(df.head)
#df.to_csv('../fielddata_files/reportdf.csv', header=True, decimal=',', sep=';', float_format='%.3f')

df_H.to_csv('../fielddata_files/reportdf_H.csv', header=True, decimal=',', sep=';', float_format='%.3f')
