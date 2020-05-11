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


tags = ['1310-FX-1800','1310-FX-1801','1310-FX-1809','1310-FX-1810',
        '1310-TT-1800','1310-TT-1801','1310-TT-1809','1310-TT-1810',
        '1310-PT-1800','1310-PT-1801','1310-PT-1802','1310-PT-1803',
        '1310-PT-1807','1310-PT-1808','1310-PT-1809','1310-PT-1810']


tabla_df = {}
tabla_dfH = {}
df_destino = pd.DataFrame()

for t in tags:
    name = './data/'+t+'.txt'
    tabla_df[t] = pd.read_csv(
        name, sep=';',
        decimal=',',
        index_col='Date/Time')
    tabla_df[t].index = pd.to_datetime(tabla_df[t].index)
    print(t)
    tabla_df[t].rename(columns={'Value':t}, inplace=True)
    tabla_dfH[t] = tabla_df[t].resample('1H').mean()
    df_destino=pd.merge(df_destino,
                        tabla_dfH[t],
                        how='outer',
                        left_index=True,
                        right_index=True)
    print(df_destino)

df_destino.to_csv('./data/tabla.csv', sep=';', decimal=',')


