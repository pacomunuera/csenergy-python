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

t1877 = ['1300-AXT-1877R','1300-BXT-1877R','1300-CXT-1877R']
t1873 = ['1300-AXT-1873R','1300-BXT-1873R','1300-CXT-1873R']
t1874 = ['1300-AXT-1874R','1300-BXT-1874R','1300-CXT-1874R']
t1871 = ['1300-AXT-1871R','1300-BXT-1871R','1300-CXT-1871R']
t1879 = ['1300-AXT-1879R','1300-BXT-1879R','1300-CXT-1879R']
t1880 = ['1300-AXT-1880R','1300-BXT-1880R','1300-CXT-1880R']



df = pd.read_csv('./data/annual data.csv', index_col='date', sep=';',
                  decimal=',')

df.index = pd.to_datetime(df.index)

# DataFrame.median(self, axis=None, skipna=None, level=None, numeric_only=None, **kwargs)[source]

df['DNI'] = df[t1877].median(axis=1)
df['DryBulb'] = df[t1873].median(axis=1)
df['Hum'] = df[t1874].median(axis=1)
df['Wspd'] = df[t1871].median(axis=1)
df['Elev'] = df[t1879].median(axis=1)
df['Azimuth'] = df[t1880].median(axis=1)

columns_to_drop = t1877 + t1873 + t1874 + t1871 + t1879 + t1880

df.drop(columns = columns_to_drop, inplace = True)

df.to_csv('./data/annual data 2.csv', sep=';', decimal=',')

