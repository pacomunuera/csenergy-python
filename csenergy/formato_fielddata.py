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


df = pd.read_csv('./fielddata_files/prueba.csv', skiprows=0, sep=';', decimal=b',', index_col=0)
df.index = pd.to_datetime(df.index)
df = df.resample('1H').mean()

print(df.head)

df.to_csv('./fielddata_files/1Hdata.csv', index=True,
          header=True, decimal=',', sep=';', float_format='%.3f')


