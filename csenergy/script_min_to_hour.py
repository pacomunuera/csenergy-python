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



xls = pd.ExcelFile('./data/AÑO METEOROLOGICO MINUTAL.xlsx')

df = xls.parse(skiprows=2, index_col='date')

df.index = pd.to_datetime(df.index)

dfH = df.resample('1H').mean()

dfH.to_csv('./data/annual data.csv', sep=';', decimal=',')

