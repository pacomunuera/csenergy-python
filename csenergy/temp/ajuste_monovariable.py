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

rowT =[]

for temp in range(288,670,10):
    for pressure in range (1100000, 2000000, 100000):
        rowT.append([PropsSI('D','T',temp,'P', pressure,'Water'),
                    temp, pressure])

dt = pd.DataFrame(rowT,columns=['y', 'temp','pressure'])
print(dt)

fig = pyplot.figure(figsize=(10, 10))       # Ajustes del gráfico
ax = Axes3D(fig)

x1 = dt['temp']  # Datos eje X
x2 = dt['pressure']  # Datos eje Y
y = dt['y']  # Datos eje Z (Var. Respuesta)

ax.scatter(x1, x2, y, marker='.', c='r')
ax.set_xlabel('Temperature')        # Etiqueta del eje X
ax.set_ylabel('Pressure')       # Etiqueta del eje Y
ax.set_zlabel('Cp(T,P');

mod = smf.ols('y ~ x1 + x2', data=dt).fit()  # Ajusta el modelo usando el registro natural de uno de los regresores
print(mod.summary())

