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
import CoolProp.CoolProp as CP
from datetime import datetime
import os.path

import pandas as pd                    ## Este proporciona una estructura similiar a los data.frame
import statsmodels.api as sm           ## Este proporciona funciones para la estimación de muchos modelos estadísticos
import statsmodels.formula.api as smf

from sklearn import datasets, linear_model
from matplotlib import pyplot             # Permite la generación de gráficos
from mpl_toolkits.mplot3d import Axes3D   # Permite agregar eje tridimensionales
import random
'''
A call to the top-level function PropsSI can provide: temperature, pressure,
density, heat capacity, internal energy, enthalpy, entropy, viscosity and
thermal conductivity. Hence, the available output keys are:
    T, P, D, C, U, H, S, V, L, Tmin and Tmax.
'''
#PropsSI('T','P',101325,'Q',0,'Water')

def get_deltaH(t, p, coolpropID):
    
    CP.set_reference_state(coolpropID,'ASHRAE')
    
    deltaH = PropsSI('H','T',t ,'P', p, coolpropID)
    
    CP.set_reference_state(coolpropID, 'DEF')

    return deltaH

def get_T(H, p, coolpropID):
    
    CP.set_reference_state(coolpropID,'ASHRAE')
    
    temperature = PropsSI('T', 'H', H, 'P', p, coolpropID)
    
    CP.set_reference_state(coolpropID, 'DEF')
    
    return temperature


rowT =[]
pressure = 20000000
temp = 623.15
for temp in range(288,670,10):
    for pressure in range(2000000,20000000,1000000):
        rowT.append([get_deltaH(temp, pressure, 'INCOMP::TVP1'),
                     temp, pressure])

dt = pd.DataFrame(rowT,columns=['y', 'temp','pressure'])
print(dt)

fig = pyplot.figure(figsize=(10, 10))       # Ajustes del gráfico
ax = Axes3D(fig)

x1 = dt['temp']  # Datos eje X
x2 = dt['pressure']
y = dt['y']  # Datos eje Z (Var. Respuesta)

ax.scatter(x1, x2, y, marker='.', c='r')
ax.set_xlabel('Temperature')        # Etiqueta del eje X
ax.set_ylabel('Pressure')       # Etiqueta del eje Y
ax.set_zlabel('H(T,P');

mod = smf.ols('y ~ x1 + x1^2 + x1^3 + x1^4 + x1^5', data=dt).fit()  # Ajusta el modelo usando el registro natural de uno de los regresores
print(mod.summary())

