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
import CoolProp
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
http://www.coolprop.org/fluid_properties/HumidAir.html

Input/Output parameters
Parameter	Units	Input/Output	Description
B, Twb, T_wb, WetBulb	K	Input/Output	Wet-Bulb Temperature
C, cp	J/kg dry air/K	Output	Mixture specific heat per unit dry air
Cha, cp_ha	J/kg humid air/K	Output	Mixture specific heat per unit humid air
D, Tdp, DewPoint, T_dp	K	Input/Output	Dew-Point Temperature
H, Hda, Enthalpy	J/kg dry air	Input/Output	Mixture enthalpy per dry air
Hha	J/kg humid air	Input/Output	Mixture enthalpy per humid air
K, k, Conductivity	W/m/K	Output	Mixture thermal conductivity
M, Visc, mu	Pa-s	Output	Mixture viscosity
psi_w, Y	mol water/mol humid air	Input/Output	Water mole fraction
P	Pa	Input	Pressure
P_w	Pa	Input	Partial pressure of water vapor
R, RH, RelHum	 	Input/Output	Relative humidity in [0, 1]
S, Sda, Entropy	J/kg dry air/K	Input/Output	Mixture entropy per unit dry air
Sha	J/kg humid air/K	Input/Output	Mixture entropy per unit humid air
T, Tdb, T_db	K	Input/Output	Dry-Bulb Temperature
V, Vda	m 3 /kg dry air	Input/Output	Mixture volume per unit dry air
Vha	m 3 /kg humid air	Input/Output	Mixture volume per unit humid air
W, Omega, HumRat	kg water/kg dry air	Input/Output	Humidity Ratio
Z	 	Output	Compressibility factor (Z=pv/(RT))
'''

#print(CP.HAProps('D', 'T', 300, 'P', 101.325, 'R', 0.5))

rowT =[]
for temp in range(288, 363 , 5):
    for press in range(90, 102, 1):
        rowT.append([CP.HAProps('D', 'T', temp, 'P', press, 'R', 0.5),temp, press])

dt = pd.DataFrame(rowT,columns=['y', 'temp','pressure'])
print(dt)

fig = pyplot.figure(figsize=(10, 10))       # Ajustes del gráfico
ax = Axes3D(fig)

x1 = dt['temp']  # Datos eje X
x2 = dt['pressure']
y = dt['y']  # Datos eje Z (Var. Respuesta)

ax.scatter(x1, y, marker='.', c='r')
ax.set_xlabel('Temperature')        # Etiqueta del eje X
ax.set_ylabel('V(T)');

mod = smf.ols('y ~ x1 + x1^2 + x2^1 + x2^2', data=dt).fit()
# Ajusta el modelo usando el registro natural de uno de los regresores
print(mod.summary())

