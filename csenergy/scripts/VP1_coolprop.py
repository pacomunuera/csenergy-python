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
import numpy as np
import matplotlib.pyplot as plt

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

def get_density(t, p, coolpropID):

    CP.set_reference_state(coolpropID,'ASHRAE')
    rho = PropsSI('D','T',t,'P', p, coolpropID)
    CP.set_reference_state(coolpropID, 'DEF')
    return rho

def get_dynamic_viscosity(t, p, coolpropID):

    CP.set_reference_state(coolpropID,'ASHRAE')
    mu = PropsSI('V','T',t,'P', p, coolpropID)
    CP.set_reference_state(coolpropID, 'DEF')
    return mu

def get_cp(t, p, coolpropID):

    CP.set_reference_state(coolpropID,'ASHRAE')
    cp = PropsSI('C','T',t,'P', p, coolpropID)
    CP.set_reference_state(coolpropID, 'DEF')
    return cp

def get_kt(t, p, coolpropID):

    CP.set_reference_state(coolpropID,'ASHRAE')
    kt = PropsSI('L','T',t,'P', p, coolpropID)
    CP.set_reference_state(coolpropID, 'DEF')
    return kt

def get_T(H, p, coolpropID):

    CP.set_reference_state(coolpropID,'ASHRAE')
    temperature = PropsSI('T', 'H', H, 'P', p, coolpropID)
    CP.set_reference_state(coolpropID, 'DEF')
    return temperature

T = []
rowT =[]
pressure = 2000000
rowH =[]
rowkt = []
rowrho = []
rowcp = []
rowmu = []
rownu = []

for temp in range(288,670,2):
    T.append(temp)
    rowH.append(get_deltaH(temp, pressure, 'INCOMP::TVP1'))
    rowkt.append(get_kt(temp, pressure, 'INCOMP::TVP1'))
    rowcp.append(get_cp(temp, pressure, 'INCOMP::TVP1'))
    rowrho.append(get_density(temp, pressure, 'INCOMP::TVP1'))
    rowmu.append(get_dynamic_viscosity(temp, pressure, 'INCOMP::TVP1'))
    rownu.append(get_dynamic_viscosity(temp, pressure, 'INCOMP::TVP1') /
                 get_density(temp, pressure, 'INCOMP::TVP1'))
for h in rowH:
    rowT.append(get_T(h, pressure, 'INCOMP::TVP1'))

datos = {'T': T, 'cp': rowcp,'rho': rowrho, 'mu': rowmu,
         'nu': rownu, 'kt': rowkt, 'Hcalc': rowH, 'Tcalc': rowT}

sols = {}
for grado in range(10,12):
  z = np.polyfit(T, rowmu, grado, full=True)
  sols[grado] = z

# Pintar datos
plt.plot(T, rowmu, 'o')

# Pintar curvas de ajuste
xp = np.linspace(rowT[0], rowT[-1]+100, 10000)
for grado, sol in sols.items():
  coefs, error, *_ = sol
  p = np.poly1d(coefs)
  plt.plot(xp, p(xp), "-", label="Gr: %s. Error %.30f" % (grado, error) )
  lista = list(coefs)
  lista.reverse()
  print(lista)

plt.legend()

dt = pd.DataFrame(datos,
                  columns=['T', 'cp','rho', 'mu', 'nu', 'kt', 'Hcalc', 'Tcalc'])

dt.to_csv('t.csv',decimal=',', index=False, sep=';')
