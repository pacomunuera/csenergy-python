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



T = []
rowmu = []

df = pd.read_csv('VP1 - MU.csv', sep='\t',
                        decimal= ',')

rowmu = df['mu'].to_list()
T = df['T'].to_list()

sols = {}
for grado in range(8,9):
  z = np.polyfit(T, rowmu, grado, full=True)
  sols[grado] = z

# Pintar datos
plt.plot(T, rowmu, 'o')

# Pintar curvas de ajuste
xp = np.linspace(T[0], T[-1]+100, 10000)
for grado, sol in sols.items():
  coefs, error, *_ = sol
  p = np.poly1d(coefs)
  plt.plot(xp, p(xp), "-", label="Gr: %s. Error %.30f" % (grado, error) )
  lista = list(coefs)
  lista.reverse()
  print(lista)

plt.legend()


