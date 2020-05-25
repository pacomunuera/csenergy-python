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
import json
from json import encoder
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
with open("fluid_files/fluids_lib.json") as file:
    fluids = json.load(file)

    fluid_lib = []

    for fs in fluids:

        f = {}

        f['name'] = fs['name']
        f['tmin'] = fs['tmin']
        f['tmax'] = fs['tmax']
        # f['cp'] = []
        # f['rho'] = []
        # f['mu'] = []
        # f['kt'] = []
        # f['h'] = []
        # f['t'] = []
        # for c in fs['coefficients']:
        #     .extend([*coeff.values()])

        for parameter in fs['coefficients']:
            coeff_list = list(parameter.keys())
            for coeff in parameter.values():
                coeff_list.extend([*coeff.values()])
            f[coeff_list[0]] = coeff_list[1:]
            print(coeff_list[1:])

        fluid_lib.append(f)

#print(fluid_lib)



