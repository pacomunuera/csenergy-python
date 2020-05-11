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
import os

import pandas as pd                    ## Este proporciona una estructura similiar a los data.frame
import statsmodels.api as sm           ## Este proporciona funciones para la estimación de muchos modelos estadísticos
import statsmodels.formula.api as smf

from sklearn import datasets, linear_model
from matplotlib import pyplot             # Permite la generación de gráficos
from mpl_toolkits.mplot3d import Axes3D   # Permite agregar eje tridimensionales
import random
import codecs

tags = ['1310-FX-1800','1310-FX-1801','1310-FX-1809','1310-FX-1810',
        '1310-TT-1800','1310-TT-1801','1310-TT-1809','1310-TT-1810',
        '1310-PT-1800','1310-PT-1801','1310-PT-1802','1310-PT-1803',
        '1310-PT-1807','1310-PT-1808','1310-PT-1809','1310-PT-1810']



files = os.listdir('./data/')
print(files)

f_destino = {}
for t in tags:
    name = './data/'+t+'.txt'
    f_destino[t] = open(name, "a+")
    f_destino[t].write('Date/Time; Point Name; Value')


for f in files:
    if 'PT.txt' in f:
        print("procesando "+ f)
        doc = codecs.open('./data/'+f,'rU','UTF-16') #open for reading with "universal" type set
        df = pd.read_csv(doc, sep='\t')
        df = df[['Date/Time','Point Name', 'Value']]
        for row in df.iterrows():
            if row[0] % 1000 == 0:
                print("row ", row[0])
            f_destino[row[1]['Point Name']].write(
                row[1]['Date/Time']+';'+
                row[1]['Value']+ os.linesep)

for t in tags:
    f_destino[t].close()

