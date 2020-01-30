# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 12:27:20 2019

@author: paco
"""

import pandas as pd
from datetime import datetime
import numpy as np
'''resample(self, rule, how=None, axis=0, fill_method=None, closed=None,
    label=None, convention='start', kind=None, loffset=None, limit=None,
    base=0, on=None, level=None)
'''

file1 = "datos_prueba.csv"
#date_rng = pd.date_range(start='1/1/2014',end='31/12/2014',freq='H')
datos = pd.read_csv(file1, sep=';', decimal=',', index_col=0)
datos.index = pd.to_datetime(datos.index)
datos = datos.apply(pd.to_numeric, errors='coerce')
print(datos)
robj = datos.resample('10T').mean()
print(robj)

for i in range(3):
    print(i)
