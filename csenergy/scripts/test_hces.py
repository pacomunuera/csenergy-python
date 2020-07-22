# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 11:29:44 2020

@author: paco
"""

# -*- coding: utf-8 -*-
import csenergy as cs
import pandas as pd
from tkinter import *
import json
import copy
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import rc


model_settings = {'max_err_t': 0.1,
                  'max_err_tro': 0.1,
                  'max_err_pr': 0.01}

model = cs.ModelBarbero4thOrder(model_settings)

with open("./fluid_files/fluids_lib_rev01.json") as fluid_files:
    set_fluids = json.load(fluid_files)

htf_settings = set_fluids[0]

htf = cs.FluidTabular(htf_settings)

tin = 616.16
pin = 2000000
mf = 5

row = []
row.append({'none': None})
row.append({'DryBulb': 298, 'Wspd': 0})


with open("./hce_files/HCE_library.json") as hces_files:
    set_hces = json.load(hces_files)


for h in set_hces:
    hce = cs.HCE(0,0, h)
    hce.tin = tin
    hce.pin = pin
    hce.massflow = mf
    for qabs in range(15000, 45000, 5):
        model.calc_pr(hce, htf, row)
        print('hce.pr',hce.pr, 'tout', tout)



