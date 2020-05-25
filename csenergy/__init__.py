# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 11:29:44 2020

@author: paco
"""

from tkinter import *
import json
import copy
import sys
from datetime import datetime
import pandas as pd
import csenergy as cs

with open("./saved_configurations/TEST_SAM_DOW_A.json") as simulation_file:
    SIMULATION_SETTINGS = json.load(simulation_file)

SIM = cs.SolarFieldSimulation(SIMULATION_SETTINGS)

FLAG_00 = datetime.now()
SIM.runSimulation()
FLAG_01 = datetime.now()
DELTA_01 = FLAG_01 - FLAG_00
print("Total runtime: ", DELTA_01.total_seconds())
