# -*- coding: utf-8 -*-

"""
__init__.py: simulation launcher
@autor: pacomunuera
2020
"""

from tkinter import *
import json
from datetime import datetime
import csenergy as cs

with open("./saved_configurations/TEST_2016_DOWA.json") as simulation_file:
    SIMULATION_SETTINGS = json.load(simulation_file)

SIM = cs.SolarFieldSimulation(SIMULATION_SETTINGS)

FLAG_00 = datetime.now()
SIM.runSimulation()
FLAG_01 = datetime.now()
DELTA_01 = FLAG_01 - FLAG_00
print("Total runtime: ", DELTA_01.total_seconds())
