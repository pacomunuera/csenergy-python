import csenergy as cs
import pandas as pd
from tkinter import *
import json
import copy
import sys
from datetime import datetime

with open("./saved_configurations/simulation_S_F_UVAC.json") as simulation_file:
    simulation_settings = json.load(simulation_file)

sim = cs.SolarFieldSimulation(simulation_settings)


flag_00 = datetime.now()
sim.runSimulation()
flag_01 = datetime.now()
delta_01 = flag_01 - flag_00
print("Total runtime: ", delta_01.total_seconds())






