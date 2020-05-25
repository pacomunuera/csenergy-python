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


with open("./saved_configurations/loop.json") as simulation_file:
    simulation_settings = json.load(simulation_file)

simulation = cs.LoopSimulation(simulation_settings)
simulation.runSimulation()



