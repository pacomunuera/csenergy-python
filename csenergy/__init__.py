import csenergy as cs
import pandas as pd
from tkinter import *
import interface
import json

#  TO-DO: PASO PREVIO -> INTERFAZ GRAFICA PARA CREAR EL simulation_file

with open("./saved_configurations/simulation.json") as simulation_file:
    simulation_settings = json.load(simulation_file)

site = cs.Site(simulation_settings)

solarplant = cs.SolarPlant(simulation_settings['solar_plant'],
                           simulation_settings['sca'],
                           simulation_settings['hce'],
                           simulation_settings['hce_model_settings'])

hotfluid = cs.HotFluid(simulation_settings['hot_fluid'])
coldfluid = cs.ColdFluid(simulation_settings['cold_fluid'])
weather = cs.Weather(simulation_settings['weather'])

#while  not hasattr(weather, 'weatherdata'):
#    weather.openWeatherDataFile()

site.name = weather.get_weather_data_site()['City']

hcemask = cs.HCEScatterMask(simulation_settings['solar_plant'],
                            simulation_settings['hce_scattered_params'])
scamask = cs.SCAScatterMask(simulation_settings['solar_plant'],
                            simulation_settings['sca_scattered_params'])

cycle = cs.ThermodynamicCycle(simulation_settings['cycle'])

# TO-DO: MÁSCARAS PARA INTRODUCIR DISPERSIÓN EN LOS PARÁMETROS.
#hcemask.applyMask(solarplant)
#scamask.applyMask(solarplant)

model = cs.ModelBarbero4grade(simulation_settings['hce_model_settings'])

#model.simulateSolarPlant(solarplant, site, weather, hotfluid)

simulation = cs.Simulation('sim1')
simulation.simulateSolarPlant(model, solarplant, site, weather, hotfluid)


# Crea aplicación
#root = Tk()
## Creamos una ventana
#a = interface.main(root, simulation)
#
##Lanza la aplicación en continuo
#root.mainloop()

'''
pandas.merge_asof(left, right, on=None, left_on=None, right_on=None, left_index=False, right_index=False, by=None, left_by=None, right_by=None, suffixes=('_x', '_y'), tolerance=None, allow_exact_matches=True, direction='backward')[source]¶
'''

' Campos del dataframe inicial: '
'    | timestamp | Tamb | Wind | Winddir | DNI | Vector Solar | MassFlow | Tin | Tout |'

#  HCE_status = {"Hydrogen ": 0, "Vacuum ": 0, "Broken ": False}



