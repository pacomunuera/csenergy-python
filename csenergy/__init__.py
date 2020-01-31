#from pvlib.version import __version__
#from pvlib import tools
#from pvlib import atmosphere
#from pvlib import clearsky
## from pvlib import forecast
#from pvlib import irradiance
#from pvlib import location
#from pvlib import solarposition
#from pvlib import iotools
#from pvlib import ivtools
#from pvlib import tracking
#from pvlib import spa
import csenergy as cs
import pandas as pd
from tkinter import *
import interface
import json

'''
pandas.merge_asof(left, right, on=None, left_on=None, right_on=None, left_index=False, right_index=False, by=None, left_by=None, right_by=None, suffixes=('_x', '_y'), tolerance=None, allow_exact_matches=True, direction='backward')[source]¶
'''
# _DC_MODEL_PARAMS = {
#     'Barbero': set(['sigma', 'eext', 'hext', 'urec']),
#     'NaumFraidenraich': set(['urec', 'uexts', 'cp', 'w']),
#     'Patnode': set(['a0', 'a1', 'a2', 'a3', 'b0', 'b2']),
#     'ASHRAE': set(['a', 'b', 'c', 'd', 'e', 'f']),
#     'Montes': set(['a0', 'a1', 'a2', 'a3', 'b0', 'b1', 'b2']),
#     'Price': set(['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6'])
#     }

' Campos del dataframe inicial: '
'    | timestamp | Tamb | Wind | Winddir | DNI | Vector Solar | MassFlow | Tin | Tout |'

#  HCE_status = {"Hydrogen ": 0, "Vacuum ": 0, "Broken ": False}

with open("./saved_configurations/simulation.json") as simulation_file:
    simulation_settings = json.load(simulation_file)

site = cs.Site(simulation_settings)

solarplant = cs.SolarPlant(simulation_settings)

weather = cs.Weather(simulation_settings)

while  not hasattr(weather, 'weatherdata'):
    weather.openWeatherDataFile()

site.name = weather.get_weather_data_site()['City']

hcemask = cs.HCEScatterMask(simulation_settings)
scamask = cs.SCAScatterMask(simulation_settings)

hot_fluid = cs.HotFluid(simulation_settings)
cold_fluid = cs.ColdFluid(simulation_settings)
cycle = cs.ThermodynamicCycle(simulation_settings)

hcemask.applyMask(solarplant)
scamask.applyMask(solarplant)


model = cs.ModelBarbero4grade(simulation_settings)

model.simulateSolarPlant(solarplant, site, weather, hot_fluid)

# Crea aplicación
#root = Tk()
## Creamos una ventana
#a = interface.main(root, simulation)
#
##Lanza la aplicación en continuo
#root.mainloop()

# for sf in range(n_solarfields):
#     plant.solarfields.append(cs.SolarField(plant, sf))
#     for l in range(n_loops):
#         plant.solarfields[sf].loops.append(
#                 cs.Loop(plant.solarfields[sf], l))
#         for s in range(n_sca):
#             plant.solarfields[sf].loops[l].scas.append(
#                     cs.SCA(plant.solarfields[sf].loops[l], s))
#             for h in range(n_hce):
#                 if hce_type == 'Barbero':
#                     plant.solarfields[sf].loops[l].scas[s].hces.append(
#                          cs.HCE_Barbero(parameters, plant.solarfields[sf].loops[l].scas[s],h))
#                 elif hce_type == 'NaumFraidenraich':
#                     pass

#if model in _DC_MODEL_PARAMS.keys():
#                # validate module parameters
#                missing_params = _DC_MODEL_PARAMS[model] - \
#                                 set(self.system.module_parameters.keys())
#                if missing_params:  # some parameters are not in module.keys()
#                    raise ValueError(model + ' selected for the DC model but '
#                                     'one or more required parameters are '
#                                     'missing : ' + str(missing_params))
#                if model == 'sapm':
#                    self._dc_model = self.sapm
#                elif model == 'desoto':
#                    self._dc_model = self.desoto
#                elif model == 'cec':
#                    self._dc_model = self.cec
#                elif model == 'pvsyst':
#                    self._dc_model = self.pvsyst
#                elif model == 'pvwatts':
#                    self._dc_model = self.pvwatts_dc
#            else:
#                raise ValueError(model + ' is not a valid DC power model')
