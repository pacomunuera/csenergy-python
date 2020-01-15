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




# Crear interface

# La interfaz nos ayuda a crear el objeto simulation

# 




'''

pandas.merge_asof(left, right, on=None, left_on=None, right_on=None, left_index=False, right_index=False, by=None, left_by=None, right_by=None, suffixes=('_x', '_y'), tolerance=None, allow_exact_matches=True, direction='backward')[source]¶
'''

# PI = 3.141592653589793

# _DC_MODEL_PARAMS = {
#     'Barbero': set(['sigma', 'eext', 'hext', 'urec']),
#     'NaumFraidenraich': set(['urec', 'uexts', 'cp', 'w']),
#     'Patnode': set(['a0', 'a1', 'a2', 'a3', 'b0', 'b2']),
#     'ASHRAE': set(['a', 'b', 'c', 'd', 'e', 'f']),
#     'Montes': set(['a0', 'a1', 'a2', 'a3', 'b0', 'b1', 'b2']),
#     'Price': set(['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6'])
#     }

#TODO: 1 Preparación dataframe inicial. A partir de CSV como punto de inicio.
#TODO: 2 Preparación de datos del HCE

' Campos del dataframe inicial: '
'    | timestamp | Tamb | Wind | Winddir | DNI | Vector Solar | MassFlow | Tin | Tout |'

#Dic for HTC characteristics 
# HCE_type = {"Model ": 'RTC70', "Longitude ": 4, "Din ":0.06, "Dout ":0.12, "e_ext ": 1,
#      "h_ext ": 1, "U_rec ": 1, "Sigma ": 1}

# HCE_status = {"Hydrogen ": 0, "Vacuum ": 0, "Broken ": False}

# HCE_operation = {"Tin ": 293, "Clean ": 1}


    
"""
El usuario selecciona:
    
    - Fichero Weather, se crea dataframe con registros a intervalos temporales
    - Configuración de planta (N_campos, N_lazos, N_SCA_lazo, N_HCE_SCA)
    - Modelo de planta:
        - Modelo 1: Parámetros sin dispersión
        - Modelo 2: Parámetros con dispersión
    - Tipo de simulacion:
        - Simulacion 1: Comparación con datos de generacion
        - Simulacion 2: Simulación de planta
    
    - Se genera la planta:
        - Modelo 1: Por ser todos iguales no sería necesario crear N_lazosxN_SCA_lazoxN_HCE_SCA
        sino que se va a recurrir a bucles para generar las salidas de temperatura y guardarlos
        en el dataframe
        
        - Modelo 2: Se crean los N_lazosxN_SCA_lazoxN_HCE_SCA HCE con sus atributos cargados según 
        las funciones de dispersión que toque. 
        
    - Simulación:
        - Simulación 1: para cada registro temporal se toma el caudal real y se calcula la temperatura de
        salida con el modelo seleccionado. Hay que manejar las temperaturas por encima del setpoint como pérdidas de energía por 
        culpa del desenfoque
        
        - Simulación 2:  para cada registro temporal se calcula el caudal para conseguir la temperatura
        de salida.
        
        datos que debe almacenarse en cada HCE
        - Rendimiento, tin, tout, massflow
        
        datos que debe almacenarse en cada SCA
        - Rendimiento (promedio de sus HCE), tin, tout, t_probe,  massflow, status (desenfoque, etc...)
        
        datos que debe almacenarse en cada LOOP
        - Rendimiento (promedio de sus SCA), tin, tout, massflow
"""        
plant_model = 'unscattered'
simulation_type = 'comparative' # Compared to real data
tin = 291 
tout_setpoint = 391 #ºC
total_mass_flow = 240 # kg/s 
#massFlow = total_mass_flow / n_loops


#weather = cs.Weather()
#simulation = cs.Simulation()

parameters = {'eext': 1,
              'hext': 1,
              'hint': 1,
              'urec': 1,
              'sigma' : 1,
              'cg' : 1,
              'pr_shw': 1,
              'pr_geo': 1,
              'pr_opt': 1
              }

hce_type = 'Barbero'

with open ("./saved_configurations/simulation.json") as simulation_file:
    simulation_settings = json.load(simulation_file)

site = cs.Site(simulation_settings.get('site'))
    
plant = cs.Plant(site)

mask = cs.ScatterMask(simulation_settings)

for sf in simulation_settings.get('plant').get('solarfields'):
    plant.solarfields.append(cs.SolarField(plant, sf))
    for l in range(sf.get('loops')):        
        plant.solarfields[-1].loops.append(
                cs.Loop(plant.solarfields[-1], l))
        for s in range(sf.get('scas')):
            plant.solarfields[-1].loops[-1].scas.append(
                    cs.SCA(plant.solarfields[-1].loops[-1], s))
            for h in range (sf.get('hces')):
                    plant.solarfields[-1].loops[-1].scas[-1].hces.append(
                        cs.HCE(plant.solarfields[-1].loops[-1].scas[-1],
                               h, 
                               simulation_settings.get('model_settings')))
                                      
cs.Model.applyMask(mask, plant)
   
#for sf in plant.solarfields:
#    print(type(sf),'\n')
#    for l in sf.loops:
#        print(type(l),'\n')
#        for s in l.scas:
#            print(type(s),'\n')
#            for h in s.hces:
#                print(h,"-->", h.parameters['eext'])


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





#mask = cs.ScatterMask(plant, simulation)
               
#for sf in plant.solarfields:
#    for l in sf.loops:
#        for s in l.scas:
#            for h in s.hces:
#                h.applyMaskToHCE(mask)
#                
#
#for sf in plant.solarfields:
#    for l in sf.loops:        
#        l.set_massflow()
#        l.set_tin()

                
                

               

#simulation = cs.Simulation(site, plant, mask, df, simulation_type)
              
#se aplicaca la máscara si procede            
   
#weather_file = Weather()
#fluid = HTF()
#site = Site()

#sf = solarField()

#mask = ScatterMask()

#plant = Plant()
#
#simulation = Simulation(self, weather_file, plant, operation, model)
#simultation_output = simulation.result()

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
        

