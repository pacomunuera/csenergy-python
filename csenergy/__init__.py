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
import pandas as pd
from tkinter import * 


PI = 3.141592653589793

_DC_MODEL_PARAMS = {
    'RubenBarbero': set(['sigma', 'e_ext', 'h_ext', 'U_rec']),
    'NaumFraidenraich': set(['U_rec', 'U_exts', 'Cp', 'w']),
    'Patnode': set(['a0', 'a1', 'a2', 'a3', 'b0', 'b2']),
    'ASHRAE': set(['a', 'b', 'c', 'd', 'e', 'f']),
    'Montes': set(['a0', 'a1', 'a2', 'a3', 'b0', 'b1', 'b2']),
    'Price': set(['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6'])
    }

#TODO: 1 Preparación dataframe inicial. A partir de CSV como punto de inicio.
#TODO: 2 Preparación de datos del HCE

' Campos del dataframe inicial: '
'    | timestamp | Tamb | Wind | Winddir | DNI | Vector Solar | MassFlow | Tin | Tout |'

#Dic for HTC characteristics 
HCE_type = {"Model ": 'RTC70', "Longitude ": 4, "Din ":0.06, "Dout ":0.12, "e_ext ": 1,
     "h_ext ": 1, "U_rec ": 1, "Sigma ": 1}

HCE_status = {"Hydrogen ": 0, "Vacuum ": 0, "Broken ": False}

HCE_operation = {"Tin ": 293, "Clean ": 1}

class HCE:   
    
    def __init__(self, sca, hce_order) -> None:
        self.long = 1
        self.Dro = 1
        self.Dri = 1
        self.massflow = 1
        self.eext = 1
        self.hext = 1
        self.hint = 1
        self.Urec = 1
        self.sigma = 1
        self.tin = 1
        self.tout = 1
        self.sca = sca
        self.hce_order = hce_order       
        
    
    
    def __PRBarbero4grade__() -> float:
        return 0
    
    def __PRBarbero1grade__() -> float:
        return 0
    
    def __PRBarbero0grade__() -> float:
        return 0
   
    def calcTempOut() -> float:
        '''
        HTF output temperature [ºC] Ec. 3.24 Barbero

        Returns
        -------
        float [ºC]
            DESCRIPTION.

        '''
        return 0      
    
    def q_abs() -> float:
        '''
        Ec. 3.20 Barbero

        Returns
        -------
        float [W]
            DESCRIPTION.
            Thermal power absorbed by the HCE
        '''
        return 0       
 




   
class SCA(object):
    
    SCA_configuration ={"HCE Number ": 12, "Position in Loop": 1,
                    "Temperature probes number": 1,
                    "Temp. probes position": 'Middle', "Defocus Order": 1}
    
    #Uso de __slots__ permite ahorrar memoria RAM
#    __slots__= ['SCA_configuration']  
    def __init__(self, loop, sca_order):
        
#        self.SCA_type = {    
#                'HCEperSCA': 24,
#                'Loop': 1,
#                'positionInLoop': 1,
#                'operationStatus': 0,
#                'SCAtemperature': 0
#                }
#        self.SCA_tin = SCA.firstHCEinSCA(self).HCE.tin
        self.loop = loop
        self.sca_order = sca_order
        self.hces = []
              
    def get_tin(self, sca):
         
        if not sca:
            sca = self
        if sca.sca_order > 0:
            return SCA.get_tout(sca.loop.scas[sca.sca_order-1])
        else:
            return Loop.get_tin(sca.loop)
     
    def get_tout(self, sca):
        if not sca:
            sca = self
            
        return HCE.get_tout(sca.hces[-1])
    
#    def __gettitem__(self, item):
#        return self.SCA_type[item]
#    
#    def setOperationStatus():
#        
#        if SCAtemperature > SCATempMax():
#            operationStatus = 'defocused'
#        else:
#            soperationStatus = 'tracking'
            
            
class Loop(object):
    def __init__(self, solarfield, loop_order):
        self.solarfield_loop = solarfield
        self.scas = []
        self._SCAperLOOP = 4
        
    def tempOut(self, weatherstatus, plantstatus):
        '''
        Calculation of the output temperature Tout of the HTF
        according to the input temperature, mass flow, wheather
        conditions, plant performance

        Returns
        -------
        Output Temperature, Tout

        '''
        pass
    
    def massFlow():
        '''
        Calculation of the current mass flow of HTF in the loop
        by dividing the total mass flow by the number of loops
        in the plant 
        

        Returns
        -------
        current HTF massflow, massFlowLoop

        '''

        pass
    
    def requiredMassFlow():
        '''
        
        Calculation of the massFlow necessary to obtain the output temperature
        required by the operator

        Returns
        -------
        required HTF mass flow, reqMassFlowLoop

        '''
        pass

class PTSolarField(object):
    '''
    Parabolic Trough Solar Field
    
    '''
    
    
    def __init__(self, plant):
        
        self.solarfield_plant = plant
        self.loops = []
        self.configuration = {'N_loops_per_solar_field': 30,
                              'N_SCA_per_loop': 4,
                              'N_HCE_per_SCA': 24
                              }
        self.massflow =1
        self.tin= 1
        self.tout = self.tin
    
    def get_configuration():
        return self.configuration
    
    
    def get_massflow(self):
        '''
        Calculates total massflow throughout the solar field as 
        a sum of each massflow in a loop belonging the solar field
        
        '''
        
    def get_tout(self):
        '''
        Calculates HTF output temperature throughout the solar field as a 
        weighted average based on the enthalpy of the mass flow in each 
        loop belonging the solar field
        
        '''       
        
        
        
class PTPlant(object):
    '''
    Parabolic Trough Concentrated Solar Power Plant
    
    '''
    
    def __init__(self, site):

        self.site = site
        self.solarfields = []
        pass
    
    
    
        
class Site(object):
    def __init__(self, lat=39, long=-3):
        
        self._lat = 39
        self._long = 3
        
class HTF(object):
    
    def __init__(self, *argv):
        self._name = "Therminol VP-1"

class Weather(object):   
    def __init__(self):
        pass
        
    def loadWheatherData(path):
        pass

class ScatterMask(object):
    
    _PARAMETERS = {'sigma': 1, 'clean_mirror': 1}
    
    def __init__(self, plant):
        parameterMatrix = plant.configuration()
        pass
    
    def define_scattered_parameters_model():
        pass
        return scatterd_parameters_model
    
    
    def fill_mask(self, scattered_parameters_model):
        
        pass
    
    
    pass


# format = '%Y-%m-%d %H:%M:%S'
# df['Datetime'] = pd.to_datetime(df['date'] + ' ' + df['time'], format=format)
# df = df.set_index(pd.DatetimeIndex(df['Datetime']))
        
    

class Simulation(object):
    ''' 
    Definimos la clase simulacion para representar las diferentes
    pruebas que lancemos, variando el archivo TMY, la configuracion del
    site, la planta, el modo de operacion o el modelo empleado. 
    
    '''
    
    def __init__(self, weather, site, plant, mask, model,df):
        self.t_initial = 1
        self.t_end = 1
        
    def run():
        pass
        
    def result():
        pass
        return df
    
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

#print("Definir el  número de subcampos", end='\n')
#input_solarfields = int(input())
#print("Definir el número de lazos en cada subcampo", end='\n')
#input_loops = int(input())
#print("Definir el número de SCA en cada lazo", end='\n')
#input_sca = int(input())
#print("Definir el número de HCE en cada SCA", end='\n')
#input_HCE = int(input())


wd = Tk()
wd.title("Data input")
#wd.resizable( True, False)
#wd.iconbitmap("appicon.ico")
wd.geometry()
#wd.config(bg="blue")
fm = Frame()

fm.pack(side="right", anchor="s", fill="both", expand="true")
fm.config(bg="grey")
fm.config(width="650", height="350")
#tb = Entry(miframe).place(50,200)
#tb.pack()
#lb = Label(fm, text="Datos Meteorológicos", fg="grey", font=("Comic Sans MS", 12)).place(x=100, y=200)
lb = Label(fm, text="Datos").place(x=100, y=200)
tb = Entry(fm)
tb.place(x=100,y=100)

wd.mainloop()
    
#wd.destroy



site = Site(39, -3)
plant = PTPlant(site)

input_solarfields = 4
input_loops = 30
input_sca = 4
input_hce = 12
input_plant_model = 'unscattered'
input_simulation_type = 'comparative' # Compared to real data
input_tout_setpoint = 391 #ºC
input_total_mass_flow = 240 # kg/s 


for sf in range(input_solarfields):
    plant.solarfields.append(PTSolarField(plant))
    for l in range(input_loops):
        plant.solarfields[sf].loops.append(Loop(plant.solarfields[sf], l))
        for s in range(input_sca):
            plant.solarfields[sf].loops[l].scas.append(SCA(plant.solarfields[sf].loops[l], s))
            for h in range(input_hce):
                 plant.solarfields[sf].loops[l].scas[s].hces.append(HCE(plant.solarfields[sf].loops[l].scas[s],h))

   
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
        

