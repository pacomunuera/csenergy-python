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

#TODO: 1 Preparación dataframe inicial. A partir de CSV como punto de inicio.
#TODO: 2 Preparación de datos del HCE

' Campos del dataframe inicial: '
'    | timestamp | Tamb | Wind | Winddir | DNI | Vector Solar | MassFlow | Tin | Tout |'

class HCE:
    
    def __init__(self):
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
        self.SCA = 1
        self.SCAposition = 1
        
        
    def __PRBarbero4grade__():
        return 0
    
    def __PRBarbero1grade__():
        return 0
    
    def __PRBarbero0grade__():
        return 0
   
    def calcTempOut():
        '''
        HTF output temperature [ºC] Ec. 3.24 Barbero

        Returns
        -------
        float [ºC]
            DESCRIPTION.

        '''
        return 0
    
    def qabs():
        '''
        Ec. 3.20 Barbero

        Returns
        -------
        float [W]
            DESCRIPTION.
            Thermal power absorbed by the HCE
        '''
        return 0       
    
    
h1 = HCE()

print(h1.__PRBarbero1grade__())
    
class SCA(object):
    
    self.HCEperSCA = 24
    self.HCEArray = ()
    self.Loop = 1
    self.positionInLoop = 1
    self.operationStatus = 0
    self.SCAtemperature = 0
    
    def __init__(self, pIL, HCEperSCA=24):
        self.positionInLoop = pIL
    
    def setOperationStatus():
        
        if self.SCAtemperature > self.SCATempMax():
            self.operationStatus = 'defocused'
        else:
            self.operationStatus = 'tracking'
            
            
class Loop(object):
    def __init__(self, SCAperLOOP = 4):
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
    
class Plant(object):
    def __init__(self, LOOPperPLANT=120):
        self._LOOPperPLANT = 120
        
class Site(object):
    def __init__(self, lat=39, long=-3):
        
        self._lat = 39
        self._long = 3
        
class HTF(object):
    
    def __init__(self, *argv):
        self._name = "Therminol VP-1"

class Weather(object):   
    def __init__(self, file):
        pass
        
    def loadWheatherData(path):
        self.wheatheData = read_tmy2(path)

# format = '%Y-%m-%d %H:%M:%S'
# df['Datetime'] = pd.to_datetime(df['date'] + ' ' + df['time'], format=format)
# df = df.set_index(pd.DatetimeIndex(df['Datetime']))
        
    

class Simulation(object):
    ''' 
    Definimos la clase simulacion para representar las diferentes
    pruebas que lancemos, variando el archivo TMY, la configuracion del
    site, la planta, el modo de operacion o el modelo empleado. 
    
    '''
    def __init__(self, weather, site, plant, operation, model):
        self.t_init = 1
        self.t_end = 1
        
    def result():
        pass
        return df
        
    
weather_file = Weather()
fluid = HTF()
site = Site()
plant = Plant()

simulation = Simulation(self, weather_file, plant, operation, model)
simultation_output = simulation.result()


        

