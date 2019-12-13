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

class HCE:
    
    def __init__(self):
        self.long = 1
        self.diameter = 1
        self.massflow = 1
        self.eext = 1
        self.hext = 1
        self.urec = 1
        self.sigma = 1
        self.tin = 1
        self.tout = 1
        
        
    def __PRBarbero4grade__(self):
        return 0
    
    def __PRBarbero1grade__(self):
        return 0
    
h1 = HCE()

print(h1.__PRBarbero1grade__())
    
class SCA(object):
    def __init__(self, HCEperSCA=24):
        self._HCEperSCA = 24
        
        
class Loop(object):
    def __init__(self, SCAperLOOP = 4):
        self._SCAperLOOP = 4
        
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
        self.file = 1
        self.data = pd.read_csv("filename.csv")
        
''' Definimos la clase simulacion para representar las diferentes
    pruebas que lancemos, variando el archivo TMY, la configuracion del
    site, la planta, el modo de operacion o el modelo empleado. '''

class Simulation(object):
    def __init__(self, weather, site, plant, operation, model):
        self.t_init = 1
        self.t_end = 1
        
    
#         

