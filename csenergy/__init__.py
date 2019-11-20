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

class HCE:
    
    def __init__(self):
        self.long = 1
        self.diameter = 1
        self.massflow = 1
        self.tin = 1
        self.tout = 1
        
        
    def __PRBarbero4grade__(self):
        return 0
    
    def __PRBarbero1grade__(self):
        return 0
    
h1 = HCE()

print(h1.__PRBarbero1grade__())
    
class SCA:
    def __init__(self, HCEperSCA=24):
        self._HCEperSCA = 24
        
        
class Loop:
    def __init__(self, SCAperLOOP = 4):
        self._SCAperLOOP = 4
        
class Plant:
    def __init__(self, LOOPperPLANT=120):
        self._LOOPperPLANT = 120
        
class Site:
    def __init__(self, lat=39, long=-3):
        self._lat = 39
        self._long = 3
        
class HTF:
    def __init__(self, *argv):
        self._name = "Therminol VP-1"
        
        

