# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 08:26:43 2019

@author: fmunuera
"""

import numpy as np
import scipy as scy


def f0() -> float:
    return qabs/(urec*(tfe-text))

def f1() -> float:
    return ((4*sigma*eext*text**3)+hext)/urec 

def f2() -> float:
    return 6*text**2*(sigma*eext/urec)*(qabs/urec)

def f3() -> float:
    return 4*text*(sigma*eext/urec)*(qabs/urec)**2

def f4() -> float: 
    return (sigma*eext/urec)*(qabs/urec)

def e3_36(pr0,f0,f1,f2,f3,f4) -> float:

    return  (1-pr0-
             f1*(pr0+1/f0)+
             f2*(pr0+1/f0)**2+
             f3*(pr0+1/f0)**3+
             f4*(pr0+1/f0)**4
             )

def e3_36prime(pr0,f0,f1,f2,f3,f4):
    
    return (-1-f1+
            2*f2*(pr0+1/f0)+
            3*f3*(pr0+1/f0)**2+
            4*f4*(pr0+1/f0)**3
            )

def pr0() -> float:
    
    return scy.optimize.newton(e3_36, pr0, fprime=e3_36prime, args=(f0,f1,f2,f3,f4))

def g(z, f0, f1, f2, f3, f4) -> float:
    
    return -(1+1/f0)+(1+f1)*z+f2*z**2+f3*z**3+f4*z**4

def g1(z, f0, f1, f2, f3, f4) -> float:
    
    return f1+2*f2*z+3*f3*z**2+4*f4*z**3
        
def g2(z, f0, f1, f2, f3, f4) -> float:
    
    return 2*f2+6*f3*z+12*f4*z**2

def g3(z, f0, f1, f2, f3, f4) -> float:
    
    return 6*f3+24*f4*z

def pr(pr0, g1, g2, g3, NTU, x) -> float:
    
    return ((pr0*g1/(1-g1))*(1/(NTU*x))*(np.exp((1-g1)*NTU*x/g1)-1) -
            (g2/(6*g1))*(pr0*NTU*x)**2 -
            (g3/(24*g1)*(pr0*NTU*x)**3)
            )
    
            


class HCE(object): 
  
    def __init__(self, sca, hce_order):
        self.long = 1
        self.dro = 1
        self.dri = 1
        self.massflow = sca.loop.massFlow
        self.eext = 1
        self.hext = 1
        self.hint = 1
        self.urec = 1
        self.sigma = 1
        self.pr = 1
        self.tin = 1
        self.tout = self.tin
        self.sca = sca
        self.hce_order = hce_order
        
    def qabs(pr_opt, cg, dni, pr_shw, pr_geo) -> float:
        return pr_opt*cg*dni*pr_shw*pr_geo
    
    def qu(urec, tro, tf) -> float:
        return urec*(tro-tf) #Ec. 3.21
    
    def urec(hint, dro, dri, krec) -> float:
        return 1/((1/hint) + (dro*np.log(dro/dri))/(2*krec)) #Ec. 3.22
        
    def qu_from_pr(pr, qabs) -> float: 
        return pr*qabs #Ec. 3.26
    
    
    def __PRBarbero4grade__() -> float:
        return 0
    
    def __PRBarbero1grade__() -> float:
        return 0
    
    def __PRBarbero0grade__() -> float:
        return 0
   
    def calcTempOut(tfe, dro, x, qabs, pr, massflow, cp) -> float:
        '''
        HTF output temperature [ºC] Ec. 3.24 Barbero

        Returns
        -------
        float [ºC]
            DESCRIPTION.

        '''
        return tfe + PI*dro*x*qabs*pr/(massflow*cp)   
    
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
    def __init__(self, solarfield, loop_order, massFlow):
        self.solarfield_loop = solarfield
        self.scas = []
        self._SCAperLOOP = 4
        self.massFlow = massFlow
        
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
    
    def calcRequiredMassFlow(self):
        '''
        
        Calculation of the massFlow necessary to obtain the output temperature
        required by the operator

        Returns
        -------
        required HTF mass flow, reqMassFlowLoop

        '''
        self.requiredMassFlow = 1
        


class SolarField(object):
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
    def calcRequiredMassFlow(self):
        req_massflow = 0
        for l in self.loops:
            req_massflow += l.requiredMassFlow
        self.requieredMassFlow = req_massflow/self.loops.len()
        
        
        
class Plant(object):
    '''
    Parabolic Trough Concentrated Solar Power Plant
    
    '''
    
    def __init__(self, site):

        self.site = site
        self.solarfields = []
        self.solarfield_to_exchanger_thermal_lost = 0.1
        self.exchanger_performance = 0.9
        self.exchanger_to_turbogroup_thermal_lost = 0.1
        self.steam_cycle_performance = 0.9
        self.turbogenerator_performance = 0.9
 
        
        pass

    def calcRequiredMassFlow(self):
        req_massflow = 0
        for sf in self.solarfields:
            req_massflow += sf.calcRequiredMassFlow()
        self.reqMassFlow = req_massflow
    
        
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
    
    def __init__(self, site, plant, mask, df, simulation_type):
        self.site = site
        self.plant = plant
        self.mask = mask
        self.data = df
        self.type = simulation_type
        
        
    def calcRequiredMassFlow():
        '''calcula el caudal promedio requerido en cada lazo para alzanzar
        la temperatura deseada. Calcula el caudal en cada lazo (si son todos 
        iguales solo lo hace una vez) y después calcula el promedio en cada 
        subcampo. De esta forma tenemos el caudal '''
        
        pass
        
    def calc():
        pass
        return df