# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 08:26:43 2019

@author: pacomunuera
"""

import numpy as np
import scipy as scy
import math as mt
import CoolProp as CP


class HTF(object):
    
    def __init__(self, name):
        pass
    
    
class Model(object):
    
    def __init__(self, hce, model_settings, mask = None):
        self.name = model_settings.get('model_settings').get('name')
    
    @classmethod
    def set_tin(cls, hce):
        if hce.hce_order > 0:
            hce.tin = hce.sca.hces[hce.hce_order-1].tout
        elif hce.sca.sca_order > 0:
            hce.tin = hce.sca.loop.scas[hce.sca.sca_order-1].hces[-1]
        else:
            hce.tin = hce.sca.loop.tin
    
    @classmethod
    def applyMask(cls, mask, plant):
        
        #for sf in simulation.get('plant').get('solarfields')
        for sf in plant.solarfields:
            for l in sf.loops:
                for s in l.scas:
                    for h in s.hces:
                        for k in mask.matrix[sf.name][l.loop_order][s.sca_order][h.hce_order].keys():
                            h.parameters[k] = h.parameters[k] * float(mask.matrix[sf.name][l.loop_order][s.sca_order][h.hce_order][k]) 
                            print(mask.matrix[sf.name][l.loop_order][s.sca_order][h.hce_order][k])
                            

            
    
        
class modelBarbero4(Model):

    def __ini__(self, hce, model_settings, mask = None):
        
        super(Model, modelBarbero4).__init__(model_settings)
    

    
    
        if mask == None:
            hce.eext = model_settings.get('eext')
            hce.hext = model_settings.get('hext')
            hce.hint = model_settings.get('hext')
            hce.urec = model_settings.get('urec')
            hce.sigma = model_settings.get('sigma')
            hce.cg = model_settings.get('cg')
            hce.pr_shw = model_settings.get('pr_shw')
            hce.pr_geo = model_settings.get('pr_geo')
            hce.pr_opt = model_settings.get('pr_opt')
        
        

    def calcTempOut(tfe, dro, x, qabs, pr, massflow, cp) -> float:
        '''
        HTF output temperature [ºC] Ec. 3.24 Barbero

        Returns
        -------
        float [ºC]
            DESCRIPTION.

        '''
        return tfe + mt.PI*dro*x*qabs*pr/(massflow*cp)   

    def set_qabs(self, pr_opt, dni, pr_shw, pr_geo):
        '''
        Ec. 3.20 Barbero

        Returns
        -------
        float [W]
            DESCRIPTION.
            Thermal power absorbed by the HCE
        '''
        self.qabs = pr_opt*self.cg*dni*pr_shw*pr_geo        
           
    def set_qu(self, urec, tro, tf) -> float:
        self.qu = urec*(tro-tf) #Ec. 3.21
    
    def urec(hint, dro, dri, krec) -> float:
        return 1/((1/hint) + (dro*np.log(dro/dri))/(2*krec)) #Ec. 3.22
    
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
    
    def set_pr(self, qabs):
        
        
        f0 = qabs/(self.urec*(self.tfe-self.text))   
        f1 = ((4*self.sigma*self.eext*self.text**3)+self.hext)/self.urec 
        f2 = 6*self.text**2*(self.sigma*self.eext/self.urec)*(self.qabs/self.urec)
        f3= 4*self.text*(self.sigma*self.eext/self.urec)*(self.qabs/self.urec)**2            
        f4 = (self.sigma*self.eext/self.urec)*(self.qabs/self.urec)
        
        pr0 = scy.optimize.newton(e3_36, pr0, fprime = e3_36prime, args=(f0,f1,f2,f3,f4))
        
        g = -(1+1/f0)+(1+f1)*z+f2*z**2+f3*z**3+f4*z**4
        g1 = f1+2*f2*z+3*f3*z**2+4*f4*z**3
        g2 = 2*f2+6*f3*z+12*f4*z**2
        g3 = 6*f3+24*f4*z
        
        self._pr = ((pr0*g1/(1-g1))*(1/(NTU*x))*(np.exp((1-g1)*NTU*x/g1)-1) -
                (g2/(6*g1))*(pr0*NTU*x)**2 -
                (g3/(24*g1)*(pr0*NTU*x)**3)
                )
        
        
        
        
        

class HCE(object):
    
    def __init__(self, sca, hce_order, parameters):
        
        self.sca = sca
        self.hce_order = hce_order
        self.parameters = parameters
        
    def get_index(self):
        return (self.sca.loop.solarfield.solarfield_order, 
                self.sca.loop.loop_order,
                self.sca.sca_order,
                self.hce_order)     
        
    def qu_from_pr(pr, qabs) -> float: 
        return pr*qabs #Ec. 3.26  
    
    def __PRBarbero4grade__() -> float:
        return 0
    
    def __PRBarbero1grade__() -> float:
        return 0
    
    def __PRBarbero0grade__() -> float:
        return 0   


class HCE_Barbero(HCE): 
  
    def __init__(self, parameters, sca, hce_order):
        
        self.eext = parameters['eext']
        self.hext = parameters['hext']
        self.hint = parameters['hext']
        self.urec = parameters['urec']
        self.sigma = parameters['sigma']
        self.cg = parameters['cg']
        self.pr_shw = parameters['pr_shw']
        self.pr_geo = parameters['pr_geo']
        self.pr_opt = parameters['pr_opt']
        
        HCE.__init__(self, sca, hce_order)
        

    def calcTempOut(tfe, dro, x, qabs, pr, massflow, cp) -> float:
        '''
        HTF output temperature [ºC] Ec. 3.24 Barbero

        Returns
        -------
        float [ºC]
            DESCRIPTION.

        '''
        return tfe + mt.PI*dro*x*qabs*pr/(massflow*cp)   

    def set_qabs(self, pr_opt, dni, pr_shw, pr_geo):
        '''
        Ec. 3.20 Barbero

        Returns
        -------
        float [W]
            DESCRIPTION.
            Thermal power absorbed by the HCE
        '''
        self.qabs = pr_opt*self.cg*dni*pr_shw*pr_geo        
           
    def set_qu(self, urec, tro, tf) -> float:
        self.qu = urec*(tro-tf) #Ec. 3.21
    
    def urec(hint, dro, dri, krec) -> float:
        return 1/((1/hint) + (dro*np.log(dro/dri))/(2*krec)) #Ec. 3.22
    
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
    
    def set_pr(self, qabs):
        
        
        f0 = qabs/(self.urec*(self.tfe-self.text))   
        f1 = ((4*self.sigma*self.eext*self.text**3)+self.hext)/self.urec 
        f2 = 6*self.text**2*(self.sigma*self.eext/self.urec)*(self.qabs/self.urec)
        f3= 4*self.text*(self.sigma*self.eext/self.urec)*(self.qabs/self.urec)**2            
        f4 = (self.sigma*self.eext/self.urec)*(self.qabs/self.urec)
        
        pr0 = scy.optimize.newton(e3_36, pr0, fprime = e3_36prime, args=(f0,f1,f2,f3,f4))
        
        g = -(1+1/f0)+(1+f1)*z+f2*z**2+f3*z**3+f4*z**4
        g1 = f1+2*f2*z+3*f3*z**2+4*f4*z**3
        g2 = 2*f2+6*f3*z+12*f4*z**2
        g3 = 6*f3+24*f4*z
        
        self._pr = ((pr0*g1/(1-g1))*(1/(NTU*x))*(np.exp((1-g1)*NTU*x/g1)-1) -
                (g2/(6*g1))*(pr0*NTU*x)**2 -
                (g3/(24*g1)*(pr0*NTU*x)**3)
                )
        
      
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
        self.solarfield = solarfield
        self.scas = []
        self.loop_order = loop_order
        #self.massFlow = massFlow
        
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
    
    def set_massFlow():
        '''
        Calculation of the current mass flow of HTF in the loop
        by dividing the total mass flow by the number of loops
        in the plant 
        

        Returns
        -------
        current HTF massflow, massFlowLoop

        '''

        self.massflow = (self.solarfield.massflow / 
                         self.solarfield.loops.len()
                        )
    
    def set_tin(self):
        
         self.tin = self.solarfield.tin
    
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
    
    
    def __init__(self, plant, solarfield_data):
        
        self.plant = plant
        self.name = solarfield_data.get('name')
        self.loops = []
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
        loop of the solar field
        '''
        H = 0.0
        
        for l in self.loops:
            H += l.get_cp * (l.get_tout-l.get_tin)*l.massflow
        
        self.tout = self.tin + H/(self.plant.hotfluid.get_cp()*self.get_massflow())
    
    def set_tin(self):
        
        self.tin = (self.plant.heatexchanger.hotfluid_tout *
                    self.coolPipeLosses
                    )
        
        
        
        
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
        
    def initializePlant(self, hE):
        '''
        Set initial values for some parameters
        
        Returns
        -------
        Nothing
        
        '''
        
        for sf in self.solarfields:
            for l in sf.loops:        
                l.set_massflow()
                l.set_tin()
        
        
        for s in self.solarfields:
            for l in s.loops:
                l.tin = s.coolPipeLosses*hE.hotfluid_tout
                
                
class HeatExchanger(object):
    
    def __init__(self, pr_heatexchanger, 
                 coolfluid_tin, 
                 hotfluid_tin, 
                 hotfluid,
                 coolfluid):
        
        self.pr_heatexchanger = pr_heatexchanger
        
    def set_fluid_tout(self):
        
        pass  
        
        
class Fluid(object):
    
    def __init__(self, name):
        
        self.name = name
        
        if self.name == 'Therminol VP1':            
            htf_data = pd.read_csv(monsanto.csv, sep=';', decimal=',', index_col=0)
        elif self.name == "Dowtherm":
            htf_data = pd.read_csv(dowtherm.csv, sep=',', decimal=',', index_col=0)
            
            
    def get_density(self):
        
        return self.density
    
    def get_cp(self):
        
        return self.cp        
        
    def set_fluid_cp(self, pressure, temperature):
        
        self.cp = 1000 # más adelante cp debe obtenerse con CoolProp, por ejemplo
        
class Site(object):
    def __init__(self, site):
        
        self._lat = site.get('latitude')
        self._long = site.get('longitude')
        
class HTF(object):
    
    def __init__(self, *argv):
        self._name = "Therminol VP-1"

class Weather(object):   
    def __init__(self):
        pass
        
    def loadWheatherData(path):
        pass

class ScatterMask(object):
    
    
    def __init__(self, simulation):
        
        self.matrix = dict()
   
        for sf in simulation.get('plant').get('solarfields'):
            self.matrix[sf["name"]]=[]
            for l in range(sf.get('loops')):
                self.matrix[sf["name"]].append([])
                for s in range(sf.get('scas')):
                    self.matrix[sf["name"]][-1].append([])
                    for h in range (sf.get('hces')):
                        self.matrix[sf["name"]][-1][-1].append(simulation["scattered_params"])

        
# format = '%Y-%m-%d %H:%M:%S'
# df['Datetime'] = pd.to_datetime(df['date'] + ' ' + df['time'], format=format)
# df = df.set_index(pd.DatetimeIndex(df['Datetime']))
        
    

class Simulation(object):
    ''' 
    Definimos la clase simulacion para representar las diferentes
    pruebas que lancemos, variando el archivo TMY, la configuracion del
    site, la planta, el modo de operacion o el modelo empleado. 
    
    '''
    
    def __init__(self):
        pass
        
        
    def calcRequiredMassFlow():
        '''calcula el caudal promedio requerido en cada lazo para alzanzar
        la temperatura deseada. Calcula el caudal en cada lazo (si son todos 
        iguales solo lo hace una vez) y después calcula el promedio en cada 
        subcampo. De esta forma tenemos el caudal '''
        
        pass
        
    def calc():
        pass
        return df