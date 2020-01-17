# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 08:26:43 2019

@author: pacomunuera
"""

import numpy as np
import scipy as scy
import math as mt
import CoolProp as CP
import time
from math import pi
import pandas as pd

      
        
class Model(object):
    
    def __init__(self, simulation):
        self.model_settings = dict(simulation.get('model_settings'))
        self.hot_fluid = Fluid(dict(simulation.get('hot_fluid')))
        self.weather = dict(simulation.get('weather'))
    
    @classmethod
    def set_tin(cls, hce):
        if hce.hce_order > 0:
            hce.tin = hce.sca.hces[hce.hce_order-1].tout
        elif hce.sca.sca_order > 0:
            hce.tin = hce.sca.loop.scas[hce.sca.sca_order-1].hces[-1].tout
        else:
            hce.tin = hce.sca.loop.tin
    

                                 
        
class ModelBarbero4(Model):

    def __ini__(self, simulation):
        
        super(Model, self).__init__(simulation)
        
    def simulateHCE(self, hce):
    
        dni = 800
                
        
        dro = hce.parameters['Dout']
        dri = hce.parameters['Din']
        x = hce.parameters['Long']
        hint = hce.parameters['hint']
        hext = hce.parameters['hint']
        sigma = hce.parameters['sigma']
        eext = hce.parameters['eext']
        krec = hce.parameters['krec']
        massFlow = hce.sca.loop.massFlow
        #cp = self.hot_fluid.get_cp(hce.tin)
        cp = 2300.0
        tfe = hce.tin
        tf = tfe
        pr0=1
        pr=1

        text = 22.0
        
        #Ec. 3.20 Barbero
        qabs = (hce.parameters['pr_opt']*
                hce.parameters['cg'] * dni *
                hce.parameters['pr_shw'] * 
                hce.parameters['pr_geo'])  

        trec = tfe

        #Ec. 4.22
        krec = (0.0153)*(trec) + 14.77 # trec ya está en ºC
        
        #Ec. 3.22
        urec = 1/(
                (1/hint) + 
                (dro*np.log(dro/dri))/(2*krec)
                )
        
        tro = tf+qabs*pr/urec
        
        #Ec. 3.70
        #t3ro-ext = 4*tfe**3
        
        
#        tro= scy.optimize.newton(ModelBarbero4.e3_100, 
#                                  tro,
#                                  None,
#                                  args=(tf,text,hext,eext,sigma),
#                                  tol=0.0000001,
#                                  maxiter=1000)
        
        #Ec. 3.21
        qu = urec*(tro-tf) 

        #Ec. 3.30
        NTU = urec * x * pi * dro / (massFlow * cp)

        f0 = qabs/(urec*(tfe-text))   
        f1 = ((4*sigma*eext*(text**3))+hext)/urec 
        f2 = 6*(text**2)*(sigma*eext/urec)*(qabs/urec)
        f3= 4*text*(sigma*eext/urec)*((qabs/urec)**2)            
        f4 = (sigma*eext/urec)*((qabs/urec)**3)
        
        pr0 = 1
        
#        pr0 = scy.optimize.newton(ModelBarbero4.e3_36, 
#                                  pr0,
#                                  ModelBarbero4.e3_36prime,
#                                  args=(f0,f1,f2,f3,f4),
#                                  tol=0.01,
#                                  maxiter=1000)
        
        pr0 = scy.optimize.newton(ModelBarbero4.e3_36, 
                          pr0,
                          None,
                          args=(f0,f1,f2,f3,f4),
                          tol=0.01,
                          maxiter=1000)
        #fprime = ModelBarbero4.e3_36prime,

        z = pr0 + 1/f0
        
        # g = -(1+1/f0)+(1+f1)*z+f2*z**2+f3*z**3+f4*z**4
        g1 = f1+2*f2*z+3*f3*z**2+4*f4*z**3
        g2 = 2*f2+6*f3*z+12*f4*z**2
        g3 = 6*f3+24*f4*z        
        
        pr = ((pr0*g1/(1-g1))*(1/(NTU*x))*(np.exp((1-g1)*NTU*x/g1)-1) -
                (g2/(6*g1))*(pr0*NTU*x)**2 -
                (g3/(24*g1)*(pr0*NTU*x)**3)
                )
        
        hce.pr = pr
        
        hce.tout = tfe+pi*dro*x*qabs*pr/(massFlow*cp)
   
        print("hce", hce.get_index(),
              "Tout=", hce.tout,
              "tf=", tf, 
              "qabs=", qabs,
              "qu=",qu,
              "urec=",urec,
              "pr0=",pr0,
              "pr=",pr)    
    
    def e3_36(pr0,f0,f1,f2,f3,f4) -> float:
    
        return  (1-pr0-
                 f1*(pr0+(1/f0))+
                 f2*((pr0+(1/f0))**2)+
                 f3*((pr0+(1/f0))**3)+
                 f4*((pr0+(1/f0))**4)
                 )
    
    def e3_36prime(pr0,f0,f1,f2,f3,f4):
        
        return (-1-f1+
                2*f2*(pr0+1/f0)+
                3*f3*(pr0+1/f0)**2+
                4*f4*(pr0+1/f0)**3
                )
        
    def e3_100(tro,tf,text,hext,eext,sigma):
        
        return 3*tro**4-4*tf*tro**3+(text**4+hext*(text-tf)/(sigma*eext))
    
    def e3_100prime(tro,tf,text,hext,eext,sigma):
        
        return 12*tro**3-12*tf*tro**2
    
    
    

class ModelBarbero1(Model):

    def __init__(self, simulation):
        
        super(Model, self).__init__(simulation)
        
    def simulateHCE(self, hce):
    
        dni = 1000

        hce.pr = pr
        
        hce.tout = tfe+pi*dro*x*qabs*pr/(massFlow*cp)                
        
        
class ModelHottelWhilier(Model):
    
        def __ini__(self, simulation):
        
            super(Model, self).__init__(simulation)
        
        def simulateHCE(self, hce):
        

            dro = hce.parameters['Dout']
            dri = hce.parameters['Din']
            x = hce.parameters['Long']
            w = hce.parameters['x']
            hint = hce.parameters['hint']
            hext = hce.parameters['hint']
            sigma = hce.parameters['sigma']
            eext = hce.parameters['eext']
            krec = hce.parameters['krec']
            massFlow = hce.sca.loop.massFlow
            #cp = self.hot_fluid.get_cp(hce.tin)
            cp = 2300.0
            tfe = hce.tin
            tf = tfe
    
            text = 22.0
            tro = tf+100
            #Ec. 3.70
            #t3ro-ext = 4*tfe**3
            
            
            tro= scy.optimize.newton(ModelBarbero4.e3_100, 
                                      tro,
                                      None,
                                      args=(tf,text,hext,eext,sigma),
                                      tol=0.0000001,
                                      maxiter=1000)
            
            trec = tro
    
            #Ec. 4.22
            krec = (0.0153)*(trec) + 14.77 # trec ya está en ºC
                
            #Ec. 2.27
            uext = sigma*eext*(tro**2+text**2)*(tro+text)+hext
                
            #Ec. 3.22
            urec = 1/(
                    (1/hint) + 
                    (dro*np.log(dro/dri))/(2*krec)
                    )
                
            #Ec. 3.2
            Fprime = urec/(uext+urec)
            
            
            pr = ((1-uext*(tfe-text)/qabs)*
                  (massFlow*cp/(w*x*uext*(1-np.e(Fprime*w*x*uext/(massFlow*cp))))))
        
        
        
class ModelNaumFrainderaich(Model):
    
        def __ini__(self, simulation):
        
            super(Model, self).__init__(simulation)

class ModelASHRAE(Model):

        def __ini__(self, simulation):
        
            super(Model, self).__init__(simulation)    
    
        
class HCE(object):
    
    def __init__(self, sca, hce_order, simulation_settings):
        
        self.sca = sca
        self.hce_order = hce_order
        self.parameters = dict(simulation_settings.get('hce'),
                               **simulation_settings.get('model_settings'))
        self.tin = 0.0
        self.tout = 0.0
        self.pr = 0
        
    def get_index(self):
        return ([self.sca.loop.solarfield.name, 
                self.sca.loop.loop_order,
                self.sca.sca_order,
                self.hce_order])     
        
      
class SCA(object):
    
    SCA_configuration ={"HCE Number ": 12, "Position in Loop": 1,
                    "Temperature probes number": 1,
                    "Temp. probes position": 'Middle', "Defocus Order": 1}
    
# Uso de __slots__ permite ahorrar memoria RAM
#    __slots__= ['SCA_configuration']  
    def __init__(self, loop, sca_order):
        
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
            
class Loop(object):
    def __init__(self, solarfield, loop_order):
        self.solarfield = solarfield
        self.scas = []
        self.loop_order = loop_order
        self.massFlow = 0.0
        self.tin = 0.0
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
        self.massFlow =1
        self.tin= 1
        self.tout = self.tin
        
        
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
    
    def __init__(self, simulation_settings):

        self.name = simulation_settings['plant']['name']
        self.solarfields = []
        self.solarfield_to_exchanger_thermal_lost = 0.1
        self.exchanger_performance = 0.9
        self.exchanger_to_turbogroup_thermal_lost = 0.1
        self.steam_cycle_performance = 0.9
        self.turbogenerator_performance = 0.9
        
        self.tin = 300
        self.massFlow = 240
        
        for sf in simulation_settings.get('plant').get('solarfields'):
            self.solarfields.append(SolarField(self, sf))
            for l in range(sf.get('loops')):        
                self.solarfields[-1].loops.append(
                    Loop(self.solarfields[-1], l))
                for s in range(sf.get('scas')):
                    self.solarfields[-1].loops[-1].scas.append(
                        SCA(self.solarfields[-1].loops[-1], s))
                    for h in range (sf.get('hces')):
                        self.solarfields[-1].loops[-1].scas[-1].hces.append(
                            HCE(self.solarfields[-1].loops[-1].scas[-1],
                            h, 
                            simulation_settings))

        self.total_loops = 0
        
        for sf in self.solarfields:
            self.total_loops += len(sf.loops)

    def calcRequiredMassFlow(self):
        req_massflow = 0
        for sf in self.solarfields:
            req_massflow += sf.calcRequiredMassFlow()
        self.reqMassFlow = req_massflow
        
    def initializePlant(self, measures = None):
        '''
        Set initial values for some parameters
        
        Returns
        -------
        Nothing
        
        '''   
        if measures is not None:        
            for sf in self.solarfields:
                for l in sf.loops:        
                    l.set_inputs(measures)
        else:
            for sf in self.solarfields:
                sf.massFlow = len(sf.loops)*self.massFlow/self.total_loops
                #sf.tin = self.exchanger.hot_fluid_Tout()
                sf.tin = self.tin
                for l in sf.loops:
                    l.massFlow = sf.massFlow/len(sf.loops)
                    l.tin = sf.tin
                    
                
    def printPlant(self):
        
        for sf in self.solarfields:
            for l in sf.loops:
                for s in l.scas:
                    for h in s.hces:
                        print("Solarfield: ", sf.name, 
                              "Lazo: ",l.loop_order,
                              "SCA: ", s.sca_order,
                              "HCE: ", h.hce_order, 
                              "tin", "=", h.tin,
                              "tout", "=", h.tout)

        
class Site(object):
    def __init__(self, simulation_settings):
        
        self.name = simulation_settings.get('site').get('name')
        self.latitude = simulation_settings.get('site').get('latitude')
        self.longitude = simulation_settings.get('site').get('longitude')


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

    def applyMask(self, plant):
        
        for sf in plant.solarfields:
            for l in sf.loops:
                for s in l.scas:
                    for h in s.hces:
                        for k in self.matrix[sf.name][l.loop_order][s.sca_order][h.hce_order].keys():
                            h.parameters[k] *= float(self.matrix[sf.name][l.loop_order][s.sca_order][h.hce_order][k])        


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
    
           
class HeatExchanger(object):
    
    def __init__(self, pr_heatexchanger, 
                 coolfluid_tin, 
                 hotfluid_tin, 
                 hotfluid,
                 coolfluid):
        
        self.pr_heatexchanger = pr_heatexchanger
        
        self.hot_fluid = Fluid("HTF")
        self.cold_fluid = Fluid("water")
        
    def set_fluid_tout(self):
        pass  
    
    def hot_fluid_tout():
        return 

        
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


class HTF(object):
    
    def __init__(self, hot_fluid):
        
        self.name = hot_fluid.get('name')

    @classmethod
    def deltaH(cls, tin, tout, m, hot_fluid):
        
        cp_0 = hot_fluid['cp_factors'][0]
        cp_1 = hot_fluid['cp_factors'][1]
        cp_2 = hot_fluid['cp_factors'][2]
        cp_3 = hot_fluid['cp_factors'][3]
    
        return m*(cp_0*(tout-tin)+
                  (1/5)*cp_1*(tout**2-tin**2)+
                  (1/3)*cp_2*(tout**3-tin**3)+
                  (1/4)*cp_3*(tout**4-tin**4))
        

class Weather(object):   
    def __init__(self):
        pass
        
    def loadWheatherData(path):
        pass


        
# format = '%Y-%m-%d %H:%M:%S'
# df['Datetime'] = pd.to_datetime(df['date'] + ' ' + df['time'], format=format)
# df = df.set_index(pd.DatetimeIndex(df['Datetime']))
        
    

