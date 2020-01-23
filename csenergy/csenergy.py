# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 08:26:43 2019

@author: pacomunuera
"""

import numpy as np
import scipy as sc
from scipy import constants
import math as mt
import CoolProp as CP
import time
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
    @classmethod 
    def set_tout(cls, hce, qabs, pr, cp):
        hce.tout = (hce.tin +sc.pi*hce.parameters['Dro']*hce.parameters['Long']*qabs*pr/
                    (hce.sca.loop.massflow*cp))
        
class ModelBarbero4grade(Model):

    def __ini__(self, simulation):
        
        super(Model, self).__init__(simulation)
        
    @classmethod
    def simulateHCE(cls, hce):
    
        dni = 800 #provisional
        cp = 2300.0 #provisional
        sigma = constants.sigma
        dro = hce.parameters['Dro']
        dri = hce.parameters['Dri']
        L = hce.parameters['Long']
        hint = hce.parameters['hint']
        hext = hce.parameters['hext']
        eext = hce.parameters['eext']
        krec = hce.parameters['krec']
        massflow = hce.sca.loop.massflow
        x= 1
        #cp = self.hot_fluid.get_cp(hce.tin)
        Model.set_tin(hce)
         # print("HCE", hce.get_index(), "hce.tin: ", hce.tin)
        tfe = hce.tin
        tf =  tfe
        hce.tout = tfe
        text = 22.0
        #Para el cálculo de la solución del rendimiento a la entrada puedes 
        #plantear como primera aproximación la solución de primer grado y 
        #aplicar le método de Newton-Raphson para obtener la solución, por ejemplo.
        
        pr0 = ModelBarbero1grade.simulateHCE(hce) 
        
                # Ec. 3.50 Barbero
        qcrit = sigma*eext*(tfe**4-text**4)+hext*(tfe-text)
        
        #Ec. 3.51 Barbero
        ucrit = 4*sigma*eext*tfe**3+hext
        #krec = (0.0153)*(trec) + 14.77 # trec ya está en ºC
        # Ec. 3.22
        urec = 1/(
            (1/hint) + 
            (dro*np.log(dro/dri))/(2*krec)
            )
        
        #Ec. 3.20 Barbero
        qabs = (hce.parameters['pr_opt']*
                hce.parameters['cg'] * dni *
                hce.parameters['pr_shw'] * 
                hce.parameters['pr_geo']) 
        
        pr0 = (1-qcrit/qabs)/(1+ucrit/urec)    
        
        # La forma de calcular sería suponer rendimiento térmico 1 y su
        pr1 = pr0
        # temperatura en pared sería Tro= Tf+q abs*rend/Urec. 

        # Con el rendimiento obtenido calculas la temperatura en pared y recalculas rendimeinto, 
        # así hasta convergencia.
        tro1 = tf+qabs*pr1/urec
        errtro = 1.
        errpr = 1.
        step = 0
        
        while (errtro > 0.0001 or errpr > 0.000001):
            
            step += 1
            #tro1 = tf+qabs*pr1/urec
            trec = (tro1+tf)/2
            #Ec. 4.22
            krec = (0.0153)*(trec) + 14.77 # trec ya está en ºC
            #Ec. 3.22
            urec = 1/(
                (1/hint) + 
                (dro*np.log(dro/dri))/(2*krec)
                )                     
            
            NTU = urec * x * L * sc.pi * dro / (massflow * cp)  #Ec. 3.30
                      
            f0 = qabs/(urec*(tfe-text))   
            f1 = ((4*sigma*eext*text**3)+hext)/urec 
            f2 = 6*(text**2)*(sigma*eext/urec)*(qabs/urec)
            f3 = 4*text*(sigma*eext/urec)*((qabs/urec)**2)            
            f4 = (sigma*eext/urec)*((qabs/urec)**3)
               
            fx=  lambda pr0: (1-pr0-
                            f1*(pr0+(1/f0))-
                            f2*((pr0+(1/f0))**2)-
                            f3*((pr0+(1/f0))**3)-
                            f4*((pr0+(1/f0))**4)
                           )
            
            dfx = lambda pr0: (-1-f1-
                             2*f2*(pr0+(1/f0))-
                             3*f3*(pr0+(1/f0))**2-
                             4*f4*(pr0+(1/f0))**3
                            )
    
            root = sc.optimize.newton(fx, 
                  pr0,
                  fprime = dfx,
                  maxiter = 10000)
            
            pr0 = root   
            z = pr0 + (1/f0)
            g1 = 1+f1+2*f2*z+3*f3*z**2+4*f4*z**3
            g2 = 2*f2+6*f3*z+12*f4*z**2
            g3 = 6*f3+24*f4*z        
            
            pr2 = ((pr0*g1/(1-g1))*(1/(NTU*x))*(sc.exp((1-g1)*NTU*x/g1)-1) -
                    (g2/(6*g1))*(pr0*NTU*x)**2 -
                    (g3/(24*g1)*(pr0*NTU*x)**3)
                    )

            tro2 = tf+qabs*pr2/urec         
            errpr = abs(pr2-pr1)
            errtro = abs(tro2-tro1)  
            tro1 = tro2
            pr1 = pr2            
        
        Model.set_tout(hce, qabs, pr1, cp )

    
    def e3_36(x, args) -> float:
    
        f0 = args[0]
        f1 = args[1]
        f2 = args[2]
        f3 = args[3]
        f4 = args[4]
        
        return  (1-x-
                 f1*(x+(1/f0))+
                 f2*((x+(1/f0))**2)+
                 f3*((x+(1/f0))**3)+
                 f4*((x+(1/f0))**4)
                 )
    
    def e3_36prime(x, args):
        
        f0 = args[0]
        f1 = args[1]
        f2 = args[2]
        f3 = args[3]
        f4 = args[4]
        
        return (-1-f1+
                2*f2*(x+1/f0)+
                3*f3*(x+1/f0)**2+
                4*f4*(x+1/f0)**3
                )
        
    def e3_100(self, tro,tf,text,hext,eext,sigma):
        
        return 3*tro**4-4*tf*tro**3+(text**4+hext*(text-tf)/(sigma*eext))
    
    def e3_100prime(self, tro,tf,text,hext,eext,sigma):
        
        return 12*tro**3-12*tf*tro**2
 

class ModelBarbero1grade(Model):

    def __ini__(self, simulation):
        
        super(Model, self).__init__(simulation)
        
    @classmethod
    def simulateHCE(cls, hce):
    
        dni = 800 #provisional
        cp = 2300.0 #provisional
        sigma = constants.sigma
        dro = hce.parameters['Dro']
        dri = hce.parameters['Dri']
        L = hce.parameters['Long']
        hint = hce.parameters['hint']
        hext = hce.parameters['hext']
        eext = hce.parameters['eext']
        krec = hce.parameters['krec']
        massflow = hce.sca.loop.massflow
        x= 1
        #cp = self.hot_fluid.get_cp(hce.tin)
        Model.set_tin(hce)
         # print("HCE", hce.get_index(), "hce.tin: ", hce.tin)
        tfe = hce.tin
        tf =  tfe
        hce.tout = tfe
        #Para el cálculo de la solución del rendimiento a la entrada puedes 
        #plantear como primera aproximación la solución de primer grado y 
        #aplicar le método de Newton-Raphson para obtener la solución, por ejemplo.
        
        pr0 = 1. 
        text = 22.0
        
        #Ec. 3.20 Barbero
        qabs = (hce.parameters['pr_opt']*
                hce.parameters['cg'] * dni *
                hce.parameters['pr_shw'] * 
                hce.parameters['pr_geo'])    
        
        # Ec. 3.50 Barbero
        qcrit = sigma*eext*(tfe**4-text**4)+hext*(tfe-text)
        
        #Ec. 3.51 Barbero
        ucrit = 4*sigma*eext*tfe**3+hext
        #krec = (0.0153)*(trec) + 14.77 # trec ya está en ºC
        # Ec. 3.22
        urec = 1/(
            (1/hint) + 
            (dro*np.log(dro/dri))/(2*krec)
            )
        
        
        hce.pr0 = (1-qcrit/qabs)/(1+ucrit/urec)
        
        
#        aext=0.5*sc.pi*dro
#        
#        qcrit = sigma*eext*(tfe**4-text**4)+hext*(tfe-text)
#        fcrit = 1/((4*sigma*eext*tfe**3/urec)+(hext/urec)+1)
#        ntuperd = (hext+4*sigma*eext*tfe**3)*aext/(massflow*cp)
#        hce.pr = (1-qcrit/qabs)*(1/(ntuperd*x))*(1-np.exp(-ntuperd*fcrit*x))
#        print("hce.pr: ", hce.pr)
#        hce.tout = tfe+sc.pi*dro*x*qabs*hce.pr/(massflow*cp)                
#        print("hce.tout: ", hce.tout)
        
class ModelBarbero0grade(Model):

    def __ini__(self, simulation):
        
        super(Model, self).__init__(simulation)
        
    def simulateHCE(self, hce):
    
        dni = 800 #provisional
                
        
        dro = hce.parameters['Dout']
        dri = hce.parameters['Din']
        x = hce.parameters['Long']
        hint = hce.parameters['hint']
        hext = hce.parameters['hext']
        sigma = hce.parameters['sigma']
        eext = hce.parameters['eext']
        krec = hce.parameters['krec']
        massflow = hce.sca.loop.massflow
        
        text = 22.
        tfe = hce.tin
        cp = 2300.0
        Model.set_tin(hce)
        print("hce.tin: ", hce.tin)
        
        qabs = (hce.parameters['pr_opt']*
        hce.parameters['cg'] * dni *
        hce.parameters['pr_shw'] * 
        hce.parameters['pr_geo']) 
        
        #krec = (0.0153)*(trec) + 14.77 # trec ya está en ºC
        #Ec. 3.22
        urec = 1/(
            (1/hint) + 
            (dro*np.log(dro/dri))/(2*krec)
            )  
        aext=0.5*sc.pi*dro
        
        qcrit = sigma*eext*(tfe**4-text**4)+hext*(tfe-text)
        fcrit = 1/((4*sigma*eext*tfe**3/urec)+(hext/urec)+1)
        ntuperd = (hext+4*sigma*eext*tfe**3)*aext/(massflow*cp)
        hce.pr = (1-qcrit/qabs)*(1/(ntuperd*x))*(1-np.exp(-ntuperd*fcrit*x))
        print("hce.pr: ", hce.pr)
        hce.tout = tfe+sc.pi*dro*x*qabs*hce.pr/(massflow*cp)                
        print("hce.tout: ", hce.tout)
        
        
class ModelBarberoSimplified(Model):

    def __ini__(self, simulation):
        
        super(Model, self).__init__(simulation)
        
    def simulateHCE(self, hce):
    
        dni = 800 #provisional
                
        
        dro = hce.parameters['Dout']
        dri = hce.parameters['Din']
        x = hce.parameters['Long']
        hint = hce.parameters['hint']
        hext = hce.parameters['hext']
        sigma = hce.parameters['sigma']
        eext = hce.parameters['eext']
        krec = hce.parameters['krec']
        massflow = hce.sca.loop.massflow
        
        text = 22.
        tfe = hce.tin
        cp = 2300.0
        Model.set_tin(hce)
        print("hce.tin: ", hce.tin)
        
        qabs = (hce.parameters['pr_opt']*
        hce.parameters['cg'] * dni *
        hce.parameters['pr_shw'] * 
        hce.parameters['pr_geo']) 
        
        #krec = (0.0153)*(trec) + 14.77 # trec ya está en ºC
        #Ec. 3.22
        urec = 1/(
            (1/hint) + 
            (dro*np.log(dro/dri))/(2*krec)
            )  
        aext=0.5*sc.pi*dro
        
        qcrit = sigma*eext*(tfe**4-text**4)+hext*(tfe-text)
        fcrit = 1/((4*sigma*eext*tfe**3/urec)+(hext/urec)+1)
        ntuperd = (hext+4*sigma*eext*tfe**3)*aext/(massflow*cp)
        hce.pr = (1-qcrit/qabs)*(1/(ntuperd*x))*(1-np.exp(-ntuperd*fcrit*x))
        print("hce.pr: ", hce.pr)
        hce.tout = tfe+sc.pi*dro*x*qabs*hce.pr/(massflow*cp)                
        print("hce.tout: ", hce.tout)
        
class ModelHottelWhilier(Model):
    
        def __ini__(self, simulation):
        
            super(Model, self).__init__(simulation)
        
        def simulateHCE(self, hce):
        
            dni = 800 #provisional
            
            dro = hce.parameters['Dout']
            dri = hce.parameters['Din']
            x = hce.parameters['Long']
            w = hce.parameters['w']
            hint = hce.parameters['hint']
            hext = hce.parameters['hext']
            sigma = hce.parameters['sigma']
            eext = hce.parameters['eext']
            krec = hce.parameters['krec']
            # En el modelo H-W uext es constante
            # uext = hce.parameters['uext']
            massflow = hce.sca.loop.massflow
            #cp = self.hot_fluid.get_cp(hce.tin)
            cp = 2300.0

            Model.set_tin(hce)
            print("hce.tin: ", hce.tin)
            time.sleep(0.010)
            tfe = hce.tin
            tf =  tfe
    
            text = 22.0
            tro = tf
            #Ec. 3.70
            #t3ro-ext = 4*tfe**3
            
            qabs = (hce.parameters['pr_opt']*
                hce.parameters['cg'] * dni *
                hce.parameters['pr_shw'] * 
                hce.parameters['pr_geo'])
            
            trec = tro
    
            #Ec. 4.22
            krec = (0.0153)*(trec) + 14.77 # trec ya está en ºC
                
            #Ec. 2.27

            #uext = sigma*eext*(tro**2+text**2)*(tro+text)+hext
            uext = 10
                
            #Ec. 3.22
            urec = 1/(
                    (1/hint) + 
                    (dro*np.log(dro/dri))/(2*krec)
                    )
                
            #Ec. 3.2
            Fprime = urec/(uext+urec)
            
                       

            print("---------------", (1-uext*(tfe-text)/qabs))

            pr = ((1-uext*(tfe-text)/qabs)*
                  massflow*cp/(w*x*uext)*
                  (1-np.exp(-Fprime*w*x*uext/(massflow*cp))))
            
            print("qabs ", qabs, "uext ", uext, "urec ", urec, "tro ", tro, "pr ", pr) 
            hce.pr = pr
            hce.tout = tfe+sc.pi*dro*x*qabs*pr/(massflow*cp)
            
        
        
class ModelNaumFrainderaich(Model):
    
        def __ini__(self, simulation):
        
            super(Model, self).__init__(simulation)

class ModelASHRAE(Model):

        def __ini__(self, simulation):
        
            super(Model, self).__init__(simulation)
            
        def simulateHCE(self, hce):
            
            dni = 800 #provisional
            dro = hce.parameters['Dout']
            dri = hce.parameters['Din']
            x = hce.parameters['Long']
            w = hce.parameters['w']
            hint = hce.parameters['hint']
            hext = hce.parameters['hext']
            sigma = hce.parameters['sigma']
            eext = hce.parameters['eext']
            krec = hce.parameters['krec']
            
            
            

        
class HCE(object):
    
    def __init__(self, sca, hce_order, simulation_settings):
        
        self.sca = sca
        self.hce_order = hce_order
        self.parameters = dict(simulation_settings.get('hce'),
                               **simulation_settings.get('model_settings'))
        self.tin = 0.0
        self.tout = 0.0
        self.pr = 0
    

    def get_previous(self):
        
        return self.sca.hces[self.hce_order-1]
    
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
        self.tin = 0.0
        self.massflow = 0.0


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
        
        Calculation of the massflow necessary to obtain the output temperature
        required by the operator

        Returns
        -------
        required HTF mass flow, reqMassFlowLoop

        '''
        self.required_massflow = 1
        


class SolarField(object):
    '''
    Parabolic Trough Solar Field
    
    '''
    
    def __init__(self, plant, simulation_settings):
        
        self.plant = plant
        self.name = simulation_settings['name']
        self.loops = []
        self.massflow =simulation_settings['massflow']
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
         
    def calcRequired_massflow(self):
        req_massflow = 0
        for l in self.loops:
            req_massflow += l.required_massflow
        self.requiered_massflow = req_massflow/self.loops.len()       
        
        
class SolarPlant(object):
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
        
        for sf in simulation_settings['plant']['solarfields']:
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

    def calcRequired_massflow(self):
        req_massflow = 0
        for sf in self.solarfields:
            req_massflow += sf.calc_required_massflow()
        self.req_massflow = req_massflow
        
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
                #sf.tin = self.exchanger.hot_fluid_Tout()
                sf.tin = self.tin
                for l in sf.loops:
                    l.massflow = sf.massflow/len(sf.loops)
                    l.tin = sf.tin
                    
                
    def print(self):
        
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
        
    def calcRequired_massflow():
        '''calcula el caudal promedio requerido en cada lazo para alzanzar
        la temperatura deseada. Calcula el caudal en cada lazo (si son todos 
        iguales solo lo hace una vez) y después calcula el promedio en cada 
        subcampo. De esta forma tenemos el caudal '''
        
        pass
    
           
class HeatExchanger(object):
    
    def __init__(self, hotfluid, coldfluid, pr_heatexchanger, 
                 hotfluidmassflow, coldfluidmassflow):
        
        self.pr_heatexchanger = pr_heatexchanger
        
    def set_fluid_tout(self):
        pass  
    
    def hot_fluid_tout():
        return 
    
class ThermodynamicCycle (object):
    
    def __init__(self, simulation_settings):
        
        self.name = simulation_settings['cycle']['name']
    
    @classmethod
    def get_NCA_pr(cls, tf, text):
        
        return (1-(text/tf)**(1/2))
        
    @classmethod
    def get_Carnot_pr(cls, tf, text):
        
        return (1-(text/tf))    
    
    
        
class Fluid(object):
    
    def __init__(self, name):
        pass
            
            
    def get_density(self):
        
        return self.density
    
    def get_cp(self):
        
        return self.cp
    
    def get_deltaH(self):
        
        cp_0 = hot_fluid['cp_factors'][0]
        cp_1 = hot_fluid['cp_factors'][1]
        cp_2 = hot_fluid['cp_factors'][2]
        cp_3 = hot_fluid['cp_factors'][3]
    
        return m*(cp_0*(tout-tin)+
                  (1/5)*cp_1*(tout**2-tin**2)+
                  (1/3)*cp_2*(tout**3-tin**3)+
                  (1/4)*cp_3*(tout**4-tin**4))
        
    def get_Reynolds(self):
        
        self.re = 0
        
    def get_ReynoldsDRI(self):
        
        self.redri = 0
        
    def get_nusselt_Dittus_Boelter(self):
        
        self.nudb = 0.023*(redri**0.8)*(prf**0.4) 
                           
    def get_nusselt_Gnielinski(self):
        
        self.nug = ((cf/2)*(redri-1000)*prf*(prf/prri)**0.11 /     
                    (1+12.7*(cf/2)**(1/2)*(prf**(2/3)-1))
                    )


class HotFluid(Fluid):
    
    def __init__(self, simulation_settings):
        
        self.name = simulation_settings['hot_fluid']['name']
        self.parameters = simulation_settings['hot_fluid']['parameters']
        
        
class ColdFluid(Fluid):
    
    def __init__(self, simulation_settings):
        
        self.name = simulation_settings['cold_fluid']['name']
        self.parameters = simulation_settings['cold_fluid']['parameters']
        

class Weather(object):   
    def __init__(self):
        pass
        
    def loadWheatherData(path):
        pass


        
# format = '%Y-%m-%d %H:%M:%S'
# df['Datetime'] = pd.to_datetime(df['date'] + ' ' + df['time'], format=format)
# df = df.set_index(pd.DatetimeIndex(df['Datetime']))
        
    

