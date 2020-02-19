# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 08:26:43 2019
@author: pacomunuera
"""

import numpy as np
import scipy as sc
from scipy import constants
import math as mt
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
import time
import pandas as pd
import pvlib as pvlib
from pvlib import iotools
from pvlib import solarposition
from pvlib import irradiance
from pvlib import iam
from tkinter import *
from tkinter.filedialog import askopenfilename
from datetime import datetime
import os.path
import matplotlib.pyplot as plt
import seaborn as sns


#  import pint

# ASHRAE: saturated liquid (Q=0) at T=233.15 K
_T_REF = 285.856

_IAM_PARAMS = {'EuroTrough ET150': {'F0': 1.0, 'F1': 0.0506, 'F2': -0.1763,
                                    'thetamin': 0, 'thetamax': 80},
               'Luz LS-2': {'F0': 1.0, 'F1': 0.0506, 'F2': -0.1763,
                            'thetamin': 0, 'thetamax': 80},
               'Luz LS-3': {'F0': 1.0, 'F1': 0.0506, 'F2': -0.1763,
                            'thetamin': 0, 'thetamax': 80},
               'Solargenix SGX-1': {'F0': 1.0, 'F1': 0.0506, 'F2': -0.1763,
                                    'thetamin': 0, 'thetamax': 80},
               'AlbiasaTrough AT150': {'F0': 1.0, 'F1': 0.0506, 'F2': -0.1763,
                                       'thetamin': 0, 'thetamax': 80},
               'Siemens SunField 6': {'F0': 1.0, 'F1': -0.0753, 'F2': -0.03698,
                                      'thetamin': 0, 'thetamax': 80},
               'SkyFuel SkyTrough': {'F0': 1.0, 'F1': 0.0327, 'F2': -0.1351,
                                     'thetamin': 0, 'thetamax': 80},
               'FLABEG Ultimate Trough RP6 (oil)':
                   {'F0': 1.0, 'F1': -0.005, 'F2': -0.102, 'thetamin': 0,
                    'thetamax': 80},
               'FLABEG Ultimate Trough RP6 (molten salt)':
                   {'F0': 1.0, 'F1': -0.008, 'F2': -0.117, 'thetamin': 0,
                    'thetamax': 80}}


class Model(object):

    def __init__(self, settings):
        self.model_settings = settings

    @classmethod
    def set_tin(cls, hce):
        if hce.hce_order > 0:
            hce.tin = hce.sca.hces[hce.hce_order-1].tout
        elif hce.sca.sca_order > 0:
            hce.tin = hce.sca.loop.scas[hce.sca.sca_order-1].hces[-1].tout
        else:
            hce.tin = hce.sca.loop.tin


class ModelBarbero4grade(Model):

    def __ini__(self, settings):

        super(Model, self).__init__(settings)

    @classmethod
    def set_pr(cls, hce, hotfluid, dni, wind, text, aoi):

        pressure = hce.sca.loop.pin
        Model.set_tin(hce)
        sigma = sc.constants.sigma
        dro = hce.parameters['dro']
        dri = hce.parameters['dri']
        dgo = hce.parameters['dgo']
        dgi = hce.parameters['dgi']
        L = hce.parameters['long']
        A = hce.sca.parameters['Aperture']
        IAM = hce.sca.get_IAM(aoi)
        pr_opt_peak = hce.get_pr_opt_peak()
        krec = hce.parameters['krec']
        massflow = hce.sca.loop.massflow
        cg = A /(np.pi*dro)
        x = 1
        tfe = hce.tin
        tf = tfe
        hce.tout = tf
        tro = tf
        tri = tf
        cp = hotfluid.get_cp(tf, pressure)

        #  nu_air = cinematic viscosity PROVISIONAL A FALTA DE VALIDAR TABLA
        # tabla en K
        nu_air = 8.678862e-11 * text**2 + 4.069284e-08 * text + -4.288741e-06

        #  Reynols number for wind at
        reext = dgo * wind / nu_air
        hext, eext = cls.get_hext_eext(hce, reext, tro, wind)

        # Ec. 4.14
        mu = hotfluid.get_dynamic_viscosity(tf, pressure)
        rho = hotfluid.get_density(tf, pressure)
        kf = hotfluid.get_thermal_conductivity(tf, pressure)
        alpha = kf / (rho * cp)
        prf = mu / alpha  #  Prandtl = viscosidad dinámica / difusividad_termica

        redri = 4 * massflow / (mu * np.pi * dri)  # Reynolds
        if redri <= 0:
            print("mu", mu, "massflow", massflow)
            
        cf = (1.58 * np.log(redri) - 3.28)**-2

        #  Prandtl number at temperature tri
        kfpri = hotfluid.get_thermal_conductivity(tri, pressure)
        rhori =  hotfluid.get_density(tri, pressure)
        cpri = hotfluid.get_cp(tri, pressure)
        alphari = kfpri / (rhori * cpri)
        muri = hotfluid.get_dynamic_viscosity(tri, pressure)
        prfri =  muri / alphari
        nug =((cf / 2)*(redri - 1000) * prfri * (prf / prfri)**0.11 /
              (1+12.7*(cf/2)**0.5 * (prf**(2/3) - 1)))

        hint = kf * nug / dri

        #  Ec. 3.50 Barbero
        qcrit = sigma * eext * (tfe**4 - text**4) + hext * (tfe - text)
        
        # Ec. 3.23
        qperd = sigma * eext * (tro**4 - text**4) + hext * (tro - text)

        #Ec. 3.51 Barbero
        ucrit = 4 * sigma * eext * tfe**3 + hext
        
        # Ec. 3.22
        urec = 1 / ((1 / hint) + ( dro * np.log(dro / dri)) / ( 2 * krec))

        #Ec. 3.20 Barbero
        qabs = (pr_opt_peak * IAM * cg * dni * hce.get_pr_shadows() *
                hce.get_pr_geo())

        if qabs > 0:
            pr0 = (1-qcrit/qabs)/(1+ucrit/urec)
        else:
            pr0 = 0.0

        pr1 = pr0

        tro1 = tf+qabs*pr1/urec
        errtro = 1.
        errpr = 1.
        step = 0

        if hce.tin >= hotfluid.tmax:
            qabs = 0 #defocused
            pr1 = 0

        if qabs > 0:

            while (errtro > 0.0001 or errpr > 0.000001):

                step += 1

                # Ec. 4.14
                cp = hotfluid.get_cp(tf, pressure)
                mu = hotfluid.get_dynamic_viscosity(tf, pressure)
                rho = hotfluid.get_density(tf, pressure)
                kf = hotfluid.get_thermal_conductivity(tf, pressure)
                #  alpha : difusividad térmica
                alpha = kf / (rho * cp)
                prf = mu / alpha  #  Prandtl = viscosidad dinámica / difusividad_termica
                redri = 4 * massflow / (mu * np.pi * dri)  # Reynolds
                #  Prandtl number at temperature tri
                kfpri = hotfluid.get_thermal_conductivity(tri, pressure)
                rhori =  hotfluid.get_density(tri, pressure)
                cpri = hotfluid.get_cp(tri, pressure)
                alphari = kfpri / (rhori * cpri)
                muri = hotfluid.get_dynamic_viscosity(tri, pressure)
                prfri =  muri / alphari
                nug =((cf / 2)*(redri - 1000) * prfri * (prf / prfri)**0.11 /
                      (1+12.7*(cf/2)**0.5 * (prf**(2/3) - 1)))

                #nudb = 0.023 * redri**0.8 * prf**0.4
                hint = kf * nug / dri

                hext, eext = cls.get_hext_eext(hce, reext, tro1, wind)
                trec = (tro1+tf)/2
                # Ec. 4.22
                krec = (0.0153)*(trec) + 14.77  # trec ya está en ºC
                # Ec. 3.22
                urec = 1/(
                    (1/hint) +
                    (dro*np.log(dro/dri))/(2*krec)
                    )

                NTU = urec * x * L * sc.pi * dro / (massflow * cp)  # Ec. 3.30

                f0 = qabs/(urec*(tfe-text))
                f1 = ((4*sigma*eext*text**3)+hext)/urec
                f2 = 6*(text**2)*(sigma*eext/urec)*(qabs/urec)
                f3 = 4*text*(sigma*eext/urec)*((qabs/urec)**2)
                f4 = (sigma*eext/urec)*((qabs/urec)**3)

                fx = lambda pr0: (1 - pr0 -
                                  f1*(pr0+(1/f0)) -
                                  f2*((pr0+(1/f0))**2) -
                                  f3*((pr0+(1/f0))**3) -
                                  f4*((pr0+(1/f0))**4))

                dfx = lambda pr0: (-1 - f1 -
                                   2*f2*(pr0+(1/f0)) -
                                   3*f3*(pr0+(1/f0))**2 -
                                   4*f4*(pr0+(1/f0))**3)

                root = sc.optimize.newton(fx,
                                          pr0,
                                          fprime=dfx,
                                          maxiter=100000)

                pr0 = root
                z = pr0 + (1/f0)
                g1 = 1+f1+2*f2*z+3*f3*z**2+4*f4*z**3
                g2 = 2*f2+6*f3*z+12*f4*z**2
                g3 = 6*f3+24*f4*z

                pr2 = ((pr0*g1/(1-g1))*(1/(NTU*x)) *
                       (sc.exp((1-g1)*NTU*x/g1)-1) -
                       (g2/(6*g1))*(pr0*NTU*x)**2 -
                       (g3/(24*g1)*(pr0*NTU*x)**3))

                tro2 = tf+qabs*pr2/urec
                errpr = abs(pr2-pr1)
                errtro = abs(tro2-tro1)
                tro1 = tro2
                pr1 = pr2
                hce.pr = pr1
                hce.set_tout(qabs, hotfluid)
                tf = hce.tout
                tri = tf
                hce.pr = pr1
                hce.set_tout(qabs, hotfluid)

        else:
            pr1 = 0
            hce.pr = pr1
            h1 = hotfluid.get_deltaH(hce.tin, hce.sca.loop.pin)
            h2 = qperd / hce.sca.loop.massflow
            hce.tout = hotfluid.get_T(h1 + h2, hce.sca.loop.pin)
            






    @classmethod
    def get_hext_eext(cls, hce, reext, tro, wind):

        eext = 0.
        hext = 0.

        if (hce.parameters['coating'] == 'CERMET' and
            hce.parameters['annulus'] == 'VACUUM'):
            if wind > 0:
                eext = 1.69E-4*reext**0.0395*tro+1/(11.72+3.45E-6*reext)
                hext = 0.
            else:
                eext = 2.44E-4*tro+0.0832
                hext = 0.
        elif (hce.parameters['coating'] == 'CERMET' and
              hce.parameters['annulus'] == 'NOVACUUM'):
            if wind > 0:
                eext = ((4.88E-10 * reext**0.0395 + 2.13E-4) * tro +
                        1 / (-36 - 1.29E-4 * reext) + 0.0962)
                hext = 2.34 * reext**0.0646
            else:
                eext = 1.97E-4 * tro + 0.0859
                hext = 3.65
        elif (hce.parameters['coating'] == 'BLACK CHROME' and
              hce.parameters['annulus'] == 'VACUUM'):
            if wind > 0:
                eext = (2.53E-4 * reext**0.0614 * tro +
                        1 / (9.92 + 1.5E-5 * reext))
                hext = 0.
            else:
                eext = 4.66E-4 * tro + 0.0903
                hext = 0.
        elif (hce.parameters['coating'] == 'BLACK CHROME' and
              hce.parameters['annulus'] == 'NOVACUUM'):
            if wind > 0:
                eext = ((4.33E-10 * reext + 3.46E-4) * tro +
                        1 / (-20.5 - 6.32E-4 * reext) + 0.149)
                hext = 2.77 * reext**0.043
            else:
                eext = 3.58E-4 * tro + 0.115
                hext = 3.6

        return hext, eext


    @classmethod
    def get_hint(cls, hce):

        hint = nbiot * krec / e

        pass


class ModelBarbero1grade(Model):

    def __ini__(self, simulation):

        super(Model, self).__init__(simulation)

    @classmethod
    def set_pr(cls, hce, hot_fluid, dni):

        dni = dni
        cp = hot_fluid.get_cp(hce.tin)
        sigma = sc.constants.sigma
        dro = hce.parameters['dro']
        dri = hce.parameters['din']
        dgo = hce.parameters['dgo']
        dgi = hce.parameters['dgi']
        L = hce.parameters['long']
        hint = hce.parameters['hint']
        hext = hce.parameters['hext']
        eext = hce.parameters['eext']
        krec = hce.parameters['krec']
        massflow = hce.sca.loop.massflow
        x = 1
        # cp = self.hot_fluid.get_cp(hce.tin)
        Model.set_tin(hce)

        tfe = hce.tin
        tf = tfe
        hce.tout = tfe

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

        fcrit = 1/(1+ucrit/urec)

        aext=0.5*sc.pi*dro
        ntuperd = ucrit*aext/(massflow*cp)

        if qabs > 0.0:
            pr = (1-qcrit/qabs)*(1/(ntuperd*x))*(1-sc.exp(-ntuperd*fcrit*x))
        else:
            pr = 0.0
#
        Model.set_tout(hce, qabs, pr, cp )


class ModelBarberoSimplified(Model):

    def __ini__(self, simulation):

        super(Model, self).__init__(simulation)

    @classmethod
    def set_pr(cls, hce, hot_fluid, dni):

        dni = dni
        cp = hot_fluid.get_cp(hce.tin)
        sigma = sc.constants.sigma
        dro = hce.parameters['dro']
        dri = hce.parameters['din']
        dgo = hce.parameters['dgo']
        dgi = hce.parameters['dgi']
        L = hce.parameters['long']
        hint = hce.parameters['hint']
        hext = hce.parameters['hext']
        eext = hce.parameters['eext']
        krec = hce.parameters['krec']
        massflow = hce.sca.loop.massflow
        cg = hce.sca.parameters['aperture']/(np.pi*hce.parameters['Dro'])
        x = 1
        # cp = self.hot_fluid.get_cp(hce.tin)
        Model.set_tin(hce)
        tfe = hce.tin
        tf = tfe
        hce.tout = tfe

        text = 22.0

        # Ec. 3.20 Barbero
        qabs = (hce.parameters['pr_opt'] *
                cg * dni *
                hce.parameters['pr_shw'] *
                hce.parameters['pr_geo'])

        # Ec. 3.50 Barbero
        qcrit = sigma*eext*(tfe**4-text**4)+hext*(tfe-text)

        # Ec. 3.51 Barbero
        ucrit = 4*sigma*eext*tfe**3+hext
        # krec = (0.0153)*(trec) + 14.77 # trec ya está en ºC
        # Ec. 3.22
        urec = 1/(
            (1/hint) +
            (dro*np.log(dro/dri))/(2*krec)
            )

        fcrit = 1/(1+ucrit/urec)

        if qabs > 0.0:
            pr = fcrit*(1-qcrit/qabs)
        else:
            pr = 0.0

        Model.set_tout(hce, qabs, pr, cp)


class ModelHottelWhilier(Model):

        def __ini__(self, simulation):

            super(Model, self).__init__(simulation)


class ModelNaumFrainderaich(Model):

        def __ini__(self, simulation):
            super(Model, self).__init__(simulation)


class ModelASHRAE(Model):

        def __ini__(self, simulation):

            super(Model, self).__init__(simulation)



class HCE(object):

    def __init__(self, sca, hce_order, settings, model_settings):

        self.sca = sca
        self.hce_order = hce_order
#        self.parameters = dict(settings.get('hce'),
#                               **settings.get('model_settings'))
        self.parameters = dict(settings, **model_settings)
        self.tin = 0.0
        self.tout = 0.0
        self.pr = 0.0
        #self.PR = 0.0
        self.qabs = 0.0


    def set_tout(self, qabs, hotfluid):
        
        h1 = hotfluid.get_deltaH(self.tin, self.sca.loop.pin)
        h2 = (np.pi * self.parameters['dro'] * self.parameters['long'] * 
              qabs * self.pr / self.sca.loop.massflow)
        
        self.tout = hotfluid.get_T(h1 + h2, self.sca.loop.pin)

        # self.tout = (self.tin + np.pi * self.parameters['dro'] *
        #             self.parameters['long'] * qabs * self.pr /
        #             (self.sca.loop.massflow * hotfluid.get_cp(self.tin,
        #                                                       self.sca.loop.pin)))

        # if (self.hce_order == 11 and self.sca.sca_order == 3 and self.sca.loop.loop_order == 29):
        #     print(qabs)

    def get_previous(self):

        return self.sca.hces[self.hce_order-1]

    def get_index(self):
        return ([self.sca.loop.solarfield.name,
                self.sca.loop.loop_order,
                self.sca.sca_order,
                self.hce_order])


    def get_pr_opt_peak(self):

        alpha = self.get_absorptance()
        tau = self.get_transmissivity()
        rho = self.sca.parameters['Reflectance']
        gamma = self.get_solar_fraction()

        return alpha * tau * rho * gamma


    def get_pr_opt(self):
        return 1.0


    def get_pr_geo(self):
        return 1.0


    def get_absorptance(self):

        #
        alpha = 0.96

        return alpha


    def get_transmissivity(self):

        #dependiendo si hay vídro y si hay vacío
        tau = 0.963

        return tau


    def get_reflectance(self):

        return self.sca.parameters['Reflectance']


    def get_solar_fraction(self):

        return 1.0

    def get_pr_shadows(self):
        return 1.0


    def get_pr_total(self, dateindex, site, data):

        aoi = self.sca.get_aoi(dateindex, site, data)

        return (self.pr * self.get_pr_shadows() *
                self.get_pr_opt_peak() * self.get_pr_geo() *
                np.cos(np.radians(aoi)) * self.sca.get_IAM(aoi))


class SCA(object):

# Uso de __slots__ permite ahorrar memoria RAM
#    __slots__= ['SCA_configuration']
    def __init__(self, loop, sca_order, settings):

        self.loop = loop
        self.sca_order = sca_order
        self.hces = []
        self.tracking_angle = 0.0
        self.parameters = dict(settings)
        self.surface_tilt = 0.0
        self.surface_azimuth = 180.0

    def get_IAM(self, theta):

        theta = np.radians(theta)
        F0 = self.parameters['IAM Coefficient F0']
        F1 = self.parameters['IAM Coefficient F1']
        F2 = self.parameters['IAM Coefficient F2']

        return F0+(F1*theta+F2*theta**2)/np.cos(theta)

    # def get_tin(self, sca):

    #     if not sca:
    #         sca = self
    #     if sca.sca_order > 0:
    #         return SCA.get_tout(sca.loop.scas[sca.sca_order-1])
    #     else:
    #         return Loop.get_tin(sca.loop)

    # def get_tout(self, sca):
    #     if not sca:
    #         sca = self
    #     return HCE.get_tout(sca.hces[-1])

    def get_aoi(self, dateindex, site, data):

        solarposition = pvlib.solarposition.get_solarposition(dateindex,
                                                        site.latitude,
                                                        site.longitude,
                                                        site.altitude,
                                                        pressure = data['Pressure'],
                                                        temperature=data['DryBulb'])

        aoi = float(pvlib.irradiance.aoi(self.surface_tilt,
                                               self.surface_azimuth,
                                               solarposition.zenith[0],
                                               solarposition.azimuth[0]))

        return aoi


class Loop(object):

    def __init__(self, solarfield, loop_order):
        self.solarfield = solarfield
        self.scas = []
        self.loop_order = loop_order
        self.tin = 0.0
        self.massflow = 0.0
        self.tout = 0.0
        self.pin = 0.0
        self.pout = 0.0
        self.cut_tout = 0.0
        self.wasted_power = 0.0

    def initialize(self):

        self.tin = self.solarfield.tin
        self.massflow = self.solarfield.massflow / len(self.solarfield.loops)
        self.pin = self.solarfield.pin
        self.pout = self.solarfield.pout

    # def get_tin(self):

    #     return self.scas[0].hces[0].tin


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


    def precalmassflow(self):
        '''

        Calculation of the massflow necessary to obtain the output temperature
        required by the operator

        Returns
        -------
        required HTF mass flow, reqMassFlowLoop

        '''
        pass


class PrototypeLoop(object):

    def __init__(self, plant_settings, sca_settings, hce_settings, model_settings):

        self.scas = []
        self.loop_order = 0
        self.plant_settings = plant_settings
        self.sca_settings = sca_settings
        self.hce_settings = hce_settings
        self.solarplant = None
        self.tin = 0.0
        self.tout = 0.0
        self.pin = 0.0
        self.pout = 0.0
        self.cut_tout = 0.0
        self.req_massflow = 0.0
        self.wasted_power = 0.0
        self.min_massflow = 0.0



        for s in range(self.plant_settings['loop']['scas']):
            self.scas.append(SCA(self, s, sca_settings))
            for h in range(self.plant_settings['loop']['hces']):
                self.scas[-1].hces.append(HCE(self.scas[-1], h ,
                         hce_settings, model_settings))

    def initialize(self, solarplant):

        self.solarplant = solarplant

    def precalcmassflow(self, row, aoi, hotfluid, model):

        self.tin = self.solarplant.tin
        self.ratedtout = self.solarplant.ratedtout
        self.pin = self.solarplant.pin
        self.pout = self.solarplant.pout
        self.massflow = self.solarplant.massflow / self.solarplant.total_loops
        self.hotfluid = hotfluid

        self.min_massflow = self.hotfluid.get_massflow_from_Reynolds(self.hce_settings['dri'],
                                   self.tin, self.pin,
                                   self.hce_settings['min_reynolds'])

        inc_ratio = 1.01
        dec_ratio = 0.99
        max_error = 1.0

        search = True

        while search:

            for s in self.scas:
                for h in s.hces:
                    model.set_pr(h, self.hotfluid, row[1]['DNI'],
                                          row[1]['Wspd'], row[1]['DryBulb'], aoi)
                    #print("HCE", h.hce_order, h.tout)
            self.tout = self.scas[-1].hces[-1].tout

            err = round(100 * abs(self.tout-self.ratedtout)/self.ratedtout,2)

            if err > max_error:

                if self.tout >= self.ratedtout:
                    self.massflow *= (1 + err/100)
                    search = True
                elif (self.massflow > self.min_massflow):
                    self.massflow *= (1 - err/100)
                    search = True
                else:
                    self.massflow = self.min_massflow
                    search = False
            else:
                search = False

        return self.massflow


class SolarField(object):
    '''
    Parabolic Trough Solar Field

    '''

    def __init__(self, plant, solarfield_settings, loop_settings):

        self.plant = plant
        self.name = solarfield_settings['name']
        self.loops = []
        self.massflow = 0.0
        self.tin= 0.0
        self.tout = 0.0
        self.cut_tout = 0.0

    def set_massflow(self):

        mf = 0.0
        for l in self.loops:
            mf += l.massflow

        self.massflow = mf


    # def get_tout(self, hotfluid):
    #     '''
    #     Calculates HTF output temperature throughout the solar field as a
    #     weighted average based on the enthalpy of the mass flow in each
    #     loop of the solar field
    #     '''
    #     H = 0.0

    #     HH = 0.0

    #     pressure = 2000000

    #     for l in self.loops:
    #         H += (hotfluid.get_cp(l.scas[-1].hces[-1].tout, pressure) *
    #               (l.scas[-1].hces[-1].tout-l.get_tin())*l.massflow)

    #         HH += PropsSI('T','H',HH ,'P', l.pin, hotfluid.coolpropID)

    #     HH /= l
    #     self.tout = (self.tin + H /
    #                  (hotfluid.get_cp(l.scas[-1].hces[-1].tout, pressure) *
    #                   self.get_massflow()))

    #     return self.tout


    def set_tout(self, hotfluid):
        '''
        Calculates HTF output temperature throughout the solar field as a
        weighted average based on the enthalpy of the mass flow in each
        loop of the solar field
        '''
        # H = 0.0
        # H2 = 0.0
        dH = 0.0
        dH_cut_tout = 0.0

        for l in self.loops:

            # H += (hotfluid.get_cp(l.tout, l.pout) *
            #       (l.tout-l.tin)*l.massflow)
            # H2  += (hotfluid.get_cp(l.cut_tout, l.pout) *
            #       (l.cut_tout-l.tin)*l.massflow)

            dH += l.massflow * hotfluid.get_deltaH(l.cut_tout, l.pout)
            dH_cut_tout += l.massflow * hotfluid.get_deltaH(l.tout, l.pout)

        dH /= self.massflow
        dH_cut_tout /= self.massflow
        self.tout =  hotfluid.get_T(dH, self.pout)
        self.cut_tout = hotfluid.get_T(dH_cut_tout, self.pout)

        # self.tout = (self.tin + H /
        #              (hotfluid.get_cp(l.tout, l.pout) *
        #               self.massflow))
        # self.cut_tout =  (self.tin + H2 /
        #              (hotfluid.get_cp(l.cut_tout, l.pout) *
        #               self.massflow))

    def set_tin(self):

        self.tin = (self.plant.heatexchanger.hotfluid_tout *
                    self.coolPipeLosses
                    )

    def apply_temp_limitation(self, hotfluid):

        for l in self.loops:
            if l.tout > hotfluid.tmax:
                l.cut_tout = hotfluid.tmax
                l.wasted_power = l.massflow * hotfluid.get_cp(
                        l.cut_tout, l.pout)*(l.tout - l.cut_tout)
            else:
                l.cut_tout = l.tout
                l.wasted_power = 0.0

    def loops_avg_out(self):

        tavg = 0.0
        cont = 0

        for l in self.loops:
            cont += 1
            tavg += l.tout

        return tavg / cont

    def initialize(self, row = None):

        if row is not None:
            self.massflow = row[1][self.name+'.mf']
            self.tin = row[1][self.name+'.tin']
            self.pin = 100000 * row[1][self.name+'.pin']
            self.pout = 100000 * row[1][self.name+'.pout']
        else:
            self.massflow = (self.plant.ratedmassflow * len(self.loops) /
                             self.plant.total_loops)
            self.tin = self.plant.tin
            self.pin = self.plant.pin
            self.pout = self.plant.pout

    def calcRequired_massflow(self):
        req_massflow = 0
        for l in self.loops:
            req_massflow += l.required_massflow
        self.requiered_massflow = req_massflow/self.loops.len()


class SolarPlant(object):
    '''
    Parabolic Trough Concentrated Solar Power Plant

    '''

    def __init__(self, plant_settings, sca_settings, hce_settings,
                 hce_model_settings):

        self.name = plant_settings['name']
        self.solarfields = []
        self.ratedtin = plant_settings['ratedtin']
        self.ratedtout = plant_settings['ratedtout']
        self.ratedpressure = plant_settings['ratedpressure']
        self.min_massflow = 0.0
        self.tout = 0.0
        self.cut_tout = 0.0
        self.tin = self.ratedtin
        self.pin = plant_settings['ratedpressure']
        self.pout = plant_settings['ratedpressure']
        self.massflow = 0.0
        self.ratedmassflow = plant_settings['ratedmassflow']
        self.hce_settings = hce_settings


        for sf in plant_settings['solarfields']:
            self.solarfields.append(SolarField(self, sf, plant_settings['loop']))
            for l in range(sf.get('loops')):
                self.solarfields[-1].loops.append(
                    Loop(self.solarfields[-1], l))
                for s in range(plant_settings['loop']['scas']):
                    self.solarfields[-1].loops[-1].scas.append(
                        SCA(self.solarfields[-1].loops[-1],
                            s,
                            sca_settings))
                    for h in range (plant_settings['loop']['hces']):
                        self.solarfields[-1].loops[-1].scas[-1].hces.append(
                            HCE(self.solarfields[-1].loops[-1].scas[-1],
                            h,
                            hce_settings,
                            hce_model_settings))

        self.total_loops = 0

        for sf in self.solarfields:
            self.total_loops += len(sf.loops)

#        proto_settings = {'name': 'prototype',
#                          'massflow': 0}
#
#        self.protosolarfield = SolarField(self, proto_settings, plant_settings['loop'])
#
#        for l in range(self.total_loops):
#                self.protosolarfield.loops.append(
#                    Loop(self.protosolarfield, l))
#                for s in range(plant_settings['loop']['scas']):
#                    self.protosolarfield.loops[-1].scas.append(
#                        SCA(self.protosolarfield.loops[-1],
#                            s,
#                            sca_settings))
#                    for h in range (plant_settings['loop']['hces']):
#                        self.protosolarfield.loops[-1].scas[-1].hces.append(
#                            HCE(self.protosolarfield.loops[-1].scas[-1],
#                            h,
#                            hce_settings,
#                            hce_model_settings))


    def set_tout(self, hotfluid):
        '''
        Calculates HTF output temperature throughout the solar plant as a
        weighted average based on the enthalpy of the mass flow in each
        loop of the solar field
        '''
        # H = 0.0
        # H2 = 0.0

        dH = 0.0
        dH_cut_tout = 0.0

        for sf in self.solarfields:

            # H += (hotfluid.get_cp(sf.tin, sf.pin) *
            #       (sf.tout - sf.tin)*sf.massflow)
            # H2 += (hotfluid.get_cp(sf.tin, sf.pin) *
            #       (sf.cut_tout - sf.tin)*sf.massflow)

            dH += sf.massflow * hotfluid.get_deltaH(sf.tout, sf.pout)
            dH_cut_tout += sf.massflow * hotfluid.get_deltaH(sf.cut_tout, sf.pout)

        dH /= self.massflow
        dH_cut_tout /= self.massflow
        self.tout = hotfluid.get_T(dH, self.pout)
        self.cut_tout = hotfluid.get_T(dH_cut_tout, self.pout)

        # self.tout = (self.tin + (H /
        #              (hotfluid.get_cp(self.tin, self.pin) * self.massflow)))

        # self.cut_tout = (self.tin + (H2 /
        #      (hotfluid.get_cp(self.tin, self.pin) * self.massflow)))

    def set_tin(self, hotfluid):
        '''
        Calculates HTF output temperature throughout the solar plant as a
        weighted average based on the enthalpy of the mass flow in each
        loop of the solar field
        '''
        # H = 0.0
        dH = 0.0

        # TT = 0.0

        for sf in self.solarfields:

            # H += (hotfluid.get_cp(sf.tin, sf.pin) *
            #       (sf.tin- 488.0)*sf.massflow)

            dH += (hotfluid.get_deltaH(sf.tin, sf.pin) * sf.massflow)

            # TT += sf.tin * sf.massflow

        # self.tin = (488.0 + (H /
        #              (hotfluid.get_cp(488.0, self.pin) * self.massflow)))

        dH /= self.massflow

        self.tin = hotfluid.get_T(dH, self.pin)
        #self.tin = 1e-3 * HH / self.massflow

    def set_massflow(self):

        mf = 0.0

        for sf in self.solarfields:
            mf += sf.massflow

        self.massflow = mf


    def calcRequired_massflow(self):
        req_massflow = 0
        for sf in self.solarfields:
            req_massflow += sf.calc_required_massflow()
        self.req_massflow = req_massflow

    def initializePlant(self, simulation, datasource = None):
        '''
        Set initial values for some parameters

        Returns
        -------
        Nothing

        '''
        if simulation.type == "type0":
            for sf in self.solarfields:
                for l in sf.loops:
                    l.set_inputs(measures)
        elif (simulation.type == "type1" and datasource is not None):
            pass
            # for sf in self.solarfields:
            #     sf.tin = datasource[sf.name+'.tin'][0]
            #     sf.massflow = datasource[sf.name+'.mf'][0]
            #     sf.pout = datasource[sf.name+'.pout'][0]
            #     sf.pin = datasource[sf.name+'.pin'][0]
            #     sf.cut_tout = self.cut_tout
            #     for l in sf.loops:
            #         l.massflow = sf.massflow/len(sf.loops)
            #         l.tin = sf.tin
            #         l.pin = sf.pin
            #         l.pout = sf.pout
            #         l.cut_tout = sf.cut_tout


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


#class SolarSystem(object):
#
#    pass

#    def __init__(self, settings):
#
#        self.parameters = settings
#
#
#    #  Duffie - Beckman calculation of the incidence angle
#    @classmethod
#    def aoi_DB(self, sca, index):
#
#        delta = 0
#        phi = 0
#        beta = 0
#        gamma_s = 0
#        omega = 0
#
#        return  np.arccos(
#                (np.sin(delta) * np.sin(phi) * np.sin(beta) -
#                np.sin(delta) * np.cos(phi) * np.cos(gamma_s) +
#                np.cos(delta) * np.cos(phi) * np.cos(beta) * cos(omega) +
#                np.cos(delta) * np.sin(phi) * np.sin(beta) *
#                np.cos(gamma_s) * np.cos(omega) +
#                np.cos(delta) * np.sen(beta) * np.sen(gamma_s) *
#                np.sen(omega)))
#
#
#    #  Rabl calculation of the incidence angle with N-S axis
#    @classmethod
#    def aoi_Rabl_NS(self, sca, index):
#
#        delta = 0
#        phi = 0
#        omega = 0
#
#        return np.arccos(np.cos(delta) *
#                          np.sqrt(
#                                  (np.cos(phi)*np.cos(omega)+
#                                   np.tan(delta)*np.sin(phi))**2 +
#                        (np.sin(omega))**2))
#
#    #  Rabl calculation of the incidence angle with E-W axis
#    @classmethod
#    def aoi_Rabl_EW(self, sca, index):
#
#        delta = 0
#        omega = 0
#
#        return np.arccos(np.sqrt(1 +
#                                 (np.cos(delta)**2) *
#                                 ((np.cos(omega))**2 - 1)))


class PowerSystem(object):
    '''
    Power Plant as a set composed by a SolarPlant, a HeatExchanger, a PowerCycle and a BOPSystem

    '''

    # Calculo de potencia de salida (positiva o negativa)
    #y potencia derivada a BOB en base al estado de la planta

    def __init__(self, settings):

        self.ratedpower = settings['ratedpower']
        self.solarfield_to_exchanger_pr = settings['solarfield_to_exchanger_pr']
        self.exchanger_pr = settings['exchanger_pr']
        self.exchanger_to_turbogroup_pr = settings['exchanger_to_turbogroup_pr']
        self.steam_cycle_pr = settings['steam_cycle_pr']
        self.turbogenerator_pr = settings['turbogenerator_pr']


    def get_poweroutput(self):
        pass

    def get_powerinput(self):
        pass

    def get_thermalpower(self, ratedpower):

        return ratedpower / (self.solarfield_to_exchanger_pr *
                                  self.exchanger_pr *
                                  self.exchanger_to_turbogroup_pr *
                                  self.steam_cycle_pr *
                                  self.turbogenerator_pr)

    def get_ratedmassflow(self, ratedpower, solarplant, hotfluid, coldfluid = None):

        cp_avg = 0.5 * (hotfluid.get_cp(solarplant.ratedtin, solarplant.ratedpressure) +
                        hotfluid.get_cp(solarplant.ratedtout, solarplant.ratedpressure))

        ratedmassflow = (self.get_thermalpower(self.ratedpower) /
               (cp_avg * (solarplant.ratedtout - solarplant.ratedtin)))

        return ratedmassflow


class BOPSystem(object):
    '''
    BOP: Balance Of Plant System.
    The BOP encompasses the electrical consumptions for auxiliary systems
    (pumping, tracing, compressed air)

    '''
    #BOP debe calcular con consumo de potencia según datos de
    #campo solar, iintercambiador y ciclo (bombeos) y datos
    #meteorologicos y de operación (consumo de traceados, aire comprimido...)

    def __init__(self, settings):

        self.name = settings['bopsystem']['name']

    def get_powerinput_from_PowerPlant(self):
        pass

    def get_powerinput_from_PowerGrid(self):
        pass

class PowerCycle(object):

    def __init__(self, settings):

        self.name = settings['name']
        self.pin = settings['pin']
        self.tin = settings['tin']
        self.pout = 0.0
        self.tout = 0.0
        self.pr = 0.0

    def set_pr_NCA(self, tdewpoint): #Novikov-Curzon-Ahlbor

        self.pr = 1 - np.sqrt((tdewpoint+273.15)/(self.tin + 273.15))

class Generator(object):

    def __init__(self, settings):

        self.name = settings['name']
        self.pr = settings['pr']

    def set_pr(self):

        #TO-DO: Desasrrollar curvas pr-carga-temp.ext por ejemplo.
        pass



class Site(object):
    def __init__(self, settings):

        self.name = settings['name']
        self.latitude = settings['latitude']
        self.longitude = settings['longitude']
        self.altitude = settings['altitude']


class HCEScatterMask(object):


    def __init__(self, plant_settings, hce_mask_settings):

        self.matrix = dict()

        for sf in plant_settings['solarfields']:
            self.matrix[sf["name"]]=[]
            for l in range(sf.get('loops')):
                self.matrix[sf["name"]].append([])
                for s in range(plant_settings['loop']['scas']):
                    self.matrix[sf["name"]][-1].append([])
                    for h in range (plant_settings['loop']['hces']):
                        self.matrix[sf["name"]][-1][-1].append(hce_mask_settings)

    def applyMask(self, plant):

        for sf in plant.solarfields:
            for l in sf.loops:
                for s in l.scas:
                    for h in s.hces:
                        for k in self.matrix[sf.name][l.loop_order][s.sca_order][h.hce_order].keys():
                            h.parameters[k] *= float(self.matrix[sf.name][l.loop_order][s.sca_order][h.hce_order][k])


class SCAScatterMask(object):


    def __init__(self, plant_settings, sca_mask_settings):

        self.matrix = dict()

        for sf in plant_settings['solarfields']:
            self.matrix[sf["name"]]=[]
            for l in range(sf.get('loops')):
                self.matrix[sf["name"]].append([])
                for s in range(plant_settings['loop']['scas']):
                    self.matrix[sf["name"]][-1].append(sca_mask_settings)

    def applyMask(self, plant):

        for sf in plant.solarfields:
            for l in sf.loops:
                for s in l.scas:
                    for k in self.matrix[sf.name][l.loop_order][s.sca_order].keys():
                        s.parameters[k] *= float(self.matrix[sf.name][l.loop_order][s.sca_order][k])


class Simulation(object):
    '''
    Definimos la clase simulacion para representar las diferentes
    pruebas que lancemos, variando el archivo TMY, la configuracion del
    site, la planta, el modo de operacion o el modelo empleado.
    '''

    def __init__(self, settings):
        self.ID =  settings['ID']
        self.type = settings['type']
        self.solarplant = None
        self.powersystem = None
        self.hotfluid = None
        self.coldfluid = None
        self.site = None
        self.model = None
        self.datasource = None
        self.powercycle = None

    def precalc(self, hce_settings):

        loop_min_massflow = self.hotfluid.get_massflow_from_Reynolds(
                hce_settings['dri'],
                self.solarplant.ratedtin,
                self.solarplant.ratedpressure,
                hce_settings['min_reynolds'])

        print("Loop_min_mass_flow", loop_min_massflow)

        self.solarplant.min_massflow = self.solarplant.total_loops * loop_min_massflow

        self.solarplant.ratedmassflow = self.powersystem.get_ratedmassflow(
                self.powersystem.ratedpower, self.solarplant, self.hotfluid)

#        fluid_speed = 4 * loop_min_massflow / ( np.pi * hce_settings['dri']**2 *
#                self.hotfluid.get_density(self.solarplant.ratedtin,
#                                     self.solarplant.ratedpressure))

        if self.solarplant.min_massflow > self.solarplant.ratedmassflow:
            print("Too low massflow", self.solarplant.min_massflow ,">",
                  self.solarplant.ratedmassflow)
        else:
            print("Rated massflow = ", self.solarplant.ratedmassflow, ">",
                  "Min massflow=", self.solarplant.min_massflow)

    def runSimulation(self):

        if self.type == "type0":
            self.simulateSolarPlant()
        elif self.type == "type1":
            self.benchmarkSolarPlant()
        else:
            return None


    def simulateSolarPlant(self):

        if self.type != 'type0':
            return None


        self.precalc(self.solarplant.hce_settings)

        self.solarplant.massflow = self.solarplant.ratedmassflow

        r_tout = []
        r_powerout = []
        r_massflow = []
        r_netpower = []
        r_cycle_pr = []
        r_heatexchanger_pr = []


        for row in self.datasource.weatherdata[0].iterrows():

            solarpos = pvlib.solarposition.get_solarposition(row[0],
                                                        self.site.latitude,
                                                        self.site.longitude,
                                                        self.site.altitude,
                                                        pressure = row[1]['Pressure'],
                                                        temperature=row[1]['DryBulb'])

            aoi = float(pvlib.irradiance.aoi(0, 0, solarpos['zenith'][0],
                                             solarpos['azimuth'][0]))

            for sf in self.solarplant.solarfields:
                sf.initialize()
                for l in sf.loops:
                    l.initialize()
                    self.prototypeloop.initialize(self.solarplant)
                    l.massflow = self.prototypeloop.precalcmassflow(
                            row, aoi, self.hotfluid, self.model)
                    for s in l.scas:
                        for h in s.hces:
                            self.model.set_pr(h, self.hotfluid, row[1]['DNI'],
                                              row[1]['Wspd'], row[1]['DryBulb'], aoi)
                    l.tout = l.scas[-1].hces[-1].tout
                sf.set_massflow()
                sf.apply_temp_limitation(self.hotfluid)
                sf.set_tout(self.hotfluid)


            self.solarplant.set_massflow()
            self.solarplant.set_tin(self.hotfluid)
            self.solarplant.set_tout(self.hotfluid)

            self.heatexchanger.set_hotfluid_in(self.solarplant.massflow,
                                             self.solarplant.tin,
                                             self.solarplant.pin)

            #self.heatexchanger.set_coldfluid_in() PENDIENTE

            self.heatexchanger.set_thermalpowertransfered(self.solarplant.tin)
            self.powercycle.set_pr_NCA(288.15) #Provisonal temperatura agua coondensador
            self.generator.set_pr()
            r_cycle_pr.append(round(self.powercycle.pr,2))
            r_heatexchanger_pr.append(round(self.heatexchanger.pr,2))

            r_tout.append(round(self.solarplant.tout, 0))
            r_powerout.append(round((self.solarplant.tout - self.solarplant.tin) *
                          self.solarplant.massflow * 1e-6 *
                          self.hotfluid.get_cp(
                                  self.solarplant.tout, self.solarplant.pout),0))
            r_massflow.append(round(self.solarplant.massflow,0))
            r_netpower.append(round(
                    1e-6 * self.heatexchanger.thermalpowertransfered *
                    self.powercycle.pr * self.generator.pr, 2))

            print(row[0].strftime('%y/%m/%d %H:%M'), 'DNI: ', round(row[1]['DNI']),
              "MF:", round(self.solarplant.massflow), 'Calc:',
              round(self.solarplant.tin), '->', round(self.solarplant.tout))

        self.datasource.weatherdata[0]['Powerout'] = r_powerout
        self.datasource.weatherdata[0]['Tout'] = r_tout
        self.datasource.weatherdata[0]['MF'] = r_massflow
        self.datasource.weatherdata[0]['NetPower'] = r_netpower
        self.datasource.weatherdata[0]['Cycle.pr'] = r_cycle_pr
        self.datasource.weatherdata[0]['HE.pr'] = r_heatexchanger_pr

        print(self.datasource.weatherdata[0])
        self.datasource.weatherdata[0][[
                'DNI','Powerout','Tout','MF', 'NetPower', 'Cycle.pr',
                'HE.pr']].plot(figsize=(20,10), linewidth=5, fontsize=20)
        plt.xlabel('Date', fontsize=20)



    def benchmarkSolarPlant(self):

        if self.type != 'type1':
            return None

        self.solarplant.initializePlant(self, self.datasource.dataframe)

        r_tout = []
        r_powerout = []
        r_massflow = []
        r_netpower = []
        r_cycle_pr = []
        r_heatexchanger_pr = []

        for row in self.datasource.dataframe.iterrows():
            solarpos = pvlib.solarposition.get_solarposition(row[0],
                                                        self.site.latitude,
                                                        self.site.longitude,
                                                        self.site.altitude,
                                                        pressure = row[1]['Pressure'],
                                                        temperature=row[1]['DryBulb'])
            aoi = float(pvlib.irradiance.aoi(0, 0, solarpos['zenith'][0],
                                             solarpos['azimuth'][0]))
            for sf in self.solarplant.solarfields:
                sf.initialize(row)
                for l in sf.loops:
                    l.initialize()
                    for s in l.scas:
                        for h in s.hces:
                            self.model.set_pr(h, self.hotfluid,
                                                   row[1]['DNI'],
                                                   row[1]['Wspd'],
                                                   row[1]['DryBulb'], aoi)
#                            h.set_PR(hce, aoi)
#                            self.model.set_tout(h, qabs, h.pr, cp)
                    l.tout = l.scas[-1].hces[-1].tout

                sf.apply_temp_limitation(self.hotfluid)
                sf.set_tout(self.hotfluid)

            self.solarplant.set_massflow()
            self.solarplant.set_tin(self.hotfluid)
            self.solarplant.set_tout(self.hotfluid)

            self.heatexchanger.set_hotfluid_in(self.solarplant.massflow,
                                             self.solarplant.tin,
                                             self.solarplant.tout)
            self.heatexchanger.set_thermalpowertransfered(self.solarplant.tin)
            self.powercycle.set_pr_NCA(288.15) #Provisonal temperatura agua coondensador
            self.generator.set_pr()
            r_cycle_pr.append(round(self.powercycle.pr,2))
            r_heatexchanger_pr.append(round(self.heatexchanger.pr,2))

            r_tout.append(round(self.solarplant.tout, 0))
            r_powerout.append(round((self.solarplant.tout - self.solarplant.tin) *
                          self.solarplant.massflow * 1e-6 *
                          self.hotfluid.get_cp(
                                  self.solarplant.tout, self.solarplant.pout),0))
            r_massflow.append(round(self.solarplant.massflow,0))
            r_netpower.append(round(
                    1e-6 * self.heatexchanger.thermalpowertransfered *
                    self.powercycle.pr * self.generator.pr, 2))

            print(row[0].strftime('%y/%m/%d %H:%M'), 'DNI: ', round(row[1]['DNI']),
                  "MF:", round(self.solarplant.massflow), 'Real:',
                  round(row[1][sf.name+'.tin']), '->',
                  round(row[1][sf.name+'.tout']), 'Calc:', round(sf.tin), '->',
                  round(self.solarplant.tout))

        self.datasource.dataframe['Powerout'] = r_powerout
        self.datasource.dataframe['Tout'] = r_tout
        self.datasource.dataframe['MF'] = r_massflow
        self.datasource.dataframe['NetPower'] = r_netpower
        self.datasource.dataframe['Cycle.pr'] = r_cycle_pr
        self.datasource.dataframe['HE.pr'] = r_heatexchanger_pr

        print(self.datasource.dataframe)
        self.datasource.dataframe[[
                'DNI','Powerout','Tout','MF', 'NetPower', 'Cycle.pr',
                'HE.pr']].plot(figsize=(20,10), linewidth=5, fontsize=20)
        plt.xlabel('Date', fontsize=20)

        #TO-DO
        estimated_pr = 0.0
        actual_pr = 0.0
        estimated_tout = 0.0
        actual_tout = 0.0
        rejected_solar_energy = 0.0


class HeatExchanger(object):

    def __init__(self, settings, hotfluid, coldfluid):

        self.pr = settings['pr']
        self.thermalpowertransfered = 0.0

        self.hotfluidmassflow = 0.0
        self.hotfluidtin = 0.0
        self.hotfluidtout = 0.0
        self.hotfluidpin = 0.0
        self.hotfluidpout = 0.0
        self.hotfluid = hotfluid

        self.coldfluidmassflow = 0.0
        self.coldfluidtin = 0.0
        self.coldfluidtout = 0.0
        self.colfluidpin = 0.0
        self.colfluidtout = 0.0
        self.coldfluid = coldfluid

    def set_hotfluid_in(self, hotfluidmassflow, hotfluidtin, hotfluidpin):

        self.hotfluidmassflow = hotfluidmassflow
        self.hotfluidtin = hotfluidtin
        self.hotfluidpin = hotfluidpin

    def set_coldfluid_in(self, coldfluidmassflow, coldfluidtin, coldfluidpin):

        self.coldfluidmassflow = coldfluidmassflow
        self.coldfluidtin = coldfluidtin
        self.coldfluidpin = coldfluidpin

    def set_thermalpowertransfered(self, hftout):

#        self.thermalpowertransfer = (
#                self.pr *
#                self.hotfluidmassflow *
#                self.hotfluid.get_cp(self.hotfluid.tin, self.hotfluid.pin) *
#                self.hotfluidtin -self.hotfluidtout)

        self.thermalpowertransfered = (
                self.pr *
                self.hotfluidmassflow *
                self.hotfluid.get_cp(self.hotfluidtin, self.hotfluidpin) *
                self.hotfluidtin - hftout)


    def hot_fluid_tout():
        return


class HeatStorage(object):

    def __init__(self, settings):

        self.name = settings['heatstorage']

    def set_fluid_tout(self):
        pass

    def hot_fluid_tout():
        return


class Fluid_CoolProp(object):

    def __init__(self, settings = None):

        self.name = settings['name']
        self.tmax = settings['tmax']
        self.tmin = settings['tmin']
        self.coolpropID = settings['CoolPropID']

    def get_density(self, t, p):

        return PropsSI('D','T',t,'P', p, self.coolpropID)

    def get_dynamic_viscosity(self, t, p):

        return  PropsSI('V','T',t,'P', p, self.coolpropID)

    def get_cp(self, t, p):

        return PropsSI('C','T',t,'P', p, self.coolpropID)

    def get_thermal_conductivity(self, t, p):
        ''' Saturated Fluid conductivity at temperature t '''

        return PropsSI('L','T',t,'P', p, self.coolpropID)

    def get_deltaH(self, t, p):

        CP.set_reference_state(self.coolpropID,'ASHRAE')

        deltaH = PropsSI('H','T',t ,'P', p, self.coolpropID)

        CP.set_reference_state(self.coolpropID, 'DEF')

        return deltaH

    def get_T(self, h, p):

        CP.set_reference_state(self.coolpropID,'ASHRAE')
        temperature = PropsSI('T', 'H', h, 'P', p, self.coolpropID)
        CP.set_reference_state(self.coolpropID, 'DEF')

        return temperature


    def get_Reynolds(self, dri, t, p, massflow):

        return massflow * np.pi * (dri**3) /( 4 * self.get_density(t, p))

    def get_massflow_from_Reynolds(self, dri, t, p, re):

        return re * np.pi * dri * self.get_dynamic_viscosity(t, p) / 4

    def get_ReynoldsDRI(self):

        self.redri = 0

    def get_Nusselt_Dittus_Boelter(self):

        self.nudb = 0.023*(redri**0.8)*(prf**0.4)

    def get_Nusselt_Gnielinski(self):

        self.nug = ((cf/2)*(redri-1000)*prf*(prf/prri)**0.11 /
                    (1+12.7*(cf/2)**(1/2)*(prf**(2/3)-1))
                    )

class Fluid_Tabular(object):

    def __init__(self, settings=None):

        self.name = settings['name']
        self.cp = settings['cp']
        self.rho = settings['rho']
        self.mu = settings['mu']
        self.kt = settings['kt']
        self.h = settings['h']
        self.t = settings['t']
        self.tmax = settings['tmax']
        self.tmin = settings['tmin']

        self.cp += [0.] * (6 - len(self.cp))
        self.rho += [0.] * (6 - len(self.rho))
        self.mu += [0.] * (6 - len(self.mu))
        self.kt += [0.] * (6 - len(self.kt))

    def get_density(self, t, p):

        rho0, rho1, rho2, rho3, rho4, rho5 = tuple(self.rho)

        return (rho0 + rho1 * t + rho2 * t**2 + rho3 * t**3 +
                rho4 * t**4 + rho5 * t**5)

    def get_dynamic_viscosity(self, t, p):

        mu0, mu1, mu2, mu3, mu4, mu5 = tuple(self.mu)

        return (mu0 + mu1 * t + mu2 * t**2 + mu3 * t**3 +
                mu4 * t**4 + mu5 * t**5)

    def get_cp(self, t, p):

        cp0, cp1, cp2, cp3, cp4, cp5 = tuple(self.cp)

        return (cp0 + cp1 * t + cp2 * t**2 + cp3 * t**3 +
                cp4 * t**4 + cp5 * t**5)

    def get_thermal_conductivity(self, t, p):
        ''' Saturated Fluid conductivity at temperature t '''

        kt0, kt1, kt2, kt3, kt4, kt5 = tuple(self.kt)

        return (kt0 + kt1 * t + kt2 * t**2 + kt3 * t**3 +
                kt4 * t**4 + kt5 * t**5)

    def get_deltaH(self, t, p):

        h0, h1, h2, h3, h4, h5 = tuple(self.h)

        href =(h0 + h1 * _T_REF + h2 * _T_REF**2 + h3 * _T_REF**3 +
               h4 * _T_REF**4 + h5 * _T_REF**5)

        h = (h0 + h1 * t + h2 * t**2 + h3 * t**3 + h4 * t**4 + h5 * t**5)
        return (h - href)

    def get_T(self, h, p):

        t0, t1, t2, t3, t4, t5 = tuple(self.t)

        return (t0 + t1 * h + t2 * h**2 + t3 * h**3 +
                t4 * h**4 + t5 * h**5)

    def get_Reynolds(self, dri, t, p, massflow):

        return (4 * massflow /
                (np.pi * dri * self.get_dynamic_viscosity(t,p)))

    def get_massflow_from_Reynolds(self, dri, t, p, re):

        return re * np.pi * dri * self.get_dynamic_viscosity(t,p) / 4

    def get_ReynoldsDRI(self):

        self.redri = 0

    def get_Nusselt_Dittus_Boelter(self):

        self.nudb = 0.023*(redri**0.8)*(prf**0.4)

    def get_Nusselt_Gnielinski(self):

        self.nug = ((cf/2)*(redri-1000)*prf*(prf/prri)**0.11 /
                    (1+12.7*(cf/2)**(1/2)*(prf**(2/3)-1))
                    )


class Weather(object):

    def __init__(self, settings):
        self.filename = settings['filename']
        self.filepath = settings['filepath']
        self.file = self.filepath + self.filename
        self.weatherdata = None
        self.openWeatherDataFile(self.file)


        self.change_units()
        self.filter_columns()

    def openWeatherDataFile(self, path = None):


        try:
            if path is None:
                root = Tk()
                root.withdraw()
                path = askopenfilename(initialdir = ".weather_files/",
                                   title = "choose your file",
                                   filetypes = (("TMY files","*.tm2"),
                                                ("TMY files","*.tm3"),
                                                ("csv files","*.csv"),
                                                ("all files","*.*")))
                root.update()
                root.destroy()

                if path is None:
                    return
                else:
                    strfilename, strext = os.path.splitext(path)

                    if  strext == ".csv":
                        self.weatherdata = pvlib.iotools.tmy.read_tmy3(path)
                        self.file = path
                    elif (strext == ".tm2" or strext == ".tmy"):
                        self.weatherdata = pvlib.iotools.tmy.read_tmy2(path)
                        self.file = path
                    elif strext == ".xls":
                        pass
                    else:
                        print("unknow extension ", strext)
                        return

            else:
                strfilename, strext = os.path.splitext(path)

                if  strext == ".csv":
                    self.weatherdata = pvlib.iotools.tmy.read_tmy3(path)
                    self.file = path
                elif (strext == ".tm2" or strext == ".tmy"):
                    self.weatherdata = pvlib.iotools.tmy.read_tmy2(path)
                    self.file = path
                elif strext == ".xls":
                    pass
                else:
                    print("unknow extension ", strext)
                    return
        except Exception:
            raise
            txMessageBox.showerror('Error loading Weather Data File',
                                   'Unable to open file: %r', self.file)

    def change_units(self):

        for c in self.weatherdata[0].columns:
            if (c == 'DryBulb') or (c == 'DewPoint'): # From Celsius Degrees to K
                self.weatherdata[0][c] *= 0.1
                self.weatherdata[0][c] += 273.15
            if c=='Pressure': # from mbar to Pa
                self.weatherdata[0][c] *= 1e2

    def filter_columns(self):

        needed_columns = ['DNI', 'DryBulb', 'DewPoint', 'Wspd', 'Wdir', 'Pressure']
        columns_to_drop = []
        for c in  self.weatherdata[0].columns:
            if c not in needed_columns:
                columns_to_drop.append(c)
        self.weatherdata[0].drop(columns = columns_to_drop, inplace = True)

#    def get_weather_data_site(self):
#
#        return self.weatherdata[1]
#
#    def loadWeatherDataFile(self):
#
#        self.weatherdata = pvlib.iotools.tmy.read_tmy2(self.file)


class FieldData(object):

    def __init__(self, settings):
        self.filename = settings['filename']
        self.filepath = settings['filepath']
        self.file = self.filepath + self.filename
        self.tags = settings['tags']
        self.dataframe = None

        self.openFieldDataFile(self.file)
        self.rename_columns()
        self.change_units()


    def openFieldDataFile(self, path = None):

        '''
        fielddata
        '''
        dateparse = lambda x: pd.datetime.strptime(x, '%Y/%m/%d %H:%M')

        try:
            if path is None:
                root = Tk()
                root.withdraw()
                path = askopenfilename(initialdir = ".fielddata_files/",
                                   title = "choose your file",
                                   filetypes = (("csv files","*.csv"),
                                                ("all files","*.*")))
                root.update()
                root.destroy()

                if path is None:
                    return
                else:
                    strfilename, strext = os.path.splitext(path)

                    if  strext == ".csv":
                        print("csv........")

                        self.dataframe = pd.read_csv(path, sep=';',
                                                     decimal= ',',
                                                     dtype= float,
                                                     parse_dates=['datetime'],
                                                     date_parser=dateparse,
                                                     index_col=0)
                        self.file = path
                    elif strext == ".xls":
                        print("xls...")
                        self.dataframe = pd.read_excel(path)
                        self.file = path
                    else:
                        print("unknow extension ", strext)
                        return
            else:
                strfilename, strext = os.path.splitext(path)

                if  strext == ".csv":
                    print("csv...")
                    self.dataframe = pd.read_csv(path, sep=';',
                                                 decimal= ',',
                                                 dayfirst=True,
                                                 index_col=0)
                    self.file = path
                elif strext == ".xls":
                    print("xls...")
                    self.dataframe = pd.read_excel(path)
                    self.file = path
                else:
                    print("unknow extension ", strext)
                    return

        except Exception:
            raise
            txMessageBox.showerror('Error loading FieldData File',
                                   'Unable to open file: %r', self.file)

        self.dataframe.index = pd.to_datetime(self.dataframe.index)

    def change_units(self):

        for c in self.dataframe.columns:
            if ('.t' in c) or ('DryBulb' in c) or ('Dew' in c):
                self.dataframe[c] += 273.15 # From Celsius Degrees to K
            if '.p' in c:
                self.dataframe[c] *= 1e5 # From Bar to Pa

    def rename_columns(self):

        rename_dict = dict(zip(self.tags.values(), self.tags.keys()))
        self.dataframe.rename(columns = rename_dict, inplace = True)
        columns_to_drop = []
        for c in  self.dataframe.columns:
            if c not in self.tags.keys():
                columns_to_drop.append(c)
        self.dataframe.drop(columns = columns_to_drop, inplace = True)








# format = '%Y-%m-%d %H:%M:%S'
# df['Datetime'] = pd.to_datetime(df['date'] + ' ' + df['time'], format=format)
# df = df.set_index(pd.DatetimeIndex(df['Datetime']))



