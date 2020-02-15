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

#  import pint

_FLUIDS_PARAMS = {'DOWTHERM A': {'cp0': 1522.77738, 'cp1': 2.59864, 'cp2': 0.00046,
                                 'v0': 1.44502E-02, 'v1': -8.32838E-05,
                                 'v2': 1.38092E-06, 'v3': -5.82052E-09,
                                 'v4': 1.17142E-11, 'v5': -8.04330E-15,
                                 'd0': 1071.13711 , 'd1': -0.66794, 'd2': -0.00072,
                                 'k0': 1.856e-1 , 'k1': -1.600e-4 , 'k2': 5.913e-12,
                                 'tmax': 400, 'tmin': 15},
                 'SYLTHERM 800': {'cp0': 1574.36919, 'cp1': 1.71019, 'cp2': -0.00001,
                                  'v0': 1.5223e-2 , 'v1': -2.63195e-4,
                                  'v2': 2.09535e-6, 'v3': -8.615274e-9,
                                  'v4': 1.75478e-11, 'v5': -1.396471e-14,
                                  'd0': 947.65387 , 'd1': -0.74514, 'd2': -0.00061,
                                  'k0': 1.387691e-1 , 'k1':- 1.88054e-4 , 'k2':-1.14165e-10,
                                  'tmax': 400, 'tmin': -40},
                 'THERMINOL VP1': {'cp0': 1499.03529, 'cp1': 2.71399, 'cp2': -0.00009,
                                   'v0': 6.1345e-3 , 'v1': -1.1974e-4,
                                   'v2': 1.0466e-6, 'v3': -4.5758e-9,
                                   'v4': 9.6995e-12, 'v5': -7.9178e-15,
                                   'd0': 1074.61508 , 'd1': -0.67848, 'd2': -0.00063,
                                   'k0': 1.381091e-1 , 'k1':-8.7078e-5 , 'k2':-1.7298e-7,
                                   'tmax': 400, 'tmin': 15}}

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

    @classmethod
    def set_tout(cls, hce, qabs, pr, cp):

        hce.tout = (hce.tin +
                    sc.pi*hce.parameters['dro'] *
                    hce.parameters['long'] * qabs * pr /
                    (hce.sca.loop.massflow * cp))


class ModelBarbero4grade(Model):

    def __ini__(self, settings):

        super(Model, self).__init__(settings)


    @classmethod
    def simulateHCE(cls, hce, hot_fluid, dni, wind, text):

        pressure = 1500000

        Model.set_tin(hce)
        sigma = sc.constants.sigma
        dro = hce.parameters['dro']
        dri = hce.parameters['dri']
        dgo = hce.parameters['dgo']
        dgi = hce.parameters['dgi']
        L = hce.parameters['long']
        A = hce.sca.parameters['Aperture']
        krec = hce.parameters['krec']
        massflow = hce.sca.loop.massflow
        cg = A /(np.pi*dro)
        x = 1
        tfe = hce.tin
        tf = tfe
        hce.tout = tf
        tro = tf
        tri = tf
        cp = hot_fluid.get_cp(tf, pressure)

        #  nu_air = cinematic viscosity PROVISIONAL A FALTA DE VALIDAR TABLA
        nu_air = 8.67886e-11 * text**2 + 8.81055e-8 * text + 1.33019e-5

        #  Reynols number for wind at
        reext = dgo * wind / nu_air
        hext, eext = cls.get_hext_eext(hce, reext, tro, wind)

        #hint = kf * nuint / dri
        #hint = Model.get_hint(hce)
        #hint = hce.parameters['hint']

        #  mu viscosidad dinámica (Pa·s)


        # mal, reext es para el exteior
        #reext = 4 * massflow / (np.pi * dri * mu)
        #air = CP.HAProps('W', 'T', 300, 'P', 101.325, 'R', 0.5)


        # Ec. 4.14
        mu = hot_fluid.get_dynamic_viscosity(tf, pressure)
        rho = hot_fluid.get_density(tf, pressure)
        kf = hot_fluid.get_thermal_conductivity(tf, pressure)
        #  alpha : difusividad térmica
        alpha = kf / (rho * cp)
        prf = mu / alpha  #  Prandtl = viscosidad dinámica / difusividad_termica
        redri = 4 * massflow / (mu * np.pi * dri)  # Reynolds
        # nudb = 0.023 * redri**0.8 * prf**0.4
        cf = (1.58 * np.log(redri) - 3.28)**-2

        #  Prandtl number at temperature tri
        kfpri = hot_fluid.get_thermal_conductivity(tri, pressure)
        rhori =  hot_fluid.get_density(tri, pressure)
        cpri = hot_fluid.get_cp(tri, pressure)
        alphari = kfpri / (rhori * cpri)
        muri = hot_fluid.get_dynamic_viscosity(tri, pressure)
        prfri =  muri / alphari
        nug =((cf / 2)*(redri - 1000) * prfri * (prf / prfri)**0.11 /
              (1+12.7*(cf/2)**0.5 * (prf**(2/3) - 1)))

        hint = kf * nug / dri

        #  Ec. 3.50 Barbero
        qcrit = sigma * eext * (tfe**4 - text**4) + hext * (tfe - text)

        #Ec. 3.51 Barbero
        ucrit = 4 * sigma * eext * tfe**3 + hext
        #krec = (0.0153)*(trec) + 14.77 # trec ya está en ºC
        # Ec. 3.22
        urec = 1/(
            (1/hint) +
            (dro*np.log(dro/dri))/(2*krec)
            )

        #Ec. 3.20 Barbero
        qabs = (hce.get_pr_opt() *
                cg * dni *
                hce.get_pr_shadows() *
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

        if qabs > 0:

            while (errtro > 0.0001 or errpr > 0.000001):

                step += 1

                # Ec. 4.14
                cp = hot_fluid.get_cp(tf, pressure)
                mu = hot_fluid.get_dynamic_viscosity(tf, pressure)
                rho = hot_fluid.get_density(tf, pressure)
                kf = hot_fluid.get_thermal_conductivity(tf, pressure)
                #  alpha : difusividad térmica
                alpha = kf / (rho * cp)
                prf = mu / alpha  #  Prandtl = viscosidad dinámica / difusividad_termica
                redri = 4 * massflow / (mu * np.pi * dri)  # Reynolds
                #  Prandtl number at temperature tri
                kfpri = hot_fluid.get_thermal_conductivity(tri, pressure)
                rhori =  hot_fluid.get_density(tri, pressure)
                cpri = hot_fluid.get_cp(tri, pressure)
                alphari = kfpri / (rhori * cpri)
                muri = hot_fluid.get_dynamic_viscosity(tri, pressure)
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
                                          maxiter=10000)

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
                Model.set_tout(hce, qabs, pr1, cp)
                tf = hce.tout
                tri = tf

        else:
            pr = 0

        Model.set_tout(hce, qabs, pr1, cp)
        hce.pr = pr1

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
    def simulateHCE(cls, hce, hot_fluid, dni):

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
    def simulateHCE(cls, hce, hot_fluid, dni):

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
        self.pr = 0


    def get_previous(self):

        return self.sca.hces[self.hce_order-1]

    def get_index(self):
        return ([self.sca.loop.solarfield.name,
                self.sca.loop.loop_order,
                self.sca.sca_order,
                self.hce_order])


    def get_pr_opt_peak(self):

        alpha = self.get_absorptivity()
        tau = self.get_transmissivity()
        rho = self.get_reflectance()
        gamma = self.get_solar_fraction()

        return alpha * tau * rho * gamma


    def get_pr_opt(self):
        return 1.0


    def get_pr_geo(self):
        return 1.0


    def get_absorptivity(self):
        return 1.0


    def get_transmissivity(self):
        return 1.0


    def get_reflectance(self):
        return self.sca.parameters['Reflectance']


    def get_solar_fraction(self):
        return 1.0

    def get_pr_shadows(self):
        return 1.0

#    def get_IAM(self, theta):
#
#        theta = np.radians(theta)
#        F0 = self.sca.parameters['F0']
#        F1 = self.sca.parameters['F1']
#        F2 = self.sca.parameters['F2']
#        return F0+(F1*theta+F2*theta**2)/np.cos(theta)


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

    def get_tin(self):

        return self.scas[0].hces[0].tin


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


    def estimateMassFlow(cls, hce, hot_fluid, dni):
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

    def __init__(self, plant, solarfield_settings, loop_settings):

        self.plant = plant
        self.name = solarfield_settings['name']
        self.loops = []
        self.massflow =loop_settings['massflow']
        self.tin= 0.0
        self.tout = 0.0

  

    def get_massflow(self):

        totalmassflow = 0.0
        for l in self.loops:
            totalmassflow += l.massflow

        return totalmassflow


    def get_tout(self, hotfluid):
        '''
        Calculates HTF output temperature throughout the solar field as a
        weighted average based on the enthalpy of the mass flow in each
        loop of the solar field
        '''
        H = 0.0

        pressure = 2000000

        for l in self.loops:
            H += (hotfluid.get_cp(l.scas[-1].hces[-1].tout, pressure) *
                  (l.scas[-1].hces[-1].tout-l.get_tin())*l.massflow)

        self.tout = (self.tin + H /
                     (hotfluid.get_cp(l.scas[-1].hces[-1].tout, pressure) *
                      self.get_massflow()))

        return self.tout
    
    
    def set_tout(self, hotfluid):
        '''
        Calculates HTF output temperature throughout the solar field as a
        weighted average based on the enthalpy of the mass flow in each
        loop of the solar field
        '''
        H = 0.0

        pressure = 2000000

        for l in self.loops:
            H += (hotfluid.get_cp(l.scas[-1].hces[-1].tout, pressure) *
                  (l.scas[-1].hces[-1].tout-l.tin)*l.massflow)

        self.tout = (self.tin + H /
                     (hotfluid.get_cp(l.scas[-1].hces[-1].tout, pressure) *
                      self.massflow))

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

    def __init__(self, plant_settings, sca_settings, hce_settings, hce_model_settings):

        self.name = plant_settings['name']
        self.solarfields = []
        self.ratedtin = plant_settings['ratedtin']
        self.ratedtout = plant_settings['ratedtout']
        self.ratedpressure = plant_settings['ratedpressure']
        self.min_massflow = 0.0
        self.tout = 0.0
        self.tin = self.ratedtin
        self.massflow = 0.0

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

        proto_settings = {'name': 'prototype',
                          'massflow': 0}

        self.protosolarfield = SolarField(self, proto_settings, plant_settings['loop'])

        for l in range(self.total_loops):
                self.protosolarfield.loops.append(
                    Loop(self.protosolarfield, l))
                for s in range(plant_settings['loop']['scas']):
                    self.protosolarfield.loops[-1].scas.append(
                        SCA(self.protosolarfield.loops[-1],
                            s,
                            sca_settings))
                    for h in range (plant_settings['loop']['hces']):
                        self.protosolarfield.loops[-1].scas[-1].hces.append(
                            HCE(self.protosolarfield.loops[-1].scas[-1],
                            h,
                            hce_settings,
                            hce_model_settings))


    def set_tout(self, hotfluid):
        '''
        Calculates HTF output temperature throughout the solar plant as a
        weighted average based on the enthalpy of the mass flow in each
        loop of the solar field
        '''
        H = 0.0

        pressure = 2000000

        for sf in self.solarfields:
                 
            H += (hotfluid.get_cp(sf.tin, pressure) *
                  (sf.tout - sf.tin)*sf.massflow)

        print("Tin planta", self.tin)
        self.tout = (self.tin + (H /
                     (hotfluid.get_cp(self.tin, pressure) * self.massflow)))
        
    def set_tin(self, hotfluid):
        '''
        Calculates HTF output temperature throughout the solar plant as a
        weighted average based on the enthalpy of the mass flow in each
        loop of the solar field
        '''
        H = 0.0

        pressure = 2000000

        for sf in self.solarfields:
                 
            H += (hotfluid.get_cp(sf.tin, pressure) *
                  (sf.tin - 25.0)*sf.massflow)

        self.tin = (25.0 + (H /
                     (hotfluid.get_cp(25.0, pressure) * self.massflow)))
    
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
            for sf in self.solarfields:
                print(sf.name+'.tin')
                sf.tin = datasource[sf.name+'.tin'][0]
                sf.massflow = datasource[sf.name+'.mf'][0]
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
    Power Plant as a set composed by a SolarPlant, a HeatExchanger, a PowerCycle
    and a BOPSystem

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

    def get_demanded_thermalpower(self, netpower):

        return netpower / (self.solarfield_to_exchanger_pr *
                                  self.exchanger_pr *
                                  self.exchanger_to_turbogroup_pr *
                                  self.steam_cycle_pr *
                                  self.turbogenerator_pr)

    def get_demanded_rated_massflow(self, ratedpower, solarplant, hotfluid, coldfluid = None):

        cp_avg = 0.5 * (hotfluid.get_cp(solarplant.ratedtin, solarplant.ratedpressure) +
                        hotfluid.get_cp(solarplant.ratedtout, solarplant.ratedpressure))

        drm = (self.get_demanded_thermalpower(self.ratedpower) /
               (cp_avg * (solarplant.ratedtout - solarplant.ratedtin)))

        return drm


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

    def precalc(self, powersystem, solarplant,
                hotfluid, simulation, hce_settings):

        loop_min_massflow = hotfluid.get_massflow_from_Reynolds(
                hce_settings['dri'],
                solarplant.ratedtin,
                solarplant.ratedpressure,
                hce_settings['min_reynolds'])


        solarplant.min_massflow = solarplant.total_loops * loop_min_massflow
        solarplant.rated_massflow = powersystem.get_demanded_rated_massflow(
                powersystem.ratedpower, solarplant, hotfluid)

        v = 4 * loop_min_massflow / ( np.pi * hce_settings['dri']**2 *
                hotfluid.get_density(solarplant.ratedtin,
                                     solarplant.ratedpressure))

        if solarplant.min_massflow > solarplant.rated_massflow:
            print("Too low massflow", solarplant.min_massflow ,">",
                  solarplant.rated_massflow)
        else:
            print("Rated massflow = ", solarplant.rated_massflow, ">",
                  "Min massflow=", solarplant.min_massflow)

    def runSolarPlant(self, model, solarplant, site, data_source, hot_fluid):

        if self.type == "type0":
            self.simulateSolarPlant(model, solarplant, site, data_source, hot_fluid)
        elif self.type == "type1":
            self.benchmarkSolarPlant(model, solarplant, site, data_source, hot_fluid)
        else:
            return None


    def simulateSolarPlant(self, model, solarplant, site, weather, hot_fluid):

        if self.type != 'type0':
            return None

        solarplant.initializePlant(self)

        for row in weather.weatherdata[0].iterrows():
            solarpos = pvlib.solarposition.get_solarposition(row[0],
                                                        site.latitude,
                                                        site.longitude,
                                                        site.altitude,
                                                        pressure = row[1]['Pressure'],
                                                        temperature=row[1]['DryBulb'])

            aoi = float(pvlib.irradiance.aoi(0, 0, solarpos['zenith'][0],
                                             solarpos['azimuth'][0]))


            for sf in solarplant.solarfields:
                for l in sf.loops:
                    for s in l.scas:
                        for h in s.hces:
                            model.simulateHCE(h, hot_fluid, row[1]['DNI'],
                                              row[1]['Wspd'], row[1]['DryBulb'])
                sf.set_tout(hot_fluid)
            solarplant.set_tout(hot_fluid)
 



    def benchmarkSolarPlant(self, model, solarplant, site, fielddata, hot_fluid):

        if self.type != 'type1':
            return None

        solarplant.initializePlant(self, fielddata.fielddata)
        
        results = pd.DataFrame(fielddata.fielddata.index)
        
        #print(results)    

        for row in fielddata.fielddata.iterrows():            
            solarposition = pvlib.solarposition.get_solarposition(row[0],
                                                        site.latitude,
                                                        site.longitude,
                                                        site.altitude,
                                                        pressure = row[1]['Pressure'],
                                                        temperature=row[1]['DryBulb'])
            aoi = float(pvlib.irradiance.aoi(0, 0, solarposition['zenith'][0],
                                             solarposition['azimuth'][0]))
            for sf in solarplant.solarfields:
                sf.massflow = row[1][sf.name+'.mf']
                sf.tin = row[1][sf.name+'.tin']
                for l in sf.loops:
                    l.massflow = sf.massflow / len(sf.loops)
                    for s in l.scas:
                        for h in s.hces:
                            model.simulateHCE(h, hot_fluid, row[1]['DNI'],
                                              row[1]['Wspd'],
                                              row[1]['DryBulb'])
                sf.set_tout(hot_fluid)
                print(sf.name, "Tout:",sf.tout, "massflow", sf.massflow)
            
            solarplant.set_massflow()
            solarplant.set_tin(hot_fluid)
            solarplant.set_tout(hot_fluid)
            # print(row.index.strftime('%y/%m/%d %H:%M'),
            #      'PRtotal:', format( sf.loops[-1].scas[-1].hces[-1].get_pr_total(row.index, site, site), '.2f'))
            print(row[0].strftime('%y/%m/%d %H:%M'), 'DNI: ', row[1]['DNI'], 
                  "MF:", solarplant.massflow, 'Real:',row[1][sf.name+'.tin'],
                  '->',row[1][sf.name+'.tout'], 'Calc:', sf.tin, '->', 
                  round(solarplant.tout))

        #TO-DO
        estimated_pr = 0.0
        actual_pr = 0.0
        estimated_tout = 0.0
        actual_tout = 0.0
        rejected_solar_energy = 0.0






class HeatExchanger(object):

    def __init__(self, hotfluid, coldfluid, pr_heatexchanger,
                 hotfluidmassflow, coldfluidmassflow):

        self.pr_heatexchanger = pr_heatexchanger

    def set_fluid_tout(self):
        pass

    def hot_fluid_tout():
        return


class HeatStorage(object):

    def __init__(self, settings):

        self.name = pr_heatexchanger

    def set_fluid_tout(self):
        pass

    def hot_fluid_tout():
        return


class ThermodynamicCycle (object):

    def __init__(self, settings):

        self.name = settings['name']

    @classmethod
    def get_NCA_pr(cls, tf, text):

        return (1-(text/tf)**(1/2))

    @classmethod
    def get_Carnot_pr(cls, tf, text):

        return (1-(text/tf))

class Fluid_CoolProp(object):

    def __init__(self, settings = None):

        self.name = settings['name']
        self.coolpropID = settings['CoolPropID']
        print(self.coolpropID)

    def get_density(self, t, p):

        t += 273.15
        return PropsSI('D','T',t,'P', p, self.coolpropID)

    def get_dynamic_viscosity(self, t, p):

        t += 273.15
        return  PropsSI('V','T',t,'P', p, self.coolpropID)

    def get_cp(self, t, p):

        t += 273.15
        return PropsSI('C','T',t,'P', p, self.coolpropID)

    def get_thermal_conductivity(self, t, p):
        ''' Saturated Fluid conductivity at temperature t '''
        t += 273.15
        return PropsSI('L','T',t,'P', p, self.coolpropID)

    def get_deltaH(self, t, p):

        t += 273.15
        return PropsSI('H','T',t,'P', p, self.coolpropID)

    def get_Reynolds(self, dri, t, p, massflow):

        t += 273.15
        return massflow * np.pi * (dri**3) /( 4 * self.get_density(t, p))

    def get_massflow_from_Reynolds(self, dri, t, p, re):

        t += 273.15
        return re * 4 * self.get_density(t, p) / (np.pi * (dri**3))

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

    def get_deltaH(self, tin, tout, p):

        pass

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


#class HotFluid(Fluid):
#
#     def __init__(self, settings):
#
#         super().__init__(settings)
#
#
#class ColdFluid(Fluid):
#
#     def __init__(self, settings):
#
#         if self.settings['name'] in ['Water']:
#
#
#
#          super().__init__(settings)


class Weather(object):

    def __init__(self, settings):
        self.filename = settings['filename']
        self.filepath = settings['filepath']
        self.file = self.filepath + self.filename
        self.openWeatherDataFile(self.file)

    def openWeatherDataFile(self, path = None):
        '''


        Parameters
        ----------
        path : string, optional

        DESCRIPTION. URI of the weatherfile. The default is None.

        Tuple of the form (data, metadata).
        data : DataFrame
             A pandas dataframe with the columns described in the table
             below. For more detailed descriptions of each component, please
             consult the TMY3 User's Manual ([1]), especially tables 1-1
             through 1-6.
        metadata : dict
            The site metadata available in the file.
        Notes
        -----
        The returned structures have the following fields.
        ===============   ======  ===================
        key               format  description
        ===============   ======  ===================
        altitude          Float   site elevation
        latitude          Float   site latitudeitude
        longitude         Float   site longitudeitude
        Name              String  site name
        State             String  state
        TZ                Float   UTC offset
        USAF              Int     USAF identifier
        ===============   ======  ===================
        =============================       ======================================================================================================================================================
        TMYData field                       description
        =============================       ======================================================================================================================================================
        TMYData.Index                       A pandas datetime index. NOTE, the index is currently timezone unaware, and times are set to local standard time (daylight savings is not included)
        TMYData.ETR                         Extraterrestrial horizontal radiation recv'd during 60 minutes prior to timestamp, Wh/m^2
        TMYData.ETRN                        Extraterrestrial normal radiation recv'd during 60 minutes prior to timestamp, Wh/m^2
        TMYData.GHI                         Direct and diffuse horizontal radiation recv'd during 60 minutes prior to timestamp, Wh/m^2
        TMYData.GHISource                   See [1]_, Table 1-4
        TMYData.GHIUncertainty              Uncertainty based on random and bias error estimates                        see [2]_
        TMYData.DNI                         Amount of direct normal radiation (modeled) recv'd during 60 mintues prior to timestamp, Wh/m^2

        Returns
        -------
        None
        '''

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


    def get_weather_data_site(self):

        return self.weatherdata[1]

    def loadWeatherDataFile(self):

        self.weatherdata = pvlib.iotools.tmy.read_tmy2(self.file)


class FieldData(object):

    def __init__(self, settings):
        self.filename = settings['filename']
        self.filepath = settings['filepath']
        self.file = self.filepath + self.filename
        self.tags = settings['tags']
                
        self.openFieldDataFile(self.file)
        self.rename_columns()


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
                        
                        self.fielddata = pd.read_csv(path, sep=';',
                                                     decimal= ',',
                                                     dtype= float,
                                                     parse_dates=['datetime'], 
                                                     date_parser=dateparse,
                                                     index_col=0)
                        self.file = path
                    elif strext == ".xls":
                        print("xls...")
                        self.fielddata = pd.read_excel(path)
                        self.file = path
                    else:
                        print("unknow extension ", strext)
                        return
            else:
                strfilename, strext = os.path.splitext(path)

                if  strext == ".csv":
                    print("csv...")
                    self.fielddata = pd.read_csv(path, sep=';',
                                                 decimal= ',',
                                                 dayfirst=True,
                                                 index_col=0)
                    self.file = path
                elif strext == ".xls":
                    print("xls...")
                    self.fielddata = pd.read_excel(path)
                    self.file = path
                else:
                    print("unknow extension ", strext)
                    return

        except Exception:
            raise
            txMessageBox.showerror('Error loading FieldData File',
                                   'Unable to open file: %r', self.file)
            
        self.fielddata.index = pd.to_datetime(self.fielddata.index)
        

    def rename_columns(self):
                 
        rename_dict = dict(zip(self.tags.values(), self.tags.keys()))
        self.fielddata.rename(columns = rename_dict, inplace = True)
        columns_to_drop = []
        for c in  self.fielddata.columns:
            if c not in self.tags.keys():
                columns_to_drop.append(c)
        self.fielddata.drop(columns = columns_to_drop, inplace = True)


        

# format = '%Y-%m-%d %H:%M:%S'
# df['Datetime'] = pd.to_datetime(df['date'] + ' ' + df['time'], format=format)
# df = df.set_index(pd.DatetimeIndex(df['Datetime']))



