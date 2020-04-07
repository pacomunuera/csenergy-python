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
import copy
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


class Model(object):

    def __init__(self, settings):
        self.model_settings = settings


class ModelBarbero4grade(Model):


    def __ini__(self, settings):
        super(Model, self).__init__(settings)


    @classmethod
    def calc_pr(cls, hce, hotfluid, aoi, solarpos, row):

        flag_0 = datetime.now()

        pressure = hce.sca.loop.pin
        hce.set_tin()
        hce.set_pin()
        hce.tout = hce.tin
        tin = hce.tin
        tf = hce.tin
        tro = hce.tin
        tri = hce.tin

        massflow = hce.sca.loop.massflow

        dni = row[1]['DNI']
        wspd = row[1]['Wspd']
        text = row[1]['DryBulb']

        hce.set_krec(hce.tin)
        krec = hce.parameters['krec']
        sigma = sc.constants.sigma
        dro = hce.parameters['dro']
        dri = hce.parameters['dri']
        dgo = hce.parameters['dgo']
        dgi = hce.parameters['dgi']
        L = hce.parameters['long']
        A = hce.sca.parameters['Aperture']
        x = 1 #  Calculation grid fits hce longitude
        IAM = hce.sca.get_IAM(aoi)
        pr_opt_peak = hce.get_pr_opt_peak(aoi, solarpos, row)
        pr_geo = hce.get_pr_geo(aoi, solarpos, row)
        pr_shadows = hce.get_pr_shadows(aoi, solarpos, row)
        pr_shadows = 0.5

        cg = A /(np.pi*dro)

        #  nu_air = cinematic viscosity PROVISIONAL A FALTA DE VALIDAR TABLA
        # tabla en K
        nu_air = Air.get_cinematic_viscosity(text)
        wspd = row[1]['Wspd']

        #  Reynols number for wind at
        reext = dgo * wspd / nu_air
        hext, eext = cls.get_hext_eext(hce, reext, tro, wspd)

        cp = hotfluid.get_cp(tf, hce.pin)
        cpri = cp

        # Ec. 4.14
        mu = hotfluid.get_dynamic_viscosity(tf, pressure)
        muri = mu

        rho = hotfluid.get_density(tf, pressure)
        rhori = rho

        kf = hotfluid.get_thermal_conductivity(tf, pressure)
        kfpri = kf

        alpha = kf / (rho * cp)
        alphari = alpha

        prf = mu / alpha  #  Prandtl = viscosidad dinámica / difusividad_termica
        prfri = prf

        #  Correlación de Gnielinski (utilizado por Forristall [1]
        #  y verificado en [38]) según ec. 4.15
        redri = 4 * massflow / (mu * np.pi * dri)  # Reynolds
        cf = (1.58 * np.log(redri) - 3.28)**-2

        # Número de Nusselt alto, casi independiente inicialmente de tri
        # ya que estamos considerando tri = tf o casi igual

        nug = ((cf / 2) * (redri - 1000) * prf * (prf / prfri)**0.11 /
               (1 + 12.7 * (cf / 2)**0.5 * (prf**(2 / 3) - 1)))


        #  Internal transmission coefficient.
        hint = kf * nug / dri

        # Ec. 3.23
        qperd = sigma * eext * (tro**4 - text**4) + hext * (tro - text)

        # Ec. 3.50 Barbero
        qcrit = qperd

        # Ec. 3.51 Barbero
        ucrit = 4 * sigma * eext * tin**3 + hext
        # Ec. 3.22
        urec = 1 / ((1 / hint) + (dro * np.log(dro / dri)) / (2 * krec))
        NTU = urec * x * L * sc.pi * dro / (massflow * cp)  # Ec. 3.30

        # Ec. 3.20 Barbero
        qabs = (pr_opt_peak * IAM * cg * dni * pr_geo * pr_shadows)

        if qabs > 0:
            pr0 = (1 - qcrit / qabs)/(1 + ucrit / urec)

            if pr0 <0:
                pr0 = 0.0
        else:
            pr0 = 0.0


        # if hce.tin >= hotfluid.tmax or hce.sca.status == 'defocused':
        #     hce.sca.status = 'defocused'
        #     qabs = 0.0  #defocused
        #     pr0 = 0.0
        # else:
        #     hce.sca.status == 'focused'
        #     if qabs > 0.0:
        #         pr0 = (1 - qcrit / qabs)/(1 + ucrit / urec)
        #     else:
        #         pr0 = 0.0

        pr1 = pr0
        tro1 = tf + qabs * pr1 / urec

        if qabs > qcrit:
            errtro = 1.
            errpr = 1.
            step = 0

            while (errtro > 0.1 or errpr > 0.01):

                step += 1
                flag_1 = datetime.now()

                f0 = qabs/(urec*(tin-text))
                f1 = ((4*sigma*eext*text**3)+hext)/urec
                f2 = 6*(text**2)*(sigma*eext/urec)*(qabs/urec)
                f3 = 4*text*(sigma*eext/urec)*((qabs/urec)**2)
                f4 = (sigma*eext/urec)*((qabs/urec)**3)

                fx = lambda pr0: (1 - pr0 -
                                  f1 * (pr0 + (1 / f0)) -
                                  f2 * ((pr0 + (1 / f0))**2) -
                                  f3 * ((pr0 + (1 / f0))**3) -
                                  f4 * ((pr0 + (1 / f0))**4))

                dfx = lambda pr0: (-1 - f1 -
                                   2 * f2 * (pr0 + (1 / f0)) -
                                   3 * f3 * (pr0 + (1 / f0))**2 -
                                   4 * f4 * (pr0 + (1 / f0))**3)

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

                errpr = abs(pr2-pr1)
                pr1 = pr2
                hce.pr = pr1
                hce.qabs = qabs
                hce.set_tout(hotfluid)
                hce.set_pout(hotfluid)
                tf = 0.5 * (hce.tin + hce.tout)

                tri = tro1


                cp = hotfluid.get_cp(tf, pressure)
                cpri = hotfluid.get_cp(tri, pressure)

                mu = hotfluid.get_dynamic_viscosity(tf, pressure)
                muri = hotfluid.get_dynamic_viscosity(tri, pressure)

                rho = hotfluid.get_density(tf, pressure)
                rhori =  hotfluid.get_density(tri, pressure)


                kf = hotfluid.get_thermal_conductivity(tf, pressure)
                kfpri = hotfluid.get_thermal_conductivity(tri, pressure)

                #  alpha : difusividad térmica
                alpha = kf / (rho * cp)
                alphari = kfpri / (rhori * cpri)

                prf = mu / alpha  #  Prandtl = viscosidad dinámica / difusividad_termica
                prfri =  muri / alphari

                redri = 4 * massflow / (mu * np.pi * dri)  # Reynolds

                # Número de Nusselt alto, casi independiente inicialmente de tri
                # ya que estamos considerando tri = tf o casi igual
                nug = ((cf / 2) * (redri - 1000) * prf * (prf / prfri)**0.11 /
                       (1 + 12.7 * (cf / 2)**0.5 * (prf**(2 / 3) - 1)))

                #nudb = 0.023 * redri**0.8 * prf**0.4
                hint = kf * nug / dri
                hext, eext = cls.get_hext_eext(hce, reext, tro1, wspd)

                trec = (tro1+tri)/2
                # Ec. 4.22
                hce.set_krec(trec)
                krec = hce.parameters['krec']
                # Ec. 3.22
                urec = 1 / ((1 / hint) + (dro * np.log(dro / dri)) / (2 * krec))
                NTU = urec * x * L * sc.pi * dro / (massflow * cp)  # Ec. 3.30

                tro2 = tf+qabs*pr2/urec
                errtro = abs(tro2-tro1)
                tro1 = tro2
                qperd = sigma * eext * (tro1**4 - text**4) + hext * (tro1 - text)

                flag_2 = datetime.now()

                # if (hce.sca.loop.loop_order == 0 and
                #     hce.sca.sca_order == (len(hce.sca.loop.scas) - 1) and
                #     hce.hce_order == (len(hce.sca.hces)-1)):
                #     delta_2 = flag_2 - flag_0
                #     #print("hce_while:", hce.get_index(), delta_2.total_seconds())

            # if (hce.sca.loop.loop_order == 0 and
            #     hce.sca.sca_order == (len(hce.sca.loop.scas) - 1) and
            #     hce.hce_order == (len(hce.sca.hces)-1)):
            #     delta_1 = flag_1 - flag_0
                #print("hce_offwhile:", hce.get_index(), delta_1.total_seconds())
        else:
            hce.pr = pr1
            hce.qperd = qperd
            hce.qabs = qabs

            #HL = hotfluid.get_deltaH(hce.tin, hce.sca.loop.pin)
            #h = qperd / hce.sca.loop.massflow
            #hce.tout = hotfluid.get_T(HL - h, hce.sca.loop.pin)
            hce.set_tout(hotfluid)
            #FLUJO LAMINAR
            hce.set_pout(hotfluid)
            flag_3 = datetime.now()



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
    def calc_pr(cls, hce, hot_fluid, dni):

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

        tin = hce.tin
        tf = tin
        hce.tout = tin

        text = 22.0

        #Ec. 3.20 Barbero
        qabs = (hce.parameters['pr_opt']*
                hce.parameters['cg'] * dni *
                hce.parameters['pr_shw'] *
                hce.parameters['pr_geo'])

        # Ec. 3.50 Barbero
        qcrit = sigma*eext*(tin**4-text**4)+hext*(tin-text)

        #Ec. 3.51 Barbero
        ucrit = 4*sigma*eext*tin**3+hext
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
    def calc_pr(cls, hce, hot_fluid, dni):

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
        tin = hce.tin
        tf = tin
        hce.tout = tin

        text = 22.0

        # Ec. 3.20 Barbero
        qabs = (hce.parameters['pr_opt'] *
                cg * dni *
                hce.parameters['pr_shw'] *
                hce.parameters['pr_geo'])

        # Ec. 3.50 Barbero
        qcrit = sigma*eext*(tin**4-text**4)+hext*(tin-text)

        # Ec. 3.51 Barbero
        ucrit = 4*sigma*eext*tin**3+hext
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
        self.parameters = dict(settings, **model_settings)
        self.tin = 0.0
        self.tout = 0.0
        self.pr = 0.0
        self.qabs = 0.0
        self.qperd = 0.0
        self.pout = 0.0
        self.pin = 0.0

    def set_tin(self):

        if self.hce_order > 0:
            self.tin = self.sca.hces[self.hce_order-1].tout
        elif self.sca.sca_order > 0:
            self.tin = self.sca.loop.scas[self.sca.sca_order-1].hces[-1].tout
        else:
            self.tin = self.sca.loop.tin

    def set_pin(self):

        if self.hce_order > 0:
            self.pin = self.sca.hces[self.hce_order-1].pout
        elif self.sca.sca_order > 0:
            self.pin = self.sca.loop.scas[self.sca.sca_order-1].hces[-1].pout
        else:
            self.pin = self.sca.loop.pin

    def set_tout(self, hotfluid):

        HL = hotfluid.get_deltaH(self.tin, self.sca.loop.pin)

        if self.qabs > 0:

            h = (np.pi * self.parameters['dro'] * self.parameters['long'] *
                 self.qabs * self.pr / self.sca.loop.massflow)

        else:
            h = -self.qperd / self.sca.loop.massflow

        self.tout = hotfluid.get_T(HL + h, self.sca.loop.pin)


    def set_pout(self, hotfluid):

        # TO-DO CÁLCULLO PERDIDA DE CARGA:
        # Ec. Colebrook-White para el cálculo del factor de fricción de Darcy

        re_turbulent = 4000

        k = self.parameters['Inner surface roughness']
        D = self.parameters['Absorber tube inner diameter']
        re = hotfluid.get_Reynolds(D, self.tin, self.pin, self.sca.loop.massflow)


        if re < re_turbulent:

            print("laminar", re)

            darcy_factor = 64 / re

        else:


            # a = (k /  D) / 3.7
            a = k / 3.7
            b = 2.51 / re
            x = 1

            fx = lambda x: x + 2 * np.log10(a + b * x )

            dfx = lambda x: 1 + (2 * b) / (np.log(10) * (a + b * x))

            root = sc.optimize.newton(fx,
                                      x,
                                      fprime=dfx,
                                      maxiter=10000)

            darcy_factor = 1 / (root**2)


        rho = hotfluid.get_density(self.tin, self.pin)
        v = 4 * self.sca.loop.massflow / (rho * np.pi * D**2)
        g = sc.constants.g

        # Ec. Darcy-Weisbach para el cálculo de la pérdida de carga
        deltap_mcl = darcy_factor * (self.parameters['long'] / D ) * (v**2 / (2 * g))

        deltap = deltap_mcl * rho * g

        # if self.hce_order == 35 and self.sca.loop.loop_order==5:

        #     datos = [re, a, b, root, darcy_factor, v, self.sca.loop.massflow,
        #            hotfluid.get_density(self.tin, self.pin),  deltap]

        #     claves = ['re', 'a', 'b',
        #               'root', 'darcy_factor', 'v', 'mf', 'rho', 'deltap']

        #     dict_pd = dict(zip(claves,datos))

        #     datos_deltap = pd.DataFrame(datos)
        #     print(datos_deltap)

        self.pout = self.pin - deltap

        # if self.hce_order == 35 and self.sca.loop.loop_order==5:
        #     print('Pout',self.pin, self.pout)


    def set_krec(self, t):

        self.parameters['krec'] = (0.0153)*(t - 273.15) + 14.77

    def get_previous(self):

        return self.sca.hces[self.hce_order-1]

    def get_index(self):

        if hasattr(self.sca.loop, 'subfield'):
            index = [self.sca.loop.subfield.name,
                self.sca.loop.loop_order,
                self.sca.sca_order,
                self.hce_order]
        else:
            index = ['prototype',
                self.sca.loop.loop_order,
                self.sca.sca_order,
                self.hce_order]

        return index

    def get_pr_opt_peak(self, aoi, solarpos, row):

        alpha = self.get_absorptance()
        tau = self.get_transmissivity()
        rho = self.sca.parameters['Reflectance']
        gamma = self.sca.get_solar_fraction(aoi, solarpos, row)

        pr_opt_peak = alpha * tau * rho * gamma

        if pr_opt_peak > 1 or pr_opt_peak < 0:
            print("ERROR", pr_opt_peak)

        return pr_opt_peak


    def get_pr_geo(self, aoi, solarpos, row):

        # Llamado "bordes" en Tesis. Pérdidas de los HCE de cabecera según aoi
        pr_geo = 1- (self.sca.parameters["Focal Len"] * np.tan(np.radians(aoi)) /
                     (len(self.sca.hces) * self.parameters["long"]))

        if pr_geo > 1:
            pr_geo = 1.0

        if pr_geo > 1 or pr_geo < 0:
            print("ERROR", pr_geo)

        return pr_geo


    def get_pr_shadows(self, aoi, solarpos, row):

        # Llamado "sombras" en Tesis. Pérdidas por sombras. ¿modelar sobre el SCA?

        # Sombras debidas a otros lazos

        shadowing = (1 -
                     np.sin(np.radians(abs(solarpos['elevation'][0]))) *
                     self.sca.loop.row_spacing /
                     self.sca.parameters['Aperture'])

        if shadowing < 0.0:
            shadowing = 0.0

        if solarpos['elevation'][0] < 0:
            shadowing = 1.0

        pr_shadows = 1 - shadowing

        if  pr_shadows > 1 or  pr_shadows < 0:
            print("ERROR",  pr_shadows)

        return pr_shadows


    def get_absorptance(self):

        #
        alpha = self.parameters['absorptance']

        return alpha


    def get_transmissivity(self):

        #dependiendo si hay vídrio y si hay vacío
        tau = self.parameters['transmissivity']

        return tau


    def get_reflectance(self):

        return self.sca.parameters['Reflectance']

    # def get_pr_total(self, dateindex, site, data):

    #     aoi = self.sca.get_aoi(dateindex, site, data)

    #     return (self.pr * self.get_pr_shadows() *
    #             self.get_pr_opt_peak() * self.get_pr_geo() *
    #             np.cos(np.radians(aoi)) * self.sca.get_IAM(aoi))


class SCA(object):

# Uso de __slots__ permite ahorrar memoria RAM
#    __slots__= ['SCA_configuration']
    def __init__(self, loop, sca_order, settings):

        self.loop = loop
        self.sca_order = sca_order
        self.hces = []
        self.status = 'focused'
        self.tracking_angle = 0.0
        self.parameters = dict(settings)
        self.surface_tilt = 0.0
        self.surface_azimuth = 180.0


    def get_solar_fraction(self, aoi, solarpos, row):

        if self.status == 'defocused':
            solarfraction = 0.0
        elif self.status == 'focused':

            # faltaría: factor de forma del Sol, dependiente de solpos, geometría
            # absorbedor, sca y atmósfera.
            # Factor: se puede usar para tasa de espejos rotos, por ejemplo
            # Availability: podría ser un valor binario 0 o 1

            solarfraction = (self.parameters['Reflectance'] *
                             self.parameters['Geom.Accuracy'] *
                             self.parameters['Track Twist'] *
                             self.parameters['Cleanliness'] *
                             self.parameters['Dust'] *
                             self.parameters['Factor'] *
                             self.parameters['Availability'] *
                             self.get_sca_shadows(aoi, solarpos, row))
        else:
            solarfraction = 1.0

        if  solarfraction > 1 or  solarfraction < 0:
            print("ERROR",  solarfraction)

        return solarfraction

    def get_sca_shadows(self, aoi, solarpos, row):

        # Pérdidas por sombras en el sca:
        # - Sombras de otros SCA. Cálculo geometrico a partir de la distancia entre lazos.
        # - Nubes

        pr_shadows = 1.0

        return pr_shadows


    def get_IAM(self, theta):

        F0 = self.parameters['IAM Coefficient F0']
        F1 = self.parameters['IAM Coefficient F1']
        F2 = self.parameters['IAM Coefficient F2']

#        PROVISIONAL HASTA ACLARAR LA FÓRMULA DE IAM
#        kiam = (1-2.23073e-4 * theta - 1.1e-4 * theta**2 +
#                3.1859e-6 * theta**3 - 4.85509e-8 * theta** 4)
#        return kiam


        if (theta > 0 and theta < 80):
            theta = np.radians(theta)
            kiam = F0 + (F1*theta+F2*theta**2) / np.cos(theta)

            if kiam > 1:
                kiam = 2 -kiam

        else:
            kiam = 0.0

        if  kiam > 1 or  kiam < 0:
            print("ERROR",  kiam, theta)

        return kiam


    def get_aoi(self, solarpos):

        sigmabeta = 0.0
        beta0 = 0.0

        if self.parameters['Tracking Type'] == 1: # N-S single axis tracker
            if solarpos['azimuth'][0] > 0 and solarpos['azimuth'][0] <= 180:
                surface_azimuth = 90 # Surface facing east
            else:
                surface_azimuth = 270 # Surface facing west
        elif self.parameters['Tracking Type'] == 2:  # E-W single axis tracker
            surface_azimuth = 180  # Surface facing the equator

        beta0 = np.degrees(
                    np.arctan(np.tan(np.radians(solarpos['zenith'][0])) *
                              np.cos(np.radians(surface_azimuth -
                                                solarpos['azimuth'][0]))))

        if beta0 >= 0:
            sigmabeta = 0
        else:
            sigmabeta = 1

        beta = beta0 + 180 * sigmabeta
        aoi = pvlib.irradiance.aoi(beta,
                                   surface_azimuth,
                                   solarpos['zenith'][0],
                                   solarpos['azimuth'][0])
        return aoi


class Loop(object):

    def __init__(self, solarfield, loop_order, subfield = None):

        self.solarfield = solarfield
        self.subfield = subfield
        self.scas = []
        self.loop_order = loop_order

        self.tin = 0.0
        self.tout = 0.0
        self.pin = 0.0
        self.pout = 0.0
        self.tmax = 0.0
        self.massflow = 0.0
        self.req_massflow = 0.0  # Required massflow (in order to achieve the setpoint tout)
        self.pr_req_massflow = 0.0
        self.pr_act_massflow = 0.0

        self.act_tin = 0.0
        self.act_tout = 0.0
        self.act_pin = 0.0
        self.act_pout = 0.0
        self.act_massflow = 0.0  # Actual massflow measured by the flowmeter

        # self.rated_tin = 0.0
        # self.rated_tout = 0.0
        # self.rated_pin = 0.0
        # self.rated_pout = 0.0
        # self.rated_massflow = 0.0

        self.wasted_power = 0.0
        self.row_spacing = self.solarfield.solarfield_settings['row_spacing']

        self.tracking = True

    def initialize(self):

        self.massflow = self.subfield.massflow / len(self.subfield.loops)
        self.tin = self.subfield.tin
        self.pin = self.subfield.pin


    def initialize_actual(self):

        self.act_massflow = self.subfield.act_massflow / len(self.subfield.loops)
        self.act_tin = self.subfield.act_tin
        self.act_pin = self.subfield.act_pin
        self.act_tout = self.subfield.act_tout
        self.act_pout = self.subfield.act_pout

        self.massflow = self.act_massflow
        self.tin = self.act_tin
        self.pin = self.act_pin
        self.tout = self.act_tout
        self.pout = self.act_pout


    def initialize_rated(self):

        self.rated_massflow = self.subfield.rated_massflow / len(self.subfield.loops)
        self.rated_tin = self.subfield.rated_tin
        self.rated_pin = self.subfield.rated_pin
        self.rated_tout = self.subfield.rated_tout
        self.rated_pout = self.subfield.rated_pout

        self.massflow = self.rated_massflow
        self.tin = self.rated_tin
        self.pin = self.rated_pin
        self.tout = self.rated_tout
        self.pout = self.rated_pout


    def get_loop_avg_pr(self):

        pr_list = []
        for s in self.scas:
            for h in s.hces:
                pr_list.append(h.pr)

        return np.mean(pr_list)


    def load_from_prototype_loop(self, prototype_loop):

        self.massflow = prototype_loop.massflow
        self.tin = prototype_loop.tin
        self.act_tin = prototype_loop.act_tin
        self.pin = prototype_loop.pin
        self.act_pin = prototype_loop.act_pin
        self.tout = prototype_loop.tout
        self.pout = prototype_loop.pout
        self.act_tout = prototype_loop.act_tout
        self.act_pout = prototype_loop.act_pout
        self.pr_req_massflow = prototype_loop.pr_req_massflow
        self.pr_act_massflow = prototype_loop.pr_act_massflow
        self.req_massflow = prototype_loop.req_massflow
        self.act_massflow = prototype_loop.act_massflow


    def calc_loop_pr_for_massflow(self, row, solarpos, hotfluid, model):

        # We work with the minimum massflow at night in order to
        # reduce thermal losses

        for s in self.scas:
            aoi = s.get_aoi(solarpos)
            for h in s.hces:
                model.calc_pr(h, hotfluid, aoi, solarpos, row)

            self.tout = self.scas[-1].hces[-1].tout
            self.pout = self.scas[-1].hces[-1].pout

    def calc_loop_pr_for_tout(self, row, solarpos, hotfluid, model):

        min_massflow = hotfluid.get_massflow_from_Reynolds(
            self.solarfield.hce_settings['dri'],
            self.tin, self.pin,
            self.solarfield.hce_settings['min_reynolds'])

        max_error = 0.1  # % desviation tolerance
        search = True

        while search:

            self.calc_loop_pr_for_massflow(row, solarpos, hotfluid, model)

            err = 100 * (abs(self.tout-self.solarfield.rated_tout) /
                   self.solarfield.rated_tout)

            if err > max_error:

                if self.tout >= self.solarfield.rated_tout:
                    self.massflow *= (1 + err / 100)
                    search = True
                elif (self.massflow > min_massflow):
                    self.massflow *= (1 - err / 100)
                    search = True
                else:
                    self.massflow = min_massflow
                    search = False
            else:
                search = False

        self.req_massflow = self.massflow


    # def calc_loop_pr(self, row, solarpos, hotfluid, model):


    #     if solarpos['zenith'][0] > 90:

    #         # We work with the minimum massflow at night in order to
    #         # reduce thermal losses
    #         self.massflow =  (self.solarfield.recirculation_massflow /
    #                           self.solarfield.total_loops)

    #         for s in self.scas:
    #             aoi = s.get_aoi(solarpos)
    #             for h in s.hces:
    #                 model.calc_pr(h, hotfluid, aoi, solarpos, row)

    #             self.tout = self.scas[-1].hces[-1].tout

    #     else:

    #         min_massflow = hotfluid.get_massflow_from_Reynolds(
    #             self.solarfield.hce_settings['dri'],
    #             self.tin, self.pin,
    #             self.solarfield.hce_settings['min_reynolds'])

    #         max_error = 0.1  # % desviation tolerance
    #         search = True

    #         while search:

    #             for s in self.scas:
    #                 aoi = s.get_aoi(solarpos)
    #                 for h in s.hces:
    #                     model.calc_pr(h, hotfluid, aoi, solarpos, row)

    #             self.tout = self.scas[-1].hces[-1].tout

    #             err = 100 * (abs(self.tout-self.solarfield.rated_tout) /
    #                    self.solarfield.rated_tout)

    #             if err > max_error:

    #                 if self.tout >= self.solarfield.rated_tout:
    #                     self.massflow *= (1 + err / 100)
    #                     search = True
    #                 elif (self.massflow > min_massflow):
    #                     self.massflow *= (1 - err / 100)
    #                     search = True
    #                 else:
    #                     self.massflow = min_massflow
    #                     search = False
    #             else:
    #                 search = False

    #         self.req_massflow = self.massflow


    def show_parameter_vs_x(self, parameters = None):

        data = []

        for s in self.scas:
            for h in s.hces:
                data.append({'num': h.hce_order,
                             'tin': h.tin,
                             'tout': h.tout,
                             'pin': round(h.pin/100000,3),
                             'pout': round(h.pout/100000,3),
                             'pr': h.pr,
                             'qabs': h.qabs,
                             'qperd': h.qperd})

        loop_df = pd.DataFrame(data)

        if parameters:

            loop_df[parameters].plot(
                figsize=(20,10), linewidth=5, fontsize=20)
            plt.xlabel('HCE order', fontsize=20)

        else:

            loop_df[['num', 'tin', 'tout', 'pin', 'pout',
                     'pr', 'qabs', 'qperd']].plot(
                figsize=(20,10), linewidth=5, fontsize=20)

            plt.xlabel('HCE order', fontsize=20)


        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        # pd.set_option('display.max_colwidth', -1)
        print(loop_df[['num', 'tin', 'tout', 'pin', 'pout',
                       'pr', 'qabs', 'qperd']])
        #print(self.datasource.dataframe)


class PrototypeLoop(Loop):

    def __init__(self, solarfield):

        super().__init__(solarfield, 0)

        # self.solarfield = solarfield
        # self.scas = []
        # self.loop_order = loop_order

        # self.tin = 0.0
        # self.tout = 0.0
        # self.pin = 0.0
        # self.pout = 0.0
        # self.tmax = 0.0
        # self.massflow = 0.0
        # self.req_massflow = 0.0  # Required massflow (in order to achieve the setpoint tout)

        # self.act_tin = 0.0
        # self.act_tout = 0.0
        # self.act_pin = 0.0
        # self.act_pout = 0.0
        # self.act_massflow = 0.0  # Actual massflow measured by the flowmeter

        # self.rated_tin = 0.0
        # self.rated_tout = 0.0
        # self.rated_pin = 0.0
        # self.rated_pout = 0.0
        # self.rated_massflow = 0.0

        # self.wasted_power = 0.0
        # self.row_spacing = self.solarfield.solarfield_settings['row_spacing']

        for s in range(self.solarfield.solarfield_settings['loop']['scas']):
            self.scas.append(SCA(self, s, self.solarfield.sca_settings))
            for h in range(self.solarfield.solarfield_settings['loop']['hces']):
                self.scas[-1].hces.append(
                    HCE(self.scas[-1], h, self.solarfield.hce_settings,
                        self.solarfield.hce_model_settings))


    def initialize(self, values = None):

        if values is not None:

            self.massflow = values['massflow']
            self.tin = values['tin']
            self.pin = values['pin']

        else:
            self.massflow = self.solarfield.massflow / self.solarfield.total_loops
            self.tin = self.solarfield.tin
            self.pin = self.solarfield.pin


class Subfield(object):
    '''
    Parabolic Trough Solar Field

    '''

    def __init__(self, solarfield, subfield_settings, loop_settings):

        self.solarfield = solarfield
        self.name = subfield_settings['name']
        self.parameters = subfield_settings
        self.loops = []

        self.tin = 0.0
        self.tout = 0.0
        self.pin = 0.0
        self.pout = 0.0
        self.tmax = 0.0
        self.massflow = 0.0
        self.req_massflow = 0.0

        self.act_tin = 0.0
        self.act_tout = 0.0
        self.act_pin = 0.0
        self.act_pout = 0.0
        self.act_massflow = 0.0

        # self.rated_tin = self.solarfield.solarfield_settings['rated_tin']
        # self.rated_tout = self.solarfield.solarfield_settings['rated_tout']
        # self.rated_pin = self.solarfield.solarfield_settings['rated_pin']
        # self.rated_pout = self.solarfield.solarfield_settings['rated_pout']
        # self.rated_massflow = self.solarfield.solarfield_settings['rated_massflow']

    def calc_massflow(self):

        mf = 0.0
        for l in self.loops:
            mf += l.massflow

        self.massflow = mf


    def calc_req_massflow(self):

        req_mf = 0.0
        for l in self.loops:
            req_mf += l.req_massflow

        self.req_massflow = req_mf

    def calc_pout(self):

        loops_var = []

        for l in self.loops:
            loops_var.append(l.pout * l.massflow)

        self.pout = np.mean(loops_var) / self.massflow

    def calc_tout(self, hotfluid):
        '''
        Calculates HTF output temperature throughout the solar field as a
        weighted average based on the enthalpy of the mass flow in each
        loop of the solar field
        '''
        # H = 0.0
        # H2 = 0.0
        dH = 0.0
        dH_tmax = 0.0

        for l in self.loops:

            dH += l.massflow * hotfluid.get_deltaH(l.tout, l.pout)
            dH_tmax += l.massflow * hotfluid.get_deltaH(l.tmax, l.pout)

        dH /= self.massflow
        dH_tmax /= self.massflow
        self.tout =  hotfluid.get_T(dH, self.pout)
        self.tmax = hotfluid.get_T(dH_tmax, self.pout)


    def set_tin(self):

        self.tin = (self.solarfield.heatexchanger.hotfluid_tout *
                    self.coolPipeLosses
                    )

    def apply_temp_limitation(self, hotfluid):

        for l in self.loops:
            if l.tout > l.solarfield.tmax:
                HH = hotfluid.get_deltaH(l.tout, l.pout)
                HL = hotfluid.get_deltaH(l.solarfield.tmax, l.pout)
                l.tmax = l.solarfield.tmax
                l.wasted_power = l.massflow * (HH - HL)

            else:
                l.tmax = l.tout
                l.wasted_power = 0.0


    def loops_avg_out(self):

        tavg = 0.0
        cont = 0

        for l in self.loops:
            cont += 1
            tavg += l.tout

        return tavg / cont

    def load_rated(self):

        self.massflow = (self.solarfield.rated_massflow * len(self.loops) /
                         self.solarfield.total_loops)
        self.tin = self.solarfield.rated_tin
        self.pin = self.solarfield.rated_pin


    def load_actual(self, row):


        self.act_massflow = row[1][self.name+'.mf']
        self.act_tin = row[1][self.name+'.tin']
        self.act_pin = row[1][self.name+'.pin']
        self.act_tout = row[1][self.name+'.tout']
        self.act_pout = row[1][self.name+'.pout']


        self.massflow = self.act_massflow
        self.tin = self.act_tin
        self.pin = self.act_pin


    # def calcRequired_massflow(self):
    #     req_massflow = 0
    #     for l in self.loops:
    #         req_massflow += l.required_massflow
    #     self.requiered_massflow = req_massflow/self.loops.len()


class SolarField(object):
    '''
    Parabolic Trough Solar Field

    '''

    def __init__(self, solarfield_settings, sca_settings, hce_settings,
                 hce_model_settings):


        self.solarfield_settings = solarfield_settings
        self.hce_settings = hce_settings
        self.hce_model_settings = hce_model_settings
        self.sca_settings = sca_settings
        self.name = solarfield_settings['name']

        self.subfields = []
        self.total_loops = 0

        self.__build__()

        self.tin = 0.0
        self.tout = 0.0
        self.pin = 0.0
        self.pout = 0.0
        self.massflow = 0.0
        self.req_massflow = 0.0

        self.act_tin = 0.0
        self.act_tout = 0.0
        self.act_pin = 0.0
        self.act_pout = 0.0
        self.act_massflow = 0.0

        self.rated_tin = solarfield_settings['rated_tin']
        self.rated_tout = solarfield_settings['rated_tout']
        self.rated_pin = solarfield_settings['rated_pin']
        self.rated_pout = solarfield_settings['rated_pout']
        self.rated_massflow = (solarfield_settings['rated_massflow'] *
                               self.total_loops)
        self.recirculation_massflow = (
            solarfield_settings['recirculation_massflow'] * self.total_loops)

        self.tmax = solarfield_settings['tmax']
        self.tmin = solarfield_settings['tmin']

        # FUTURE WORK
        self.storage_available = False
        self.operation_mode = "subfield_heating"


    def __build__(self):

        for sf in self.solarfield_settings['subfields']:
            self.subfields.append(Subfield(self, sf, self.solarfield_settings))
            for l in range(sf['loops']):
                self.subfields[-1].loops.append(
                    Loop(self, l, self.subfields[-1]))
                for s in range(self.solarfield_settings['loop']['scas']):
                    self.subfields[-1].loops[-1].scas.append(
                        SCA(self.subfields[-1].loops[-1],
                            s,
                            self.sca_settings))
                    for h in range (self.solarfield_settings['loop']['hces']):
                        self.subfields[-1].loops[-1].scas[-1].hces.append(
                            HCE(self.subfields[-1].loops[-1].scas[-1],
                            h,
                            self.hce_settings,
                            self.hce_model_settings))

        for sf in self.subfields:
            self.total_loops += len(sf.loops)


    def load_rated(self, hotfluid):

        for sf in self.subfields:
                sf.load_rated()

    def initialize_rated(self):

        self.massflow = self. rated_massflow
        self.tin = self.rated_tin
        self.tout = self.rated_tout
        self.pin = self.rated_pin
        self.pout = self.rated_pout



    def load_actual(self, hotfluid, row):

        massflow = 0.0
        H_tin = 0.0
        H_tout = 0.0
        list_pin = []
        list_pout = []

        for sf in self.subfields:
            sf.load_actual(row)
            massflow += sf.act_massflow
            H_tin += (sf.act_massflow *
                      hotfluid.get_deltaH(sf.act_tin, sf.act_pin))
            H_tout += (sf.act_massflow *
                       hotfluid.get_deltaH(sf.act_tout, sf.act_pout))
            list_pin.append(sf.act_pin * sf.act_massflow)
            list_pout.append(sf.act_pout * sf.act_massflow)

        H_tin /= massflow
        H_tout /= massflow

        self.act_massflow = massflow
        self.act_pin = np.mean(list_pin) / massflow
        self.act_pout = np.mean(list_pout) / massflow
        self.act_tin = hotfluid.get_T(H_tin, self.act_pin)
        self.act_tout = hotfluid.get_T(H_tout, self.act_pout)


    def initialize_actual(self):

        self.massflow = self.act_massflow
        self.tin = self.act_tin
        self.tout = self.act_tout
        self.pin = self.act_pin
        self.pout = self.act_pout


    def calc_massflow(self):

        mf = 0.0
        req_mf = 0.0

        for sf in self.subfields:
            mf += sf.massflow
            req_mf += sf.req_massflow

        self.massflow = mf
        self.req_massflow = req_mf


    def calc_req_massflow(self):

        mf = 0.0

        for sf in self.subfields:
            mf += sf.req_massflow

        self.req_massflow = mf

    def calc_pout(self):

        subfields_var = []

        for sf in self.subfields:
            subfields_var.append(sf.pout * sf.massflow)

        self.pout = np.mean(subfields_var) / self.massflow

    def calc_tout(self, hotfluid):
        '''
        Calculates HTF output temperature throughout the solar plant as a
        weighted average based on the enthalpy of the mass flow in each
        loop of the solar field
        '''

        dH = 0.0
        dH_tmax = 0.0

        for sf in self.subfields:

            dH += sf.massflow * hotfluid.get_deltaH(sf.tout, sf.pout)
            dH_tmax += sf.massflow * hotfluid.get_deltaH(sf.tmax, sf.pout)

        dH /= self.massflow
        dH_tmax /= self.massflow
        self.tout = hotfluid.get_T(dH, self.pout)
        self.tmax = hotfluid.get_T(dH_tmax, self.pout)

    def set_act_tout(self, hotfluid):
        '''
        Calculates HTF output temperature throughout the solar plant as a
        weighted average based on the enthalpy of the mass flow in each
        loop of the solar field
        '''

        dH_actual = 0.0

        for sf in self.subfields:

            dH_actual += sf.act_massflow * hotfluid.get_deltaH(sf.tout, sf.pout)
            dH_actual += sf.act_massflow * hotfluid.get_deltaH(sf.act_tout, sf.act_pout)

        dH_actual /= self.act_massflow
        self.act_tout = hotfluid.get_T(dH_actual, self.act_pout)


    def set_tin(self, hotfluid):
        '''
        Calculates HTF output temperature throughout the solar plant as a
        weighted average based on the enthalpy of the mass flow in each
        loop of the solar field
        '''
        H = 0.0

        for sf in self.subfields:

            H += (hotfluid.get_deltaH(sf.tin, sf.pin) * sf.massflow)

        H /= self.massflow
        self.tin = hotfluid.get_T(H, self.rated_pin)


    def set_pin(self):

        subfields_var = []

        for sf in self.subfields:
            subfields_var.append(sf.pin * sf.massflow)

        self.pin = np.mean(subfields_var) / self.massflow


    def set_act_pin(self):

        subfields_var = []

        for sf in self.subfields:
            subfields_var.append(sf.act_pin * sf.act_massflow)

        self.act_pin = np.mean(subfields_var) / self.act_massflow


    def get_thermalpoweroutput(self, hotfluid):

        HL = hotfluid.get_deltaH(self.tin, self.pin)
        HH = hotfluid.get_deltaH(self.tout, self.pout)

        return (HH - HL) * self.massflow



    # def calcRequired_massflow(self):
    #     req_massflow = 0
    #     for sf in self.subfields:
    #         req_massflow += sf.calc_required_massflow()
    #     self.req_massflow = req_massflow

    def set_operation_mode(self, mode = None):

        if mode is not None:
            self.operation_mode = mode
        else:
            if self.storage_available == True:
                pass
            else:
                if self.tout > self.tin:
                   self.operation_mode = "solarfield_heating"
                else:
                    self.operation_mode = "solarfield_not_heating"



    def print(self):

        for sf in self.subfields:
            for l in sf.loops:
                for s in l.scas:
                    for h in s.hces:
                        print("subfield: ", sf.name,
                              "Lazo: ",l.loop_order,
                              "SCA: ", s.sca_order,
                              "HCE: ", h.hce_order,
                              "tin", "=", h.tin,
                              "tout", "=", h.tout)



class Simulation(object):
    '''
    Definimos la clase simulacion para representar las diferentes
    pruebas que lancemos, variando el archivo TMY, la configuracion del
    site, la planta, el modo de operacion o el modelo empleado.
    '''

    def __init__(self, settings):
        self.ID =  settings['ID']
        self.simulation = settings['simulation']
        self.benchmark = settings['benchmark']
        self.datatype = settings['datatype']
        self.fastmode = settings['fastmode']
        self.tracking = True
        self.solarfield = None
        self.powersystem = None
        self.hotfluid = None
        self.coldfluid = None
        self.site = None
        self.model = None
        self.datasource = None
        self.powercycle = None


    def check_min_massflow(self):

        dri = self.solarfield.hce_settings['dri']
        t = self.solarfield.tin
        p = self.solarfield.pin
        re = self.solarfield.hce_settings['min_reynolds']
        loop_min_massflow = self.hotfluid.get_massflow_from_Reynolds(
                dri, t, p , re)
        solarfield_min_massflow = self.solarfield.total_loops * loop_min_massflow
        # self.solarfield.rated_massflow = self.powersystem.get_rated_massflow(
        #         self.powersystem.ratedpower, self.solarfield, self.hotfluid)
#        fluid_speed = 4 * loop_min_massflow / ( np.pi * hce_settings['dri']**2 *
#                self.hotfluid.get_density(self.solarfield.rated_tin,
#                                     self.solarfield.rated_pressure))
        if  solarfield_min_massflow > self.solarfield.rated_massflow:
            print("Too low massflow", solarfield_min_massflow ,">",
                  self.solarfield.rated_massflow)


    def runSimulation(self):

        for row in self.datasource.dataframe.iterrows():

            solarpos = pvlib.solarposition.get_solarposition(
                    row[0],
                    self.site.latitude,
                    self.site.longitude,
                    self.site.altitude,
                    pressure=row[1]['Pressure'],
                    temperature=row[1]['DryBulb'])

            if solarpos['zenith'][0] < 90:
                print("Generation mode ON. Solar Zenith = ", solarpos['zenith'][0])
                self.tracking = True
            else:
                print("Generation mode OFF. Solar Zenith = ", solarpos['zenith'][0])
                self.tracking = False

            if self.simulation:
                print("Running simulation for", self.site.name)
                self.simulate_solarfield(solarpos, row)

            if self.benchmark and self.datatype == 'field data':
                print("Running benchmark for", self.site.name)
                self.benchmarksolarfield(solarpos, row)

            self.plantperformance(row)
            self.gatherdata(row, solarpos)

        self.showresults()
        self.save_results()


    def simulate_solarfield(self, solarpos, row):

        flag_0 = datetime.now()
        self.solarfield.load_rated(self.hotfluid)
        self.solarfield.initialize_rated()
        self.check_min_massflow()
        self.prototype_loop.initialize()


        # Force recirculation massflow at night
        if solarpos['zenith'][0] > 90:
            self.prototype_loop.massflow =  (self.solarfield.recirculation_massflow /
                          self.solarfield.total_loops)

            self.prototype_loop.calc_loop_pr_for_massflow(
                row, solarpos, self.hotfluid, self.model)

        else:
            self.prototype_loop.calc_loop_pr_for_tout(
                row, solarpos, self.hotfluid, self.model)


        self.prototype_loop.tout = self.prototype_loop.scas[-1].hces[-1].tout
        self.prototype_loop.pout = self.prototype_loop.scas[-1].hces[-1].pout
        self.prototype_loop.pr_req_massflow = self.prototype_loop.get_loop_avg_pr()

        if  self.fastmode:
            flag_1 = datetime.now()
            delta_1 = flag_1 - flag_0

            for s in self.solarfield.subfields:
                for l in s.loops:
                    l.load_from_prototype_loop(self.prototype_loop)

                s.calc_massflow()
                s.calc_req_massflow()
                s.apply_temp_limitation(self.hotfluid)
                s.calc_pout()
                s.calc_tout(self.hotfluid)

            flag_2 = datetime.now()
            delta_2 = flag_2 - flag_1

        else:
            for s in self.solarfield.subfields:
                for l in s.loops:
                    l.initialize()


                    if solarpos['zenith'][0] > 90:
                        l.massflow =  (self.solarfield.recirculation_massflow /
                                      self.solarfield.total_loops)

                        l.calc_loop_pr_for_massflow(
                            row, solarpos, self.hotfluid, self.model)
                        l.pr_req_massflow = l.get_loop_avg_pr()

                    else:

                        # Start with the previous loop massflow, for a better convergence
                        if l.loop_order > 0:
                            l.massflow = l.subfield.loops[l.loop_order-1].massflow
                        else:
                            l.massflow = self.prototype_loop.massflow

                        l.calc_loop_pr_for_tout(
                            row, solarpos, self.hotfluid, self.model)
                        l.pr_req_massflow = l.get_loop_avg_pr()

                s.calc_massflow()
                s.calc_req_massflow()
                s.apply_temp_limitation(self.hotfluid)
                s.calc_tout(self.hotfluid)
                s.calc_pout()

            flag_3 = datetime.now()
            delta_3 = flag_3 - flag_0

        self.solarfield.calc_massflow()
        self.solarfield.calc_req_massflow()
        self.solarfield.calc_tout(self.hotfluid)
        self.solarfield.calc_pout()


    def benchmarksolarfield(self, solarpos, row):

        flag_0 = datetime.now()

        self.solarfield.load_actual(self.hotfluid, row)
        self.solarfield.initialize_actual()
        #self.prototype_loop.initialize()

        # self.solarfield.massflow = (self.prototype_loop.massflow *
        #                             self.solarfield.total_loops)

        flag_1 = datetime.now()
        delta_1 = flag_1 - flag_0

        if self.fastmode:

            # self.prototype_loop.calc_loop_pr(
            #     row, solarpos, self.hotfluid, self.model)


            for s in self.solarfield.subfields:

                # sf.initialize_actual(row)
                # sf.tin = sf.act_tin

                self.prototype_loop.act_massflow = s.act_massflow / len(s.loops)
                self.prototype_loop.massflow = self.prototype_loop.act_massflow
                self.prototype_loop.act_tin = s.act_tin
                self.prototype_loop.tin = s.act_tin
                self.prototype_loop.act_tout = s.act_tout
                self.prototype_loop.act_pin = s.act_pin
                self.prototype_loop.pin = s.act_pin
                self.prototype_loop.act_pout = s.act_pout

                self.prototype_loop.calc_loop_pr_for_massflow(
                    row, solarpos, self.hotfluid, self.model)

                self.prototype_loop.pr_act_massflow = self.prototype_loop.get_loop_avg_pr()

                for l in s.loops:
                    l.load_from_prototype_loop(self.prototype_loop)

                s.calc_massflow()
                s.calc_req_massflow()
                s.apply_temp_limitation(self.hotfluid)
                s.calc_pout()
                s.calc_tout(self.hotfluid)

            flag_2 = datetime.now()
            delta_2 = flag_2 - flag_1

        else:
            for s in self.solarfield.subfields:
                # sf.initialize_actual(row)
                # sf.tin = sf.act_tin
                for l in s.loops:
                    l.initialize_actual()
                    l.tin = l.act_tin
                    l.massflow = l.act_massflow
                    l.calc_loop_pr_for_massflow(
                        row, solarpos, self.hotfluid, self.model)
                    l.pr_act_massflow = l.get_loop_avg_pr()

                s.calc_massflow()
                s.calc_req_massflow()
                s.apply_temp_limitation(self.hotfluid)
                s.calc_pout()
                s.calc_tout(self.hotfluid)

            flag_3 = datetime.now()
            delta_3 = flag_3 - flag_0

        self.solarfield.calc_massflow()
        self.solarfield.calc_tout(self.hotfluid)
        self.solarfield.calc_pout()


    def plantperformance(self, row):

        self.solarfield.set_operation_mode()

        if self.solarfield.operation_mode == "subfield_not_heating":
            massflow_to_HE = 0
        else:
            massflow_to_HE = self.solarfield.massflow

        self.heatexchanger.set_hotfluid_in(massflow_to_HE,
                                         self.solarfield.tout,
                                         self.solarfield.pout)

        self.heatexchanger.set_coldfluid_in(
                self.powercycle.massflow, self.powercycle.tout, self.powercycle.pin) #PENDIENTE

        self.heatexchanger.set_thermalpowertransfered()
        self.powercycle.calc_pr_NCA(self.powercycle.tout) #Provisonal temperatura agua coondensador
        self.generator.calc_pr()
        print(row[0].strftime('%y/%m/%d %H:%M'), 'DNI: ', round(row[1]['DNI']),
              "MF:", round(self.solarfield.massflow),
              "Tin:", round(self.solarfield.tin),
              "Tout:", round(self.solarfield.tout))


    def gatherdata(self, row, solarpos):

        #TO-DO: CALCULOS PARA AGREGAR AL DATAFRAME

        if self.datatype == 'field data':
            mf_df = 0.0
            H_df = 0.0

            for col in self.datasource.dataframe.columns:
                if '.mf' in col:
                    mf_df  += self.datasource.dataframe.at[row[0], col]

                if '.tout' in col:
                    col_pout = col.replace('.tout', '.pout')
                    col_mf = col.replace('.tout', '.mf')
                    H_df +=  (self.hotfluid.get_deltaH(
                            self.datasource.dataframe.at[row[0], col],
                            self.datasource.dataframe.at[row[0], col_pout]) *
                            self.datasource.dataframe.at[row[0], col_mf])

            solarfield_df_tout = self.hotfluid.get_T(H_df / mf_df, self.solarfield.pout)
            self.datasource.dataframe.at[row[0], 'Massflow_df'] = mf_df
            self.datasource.dataframe.at[row[0], 'Tout_df'] = solarfield_df_tout

        for sf in self.solarfield.subfields:
            for l in sf.loops:

                self.datasource.dataframe.at[row[0], sf.name + str(l.loop_order).format("000") + '_tin'] = l.tin
                self.datasource.dataframe.at[row[0], sf.name + str(l.loop_order).format("000") + '_tout'] = l.tout
                self.datasource.dataframe.at[row[0], sf.name + str(l.loop_order).format("000") + '_mf'] = l.massflow
                self.datasource.dataframe.at[row[0], sf.name + str(l.loop_order).format("000") + '_pr_req_massflow'] = l.pr_req_massflow
                self.datasource.dataframe.at[row[0], sf.name + str(l.loop_order).format("000") + '_pr_act_massflow'] = l.pr_act_massflow
                self.datasource.dataframe.at[row[0], sf.name + str(l.loop_order).format("000") + '_act_tin'] = l.act_tin
                self.datasource.dataframe.at[row[0], sf.name + str(l.loop_order).format("000") + '_act_tout'] = l.act_tout
                self.datasource.dataframe.at[row[0], sf.name + str(l.loop_order).format("000") + '_act_mf'] = l.act_massflow
                self.datasource.dataframe.at[row[0], sf.name + str(l.loop_order).format("000") + '_req_mf'] = l.req_massflow

        self.datasource.dataframe.at[row[0], 'PTLoop_mf'] = self.prototype_loop.massflow
        self.datasource.dataframe.at[row[0], 'PTLoop_tin'] = self.prototype_loop.tin
        self.datasource.dataframe.at[row[0], 'PTLoop_tout'] = self.prototype_loop.tout
        self.datasource.dataframe.at[row[0], 'Pthermal'] = round(
                self.solarfield.get_thermalpoweroutput(self.hotfluid)/1e6,1)
        self.datasource.dataframe.at[row[0], 'Tin'] = round(self.solarfield.tin, 0)
        self.datasource.dataframe.at[row[0], 'Tout'] = round(self.solarfield.tout, 0)
        self.datasource.dataframe.at[row[0], 'MF'] = round(self.solarfield.massflow,0)
        self.heatexchanger.thermalpowertransfered
        self.datasource.dataframe.at[row[0], 'Pmec']  = round(
                self.heatexchanger.thermalpowertransfered *
                self.powercycle.pr / 1e6, 2)
        self.datasource.dataframe.at[row[0], 'Cycle.pr']  = round(self.powercycle.pr,2)
        self.datasource.dataframe.at[row[0], 'HE.pr']  = round(self.heatexchanger.pr,2)
        self.datasource.dataframe.at[row[0], 'Pelec']  = round(
                self.heatexchanger.thermalpowertransfered *
                self.powercycle.pr *
                self.generator.pr / 1e6, 2)
#        self.datasource.dataframe.at[row[0], 'iam'] = self.solarfield.subfields[3].loops[29].scas[0].get_IAM(aoi)
#        self.datasource.dataframe.at[row[0], 'azimuth'] = solarpos['azimuth'][0]
#        self.datasource.dataframe.at[row[0], 'zenith'] = solarpos['zenith'][0]


        #TO-DO
        estimated_pr = 0.0
        act_pr = 0.0
        estimated_tout = 0.0
        act_tout = 0.0
        rejected_solar_energy = 0.0

    def showresults(self):

        self.datasource.dataframe[[
                'DNI','MF', 'Pthermal', 'Pmec', 'Pelec', 'Tin', 'Tout']].plot(figsize=(20,10), linewidth=5, fontsize=20)
        plt.xlabel('Date', fontsize=20)

        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        # pd.set_option('display.max_colwidth', -1)

        if self.datatype == 'weather':
            print(self.datasource.dataframe[
                ['DNI', 'Tin', 'Tout','MF','Pthermal', 'Pmec',
                 'Pelec', 'MF']])

        elif self.datatype == 'field data':

            print(self.datasource.dataframe[
                ['DNI', 'Tin', 'Tout','MF','Pthermal', 'Pmec',
                 'Pelec', 'Massflow_df', 'Tout_df']])


        #print(self.datasource.dataframe)

    def save_results(self):


        try:
            initialdir = "./simulations_outputs/"
            prefix = datetime.today().strftime("%Y%m%d %H%M%S")
            filename = str(self.ID) + "_" + str(self.datatype)
            sufix = ".csv"

            path = initialdir + prefix + filename + sufix

            self.datasource.dataframe.to_csv(path, sep=';', decimal = ',')

        except Exception:
            raise
            print('Error saving results, unable to save file: %r', path)


    def get_solarposition(self, row):

        solarpos = pvlib.solarposition.get_solarposition(
                row[0],
                self.site.latitude,
                self.site.longitude,
                self.site.altitude,
                pressure=row[1]['Pressure'],
                temperature=row[1]['DryBulb'])

        return solarpos

    def testgeo(self):


        for row in self.datasource.dataframe.iterrows():

            solarpos = pvlib.solarposition.get_solarposition(
                    row[0],
                    self.site.latitude,
                    self.site.longitude,
                    self.site.altitude,
                    pressure = row[1]['Pressure'],
                    temperature=row[1]['DryBulb'])

            for s in self.prototype_loop.scas:
                aoi = s.get_aoi(solarpos)
                self.datasource.dataframe.at[row[0], 'sol.ze'] = round(solarpos['zenith'][0],2)
                self.datasource.dataframe.at[row[0], 'sol.az'] = round(solarpos['azimuth'][0],2)
                self.datasource.dataframe.at[row[0], 'aoi'] = round(aoi,2)

        self.datasource.dataframe[['sol.ze', 'sol.az', 'aoi']].plot(figsize=(20,10), linewidth=5, fontsize=20)
        plt.xlabel('Date', fontsize=20)

        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        pd.set_option('display.max_colwidth', -1)

        print(self.datasource.dataframe)


class Air(object):

    @classmethod
    def get_cinematic_viscosity(cls, t):

        # t: air temperature [K]
        return 8.678862e-11 * t**2 + 4.069284e-08 * t + -4.288741e-06

class Fluid(object):

    def test_fluid(self, tmax, tmin, p):

        data = []

        for tt in range(int(round(tmax)), int(round(tmin)), -5):
            data.append({'T': tt,
                         'P': p,
                         'cp': self.get_cp(tt, p),
                         'rho': self.get_density(tt, p),
                         'mu': self.get_dynamic_viscosity(tt, p),
                         'kt': self.get_thermal_conductivity(tt, p),
                         'H': self.get_deltaH(tt, p),
                         'T-H': self.get_T(self.get_deltaH(tt, p), p)})
        df = pd.DataFrame(data)
        print(round(df, 6))

class Fluid_CoolProp(Fluid):

    def __init__(self, settings = None):

        self.name = settings['name']
        self.tmax = settings['tmax']
        self.tmin = settings['tmin']
        self.coolpropID = settings['CoolPropID']

    def get_density(self, t, p):

        if t > self.tmax:
            t = self.tmax

        return PropsSI('D','T',t,'P', p, self.coolpropID)

    def get_dynamic_viscosity(self, t, p):

        if t > self.tmax:
            t = self.tmax

        return  PropsSI('V','T',t,'P', p, self.coolpropID)

    def get_cp(self, t, p):

        if t > self.tmax:
            t = self.tmax
        return PropsSI('C','T',t,'P', p, self.coolpropID)

    def get_thermal_conductivity(self, t, p):
        ''' Saturated Fluid conductivity at temperature t '''

        if t > self.tmax:
            t = self.tmax

        return PropsSI('L','T',t,'P', p, self.coolpropID)

    def get_deltaH(self, t, p):

        if t > self.tmax:
            t = self.tmax

        CP.set_reference_state(self.coolpropID,'ASHRAE')

        deltaH = PropsSI('H','T',t ,'P', p, self.coolpropID)

        CP.set_reference_state(self.coolpropID, 'DEF')

        return deltaH

    def get_T(self, h, p):

        # if t > self.tmax:
        #     t = self.tmax

        CP.set_reference_state(self.coolpropID,'ASHRAE')
        temperature = PropsSI('T', 'H', h, 'P', p, self.coolpropID)
        CP.set_reference_state(self.coolpropID, 'DEF')

        return temperature


    def get_Reynolds(self, dri, t, p, massflow):

        if t > self.tmax:
            t = self.tmax

        return massflow * np.pi * (dri**3) /( 4 * self.get_density(t, p))

    def get_massflow_from_Reynolds(self, dri, t, p, re):

        if t > self.tmax:
            t = self.tmax

        return re * np.pi * dri * self.get_dynamic_viscosity(t, p) / 4


class Fluid_Tabular(Fluid):

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

        # Dowtherm A.pdf, 2.2 Single Phase Liquid Properties. pg. 8.

        rho0, rho1, rho2, rho3, rho4, rho5 = tuple(self.rho)

        return (rho0 + rho1 * t + rho2 * t**2 + rho3 * t**3 +
                rho4 * t**4 + rho5 * t**5) * (p* 1.0e-4)**1.0e-3

    def get_dynamic_viscosity(self, t, p):

        mu0, mu1, mu2, mu3, mu4, mu5 = tuple(self.mu)

        if t > self.tmax:
            t= self.tmax

        return (mu0 + mu1 * t + mu2 * t**2 + mu3 * t**3 +
                mu4 * t**4 + mu5 * t**5)

    def get_cp(self, t, p):

        cp0, cp1, cp2, cp3, cp4, cp5 = tuple(self.cp)

        return (cp0 + cp1 * t + cp2 * t**2 + cp3 * t**3 +
                cp4 * t**4 + cp5 * t**5)

    def get_thermal_conductivity(self, t, p):
        ''' Saturated Fluid conductivity at temperature t '''

        kt0, kt1, kt2, kt3, kt4, kt5 = tuple(self.kt)

        kt =(kt0 + kt1 * t + kt2 * t**2 + kt3 * t**3 +
                kt4 * t**4 + kt5 * t**5)

        if kt < 0:
            kt = 0

        return kt


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

        if t > self.tmax:
            t= self.tmax

        return (4 * massflow /
                (np.pi * dri * self.get_dynamic_viscosity(t,p)))

    def get_massflow_from_Reynolds(self, dri, t, p, re):

        if t > self.tmax:
            t= self.tmax

        return re * np.pi * dri * self.get_dynamic_viscosity(t,p) / 4


class Weather(object):

    def __init__(self, settings):
        self.filename = settings['filename']
        self.filepath = settings['filepath']
        self.file = self.filepath + self.filename
        self.dataframe = None
        self.site = None
        self.weatherdata = None

        self.openWeatherDataFile(self.file)

        self.dataframe = self.weatherdata[0]
        self.site = self.weatherdata[1]

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

        for c in self.dataframe.columns:
            if (c == 'DryBulb') or (c == 'DewPoint'): # From Celsius Degrees to K
                self.dataframe[c] *= 0.1
                self.dataframe[c] += 273.15
            if c=='Pressure': # from mbar to Pa
                self.dataframe[c] *= 1e2

    def filter_columns(self):

        needed_columns = ['DNI', 'DryBulb', 'DewPoint', 'Wspd', 'Wdir', 'Pressure']
        columns_to_drop = []
        for c in  self.dataframe.columns:
            if c not in needed_columns:
                columns_to_drop.append(c)
        self.dataframe.drop(columns = columns_to_drop, inplace = True)

    def site_to_dict(self):
        '''
        pvlib.iotools realiza modificaciones en los nombres de las columnas.

        Source,Location ID,City,State,Country,Latitude,Longitude,Time Zone,Elevation,Local Time Zone,Dew Point Units,DHI Units,DNI Units,GHI Units,Temperature Units,Pressure Units,Wind Direction Units,Wind Speed,Surface Albedo Units,Version'''

        return {"name": self.site['City'],
                "latitude": self.site['latitude'],
                "longitude": self.site['longitude'],
                "altitude": self.site['altitude']}

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
        #dateparse = lambda x: pd.datetime.strptime(x, '%YYYY/%m/%d %H:%M')

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

        self.dataframe.index = pd.to_datetime(self.dataframe.index, infer_datetime_format=True)

    def change_units(self):

        for c in self.dataframe.columns:
            if ('.t' in c) or ('DryBulb' in c) or ('Dew' in c):
                self.dataframe[c] += 273.15 # From Celsius Degrees to K
            if '.p' in c:
                self.dataframe[c] *= 1e5 # From Bar to Pa
            if 'Pressure' in c:
                self.dataframe[c] *= 1e2 # From mBar to Pa

    def rename_columns(self):

        # Replace tags  with names as indicated in configuration file
        # (field_data_file: tags)

        rename_dict = dict(zip(self.tags.values(), self.tags.keys()))
        self.dataframe.rename(columns = rename_dict, inplace = True)


        # Remove unnecessary columns

        columns_to_drop = []
        for c in  self.dataframe.columns:
            if c not in self.tags.keys():
                columns_to_drop.append(c)
        self.dataframe.drop(columns = columns_to_drop, inplace = True)


class PowerSystem(object):
    '''
    Power Plant as a set composed by a solarfield, a HeatExchanger, a PowerCycle and a BOPSystem

    '''

    # Calculo de potencia de salida (positiva o negativa)
    #y potencia derivada a BOB en base al estado de la planta

    def __init__(self, settings):

        self.ratedpower = settings['ratedpower']
        self.subfield_to_exchanger_pr = settings['solarfield_to_exchanger_pr']
        self.exchanger_pr = settings['exchanger_pr']
        self.exchanger_to_turbogroup_pr = settings['exchanger_to_turbogroup_pr']
        self.steam_cycle_pr = settings['steam_cycle_pr']
        self.turbogenerator_pr = settings['turbogenerator_pr']


    def get_poweroutput(self):
        pass

    def get_powerinput(self):
        pass

    def get_thermalpower(self, ratedpower):

        return ratedpower / (self.subfield_to_exchanger_pr *
                                  self.exchanger_pr *
                                  self.exchanger_to_turbogroup_pr *
                                  self.steam_cycle_pr *
                                  self.turbogenerator_pr)

    def get_rated_massflow(self, ratedpower, solarfield, hotfluid, coldfluid = None):

        HH = hotfluid.get_deltaH(solarfield.rated_tout, solarfield.rated_pout)
        HL = hotfluid.get_deltaH(solarfield.rated_tin, solarfield.rated_pin)
        rated_massflow = self.get_thermalpower(self.ratedpower) / (HH - HL)

        return rated_massflow

class HeatExchanger(object):

    def __init__(self, settings, hotfluid, coldfluid):

        self.pr = settings['pr']
        self.thermalpowertransfered = 0.0

        self.hotfluidmassflow = 0.0
        self.hotfluidtin = 0.0
        self.hotfluidtout = 393 + 273.15
        self.hotfluidpin = 1500000
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

    def set_thermalpowertransfered(self):

        #provisional
        pinch = 10.0

        HH = self.hotfluid.get_deltaH(self.hotfluidtin, self.hotfluidpin)
        HL = self.hotfluid.get_deltaH(self.coldfluidtin + pinch, self.hotfluidpout)

        self.thermalpowertransfered = self.pr * ((HH-HL)*self.hotfluidmassflow)
#        self.thermalpowertransfered = (
#                self.pr *
#                self.hotfluidmassflow *
#                self.hotfluid.get_cp(self.hotfluidtin, self.hotfluidpin) *
#                self.hotfluidtin - hftout)


    def hot_fluid_tout():
        return


class HeatStorage(object):

    def __init__(self, settings):

        self.name = settings['heatstorage']

    def set_fluid_tout(self):
        pass

    def hot_fluid_tout():
        return

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
        self.pin = settings['rated_pin']
        self.tin = settings['rated_tin']
        self.pout = settings['rated_pout']
        self.tout = settings['rated_tout']
        self.massflow = settings['rated_massflow']
        self.pr = 0.0

    def calc_pr_NCA(self, tout): #Novikov-Curzon-Ahlbor

        self.pr = 1 - np.sqrt((tout)/(self.tin))

class Generator(object):

    def __init__(self, settings):

        self.name = settings['name']
        self.pr = settings['pr']

    def calc_pr(self):

        #TO-DO: Desasrrollar curvas pr-carga-temp.ext por ejemplo.
        pass



class Site(object):
    def __init__(self, settings):

        self.name = settings['name']
        self.latitude = settings['latitude']
        self.longitude = settings['longitude']
        self.altitude = settings['altitude']


class HCEScatterMask(object):

# %matplotlib inline

# import matplotlib.pyplot as plt
# import numpy as np
# from scipy import stats
# import seaborn as sns

# np.random.seed(2016) # replicar random

# # parametros esteticos de seaborn
# sns.set_palette("deep", desat=.6)
# sns.set_context(rc={"figure.figsize": (8, 4)})

# mu, sigma = 0, 0.2 # media y desvio estandar
# datos = np.random.normal(mu, sigma, 1000) #creando muestra de datos

# # histograma de distribución normal.
# cuenta, cajas, ignorar = plt.hist(datos, 20)
# plt.ylabel('frequencia')
# plt.xlabel('valores')
# plt.title('Histograma')
# plt.show()


    def __init__(self, solarfield_settings, hce_mask_settings):

        self.matrix = dict()

        for sf in solarfield_settings['subfields']:
            self.matrix[sf["name"]]=[]
            for l in range(sf.get('loops')):
                self.matrix[sf["name"]].append([])
                for s in range(solarfield_settings['loop']['scas']):
                    self.matrix[sf["name"]][-1].append([])
                    for h in range(solarfield_settings['loop']['hces']):
                        self.matrix[sf["name"]][-1][-1].append(hce_mask_settings)

    def applyMask(self, solarfield):

        for sf in solarfield.subfields:
            for l in sf.loops:
                for s in l.scas:
                    for h in s.hces:
                        for k in self.matrix[sf.name][l.loop_order][s.sca_order][h.hce_order].keys():
                            h.parameters[k] *= float(self.matrix[sf.name][l.loop_order][s.sca_order][h.hce_order][k])


class SCAScatterMask(object):


    def __init__(self, solarfield_settings, sca_mask_settings):

        self.matrix = dict()

        for sf in solarfield_settings['subfields']:
            self.matrix[sf["name"]]=[]
            for l in range(sf.get('loops')):
                self.matrix[sf["name"]].append([])
                for s in range(solarfield_settings['loop']['scas']):
                    self.matrix[sf["name"]][-1].append(sca_mask_settings)

    def applyMask(self, solarfield):

        for sf in solarfield.subfields:
            for l in sf.loops:
                for s in l.scas:
                    for k in self.matrix[sf.name][l.loop_order][s.sca_order].keys():
                        s.parameters[k] *= float(self.matrix[sf.name][l.loop_order][s.sca_order][k])








