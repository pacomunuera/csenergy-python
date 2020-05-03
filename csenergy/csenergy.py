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
import time
import os.path
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
#  import pint


class Model:

    def calc_pr(self):
        pass

    def get_hext_eext(self, hce, reext, tro, wind):

        eext = 0.
        hext = 0.

        if hce.parameters['Name'] == 'Solel UVAC 2/2008':
            pass

        elif hce.parameters['Name'] == 'Solel UVAC 3/2010':
            pass

        elif hce.parameters['Name'] == 'Schott PTR70':
            pass

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


class ModelBarbero4thOrder(Model):

    def __init__(self, settings):

        self.parameters = settings
        self.max_err_t = self.parameters['max_err_t']
        self.max_err_tro = self.parameters['max_err_tro']
        self.max_err_pr = self.parameters['max_err_pr']

    def calc_pr(self, hce, htf, row, qabs=None):

        flag_0 = datetime.now()

        if qabs is not None:
            qabs = qabs
        else:
            qabs = hce.qabs

        tin = hce.tin
        tf = hce.tin  # HTF bulk temperature
        tri = hce.tin  #  Absorber tube inner surface temperature
        massflow = hce.sca.loop.massflow
        wspd = row[1]['Wspd']  #  Wind speed
        text = row[1]['DryBulb']  #  Dry bulb ambient temperature
        sigma = sc.constants.sigma  #  Stefan-Bolztmann constant
        dro = hce.parameters['Absorber tube outer diameter']
        dri = hce.parameters['Absorber tube inner diameter']
        dgo = hce.parameters['Glass envelope outer diameter']
        dgi = hce.parameters['Glass envelope inner diameter']
        L = hce.parameters['Length']
        A = hce.sca.parameters['Aperture']
        x = 1 #  Calculation grid fits hce longitude

        #  HCE wall thermal conductivity
        krec = hce.get_krec(tf)

        #  Specific Capacity
        cp = htf.get_cp(tf, hce.pin)

        #  Internal transmission coefficient.
        hint = hce.get_hint(tf, hce.pin, htf)

        #  Internal heat transfer coefficient. Eq. 3.22 Barbero2016
        urec = 1 / ((1 / hint) + (dro * np.log(dro / dri)) / (2 * krec))

        #  We suppose thermal performance, pr = 1, at first
        pr = 1.0
        tro = tf + qabs * pr / urec

        #  HCE emittance
        eext = hce.get_emittance(tro, wspd)
        #  External Convective Heat Transfer equivalent coefficient
        hext = hce.get_hext(wspd)

        #  Thermal power lost. Eq. 3.23 Barbero2016
        qperd = sigma * eext * (tro**4 - text**4) + hext * (tro - text)

        #  Thermal power lost througth  bracktets
        qlost_brackets = hce.get_qlost_brackets(tf, text)

        #  Critical Thermal power loss. Eq. 3.50 Barbero2016
        qcrit = sigma * eext * (tf**4 - text**4) + hext * (tf - text)

        #  Critical Internal heat transfer coefficient, Eq. 3.51 Barbero2016
        ucrit = 4 * sigma * eext * tf**3 + hext

        #  Transmission Units Number, Ec. 3.30 Barbero2016
        NTU = urec * x * L * np.pi * dro / (massflow * cp)

        if qabs > qcrit:

            #  We use Barbero2016's simplified model aproximation
            #  Eq. 3.63 Barbero2016
            fcrit = 1 / (1 + (ucrit / urec))

            #  Eq. 3.71 Barbero2016
            pr = fcrit * (1 - qcrit / qabs)

            errtro = 1.
            errpr = 1.
            step = 0

            while (errtro > self.max_err_tro or errpr > self.max_err_pr):

                step += 1
                flag_1 = datetime.now()

                #  Eq. 3.32 Barbero2016
                f0 = qabs / (urec * (tf - text))

                #  Eq. 3.34 Barbero2016
                f1 = ((4 * sigma * eext * text**3) + hext) / urec
                f2 = 6 * (text**2) * (sigma * eext / urec) * (qabs / urec)
                f3 = 4 * text * (sigma * eext / urec) * ((qabs / urec)**2)
                f4 = (sigma * eext / urec) * ((qabs / urec)**3)


                #  We solve Eq. 3.36 Barbero2016 in order to find the performance
                #  of the first section of the HCE (at entrance)

                pr0 = pr
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

                #  Eq. 3.37 Barbero2016
                z = pr0 + (1 / f0)

                #  Eq. 3.40, 3.41 & 3.42 Babero2016
                g1 = 1 + f1 + 2 * f2 * z + 3 * f3 * z**2 + 4 * f4 *z**3
                g2 = 2 * f2 + 6 * f3 * z + 12 * f4 *z**2
                g3 = 6 * f3 + 24 * f4 * z

                #  Eq. 3.39 Barbero2016
                pr2 = ((pr0 * g1 / (1 - g1)) * (1 / (NTU * x)) *
                       (np.exp((1 - g1) * NTU * x / g1) - 1) -
                       (g2 / (6 * g1)) * (pr0 * NTU * x)**2 -
                       (g3 / (24 * g1) * (pr0 * NTU * x)**3))

                errpr = abs(pr2-pr)
                pr = pr2
                hce.pr = pr
                hce.set_tout(htf)
                hce.set_pout(htf)
                tf = 0.5 * (hce.tin + hce.tout)

                #  HCE wall thermal conductivity
                krec = hce.get_krec(tf)

                #  Specific Capacity
                cp = htf.get_cp(tf, hce.pin)

                #  Internal transmission coefficient.
                hint = hce.get_hint(tf, hce.pin, htf)

                #  Internal heat transfer coefficient. Eq. 3.22 Barbero2016
                urec = 1 / ((1 / hint) + (dro * np.log(dro / dri)) / (2 * krec))

                #  HCE emittance
                eext = hce.get_emittance(tro, wspd)

                #  External Convective Heat Transfer equivalent coefficient
                hext = hce.get_hext(wspd)

                #  We calculate tro again.

                tro2 = tf + qabs * pr / urec
                errtro = abs(tro2-tro)
                tro = tro2

                #  Thermal power loss. Eq. 3.23 Barbero2016
                qperd = sigma * eext * (tro2**4 - text**4) + hext * (tro2 - text)

                #  Increase qperd with  the thermal power loss througth  bracktets
                qlost_brackets = hce.get_qlost_brackets(tf, text)
                qperd = qperd + qlost_brackets

                #  Critical Thermal power loss. Eq. 3.50 Barbero2016
                qcrit = sigma * eext * (tf**4 - text**4) + hext * (tf - text)

                #  Critical Internal heat transfer coefficient, Eq. 3.51 Barbero2016
                ucrit = 4 * sigma * eext * tf**3 + hext

                #  Transmission Units Number, Ec. 3.30 Barbero2016
                NTU = urec * x * L * np.pi * dro / (massflow * cp)

            hce.pr = hce.pr * (1 - qlost_brackets / qabs)
            hce.qperd = qperd  * hce.parameters['Length']
            hce.qlost_brackets = qlost_brackets * hce.parameters['Length']
            hce.qabs = qabs * hce.parameters['Length']

        else:
            hce.pr = 0.0
            hce.qperd = qperd * hce.parameters['Length']
            hce.qabs = qabs * hce.parameters['Length']
            hce.qlost_brackets =  qlost_brackets * hce.parameters['Length']
            hce.set_tout(htf)
            hce.set_pout(htf)
            flag_3 = datetime.now()


class ModelBarbero1stOrder(Model):

    def __init__(self, settings):

        self.parameters = settings
        self.max_err_t = self.parameters['max_err_t']

    def calc_pr(self, hce, htf, row, qabs = None):

        flag_0 = datetime.now()

        if qabs is not None:
            qabs = qabs
        else:
            qabs = hce.qabs

        tin = hce.tin
        tf = hce.tin  # HTF bulk temperature
        tri = hce.tin  #  Absorber tube inner surface temperature
        massflow = hce.sca.loop.massflow
        wspd = row[1]['Wspd']  #  Wind speed
        text = row[1]['DryBulb']  #  Dry bulb ambient temperature
        sigma = sc.constants.sigma  #  Stefan-Bolztmann constant
        dro = hce.parameters['Absorber tube outer diameter']
        dri = hce.parameters['Absorber tube inner diameter']
        dgo = hce.parameters['Glass envelope outer diameter']
        dgi = hce.parameters['Glass envelope inner diameter']
        L = hce.parameters['Length']
        A = hce.sca.parameters['Aperture']
        x = 1 #  Calculation grid fits hce longitude

        #  HCE wall thermal conductivity
        krec = hce.get_krec(tf)

        #  Specific Capacity
        cp = htf.get_cp(tf, hce.pin)

        #  Internal transmission coefficient.
        hint = hce.get_hint(tf, hce.pin, htf)

        #  Internal heat transfer coefficient. Eq. 3.22 Barbero2016
        urec = 1 / ((1 / hint) + (dro * np.log(dro / dri)) / (2 * krec))

        #  We suppose performance, pr = 1, at first
        pr = 1.0
        tro = tf + qabs * pr / urec

        #  HCE emittance
        eext = hce.get_emittance(tro, wspd)
        #  External Convective Heat Transfer equivalent coefficient
        hext = hce.get_hext(wspd)

        #  Thermal power lost. Eq. 3.23 Barbero2016
        qperd = sigma * eext * (tro**4 - text**4) + hext * (tro - text)

        #  Thermal power lost througth  bracktets
        qlost_brackets = hce.get_qlost_brackets(tf, text)

        #  Critical Thermal power loss. Eq. 3.50 Barbero2016
        qcrit = sigma * eext * (tf**4 - text**4) + hext * (tf - text)

        #  Critical Internal heat transfer coefficient, Eq. 3.51 Barbero2016
        ucrit = 4 * sigma * eext * tf**3 + hext

        #  Ec. 3.63
        ## fcrit = (1 / ((4 * eext * tfe**3 / urec) + (hext / urec) + 1))
        fcrit = 1 / (1 + (ucrit / urec))

        #  Ec. 3.64
        Aext = np.pi * dro * x / 2  # Pendiente de confirmar
        NTUperd = ucrit * Aext / (massflow * cp)

        if qabs > qcrit:

            hce.pr = ((1 - (qcrit / qabs)) *
                  (1 / (NTUperd * x)) *
                  (1 - np.exp(-NTUperd * fcrit * x)))

            hce.pr = hce.pr * (1 - qlost_brackets / qabs)

        else:
            hce.pr = 0

        hce.qabs = qabs * hce.parameters['Length']
        hce.qperd = qperd * hce.parameters['Length']
        hce.qlost_brackets = qlost_brackets * hce.parameters['Length']
        hce.set_tout(htf)
        hce.set_pout(htf)


class ModelBarberoSimplified(Model):

    def __init__(self, settings):

        self.parameters = settings
        self.max_err_t = self.parameters['max_err_t']

    def calc_pr(self, hce, htf, row, qabs=None):

        if qabs is not None:
            qabs = qabs
        else:
            qabs = hce.qabs

        tin = hce.tin
        tf = hce.tin  # HTF bulk temperature
        tri = hce.tin  #  Absorber tube inner surface temperature
        massflow = hce.sca.loop.massflow
        wspd = row[1]['Wspd']  #  Wind speed
        text = row[1]['DryBulb']  #  Dry bulb ambient temperature
        sigma = sc.constants.sigma  #  Stefan-Bolztmann constant
        dro = hce.parameters['Absorber tube outer diameter']
        dri = hce.parameters['Absorber tube inner diameter']
        dgo = hce.parameters['Glass envelope outer diameter']
        dgi = hce.parameters['Glass envelope inner diameter']
        L = hce.parameters['Length']
        A = hce.sca.parameters['Aperture']
        x = 1 #  Calculation grid fits hce longitude

        #  HCE wall thermal conductivity
        krec = hce.get_krec(tf)

        #  Specific Capacity
        cp = htf.get_cp(tf, hce.pin)

        #  Internal transmission coefficient.
        hint = hce.get_hint(tf, hce.pin, htf)

        #  Internal heat transfer coefficient. Eq. 3.22 Barbero2016
        urec = 1 / ((1 / hint) + (dro * np.log(dro / dri)) / (2 * krec))

        #  We suppose performance, pr = 1, at first
        pr = 1.0
        tro = tf + qabs * pr / urec

        #  HCE emittance
        eext = hce.get_emittance(tro, wspd)
        #  External Convective Heat Transfer equivalent coefficient
        hext = hce.get_hext(wspd)

        #  Thermal power loss. Eq. 3.23 Barbero2016
        qperd = sigma * eext * (tro**4 - text**4) + hext * (tro - text)

        #  Thermal power lost througth  bracktets
        qlost_brackets = hce.get_qlost_brackets(tf, text)

        #  Critical Thermal power loss. Eq. 3.50 Barbero2016
        qcrit = sigma * eext * (tf**4 - text**4) + hext * (tf - text)

        #  Critical Internal heat transfer coefficient, Eq. 3.51 Barbero2016
        ucrit = 4 * sigma * eext * tf**3 + hext

        #  Ec. 3.63
        ## fcrit = (1 / ((4 * eext * tfe**3 / urec) + (hext / urec) + 1))
        fcrit = 1 / (1 + (ucrit / urec))

        if qabs > qcrit:

            hce.pr = fcrit * (1 - (qcrit / qabs))

            hce.pr = hce.pr * (1 - qlost_brackets / qabs)

        else:
            hce.pr = 0


        hce.qabs = qabs * hce.parameters['Length']
        hce.qperd = qperd * hce.parameters['Length']
        hce.qlost_brackets = qlost_brackets * hce.parameters['Length']
        hce.set_tout(htf)
        hce.set_pout(htf)


# class ModelHottelWhilier(Model):

#         def __ini__(self, simulation):
#             super(Model, self).__init__(simulation)


# class ModelNaumFrainderaich(Model):

#         def __ini__(self, simulation):
#             super(Model, self).__init__(simulation)


# class ModelASHRAE(Model):

#         def __ini__(self, simulation):
#             super(Model, self).__init__(simulation)


class HCE(object):

    def __init__(self, sca, hce_order, settings):

        self.sca = sca  # SCA in which the HCE is located
        self.hce_order = hce_order  # Relative position of the HCE in the SCA
        self.parameters = dict(settings)  # Set of parameters of the HCE
        self.tin = 0.0  # Temperature of the HTF when enters in the HCE
        self.tout = 0.0  # Temperature of the HTF when goes out the HCE
        self.pr = 0.0  # Thermal performance of the HCE
        self.pr_opt = 0.0  # Optical performance of the HCE+SCA set
        self.qabs = 0.0  # Thermal which reach the aboserber tube
        self.qperd = 0.0  # Thermal power lost througth out the HCE
        self.qlost_brackets = 0.0  # Thermal power lost in brackets and arms
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

    def set_tout(self, htf):

        HL = htf.get_deltaH(self.tin, self.sca.loop.pin)

        if self.qabs > 0:

            h = (self.qabs * self.pr * self.pr_opt / self.sca.loop.massflow)

        else:
            h = -self.qperd / self.sca.loop.massflow

        self.tout = htf.get_T(HL + h, self.sca.loop.pin)


    def set_pout(self, htf):

        # TO-DO CÁLCULLO PERDIDA DE CARGA:
        # Ec. Colebrook-White para el cálculo del factor de fricción de Darcy

        re_turbulent = 4000

        k = self.parameters['Inner surface roughness']
        D = self.parameters['Absorber tube inner diameter']
        re = htf.get_Reynolds(D, self.tin, self.pin, self.sca.loop.massflow)


        if re < re_turbulent:
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

        rho = htf.get_density(self.tin, self.pin)
        v = 4 * self.sca.loop.massflow / (rho * np.pi * D**2)
        g = sc.constants.g

        # Ec. Darcy-Weisbach para el cálculo de la pérdida de carga
        deltap_mcl = darcy_factor * (self.parameters['Length'] / D) * \
            (v**2 / (2 * g))
        deltap = deltap_mcl * rho * g
        self.pout = self.pin - deltap

    def set_pr_opt(self, aoi):

        IAM = self.sca.get_IAM(aoi)
        pr_opt_peak = self.get_pr_opt_peak()
        self.pr_opt = pr_opt_peak * IAM

    def set_qabs(self, aoi, solarpos, row):
        ''' Total solar power that reach de absorber tube per longitude unit'''

        dni = row[1]['DNI']
        cg = (self.sca.parameters['Aperture'] /
              (np.pi*self.parameters['Absorber tube outer diameter']))

        pr_geo = self.get_pr_geo(aoi)
        pr_shadows = self.get_pr_shadows(solarpos)
        #  Ec. 3.20 Barbero
        self.qabs = (self.pr_opt * cg * dni * pr_geo * pr_shadows)

    def get_krec(self, t):

        # Ec. 4.22 Conductividad para el acero 321H, ver otras opciones.
        return  0.0153 * (t - 273.15) + 14.77

    def get_previous(self):

        return self.sca.hces[self.hce_order-1]

    def get_index(self):

        if hasattr(self.sca.loop, 'subfield'):
            index = [self.sca.loop.subfield.name,
                self.sca.loop.loop_order,
                self.sca.sca_order,
                self.hce_order]
        else:
            index = ['BL',
                self.sca.sca_order,
                self.hce_order]

        return index

    def get_pr_opt_peak(self):

        alpha = self.get_absorptance()
        tau = self.get_transmittance()
        rho = self.sca.parameters['Reflectance']
        gamma = self.sca.get_solar_fraction()

        pr_opt_peak = alpha * tau * rho * gamma

        if pr_opt_peak > 1 or pr_opt_peak < 0:
            print("ERROR", pr_opt_peak)

        return pr_opt_peak


    def get_pr_geo(self, aoi):

        if aoi > 90:
            pr_geo = 0.0

        else:
            # Llamado "bordes" en Tesis. Pérdidas de los HCE de cabecera según aoi
            sca_unused_length = (self.sca.parameters["Focal Len"] *
                                 np.tan(np.radians(aoi)))

            unused_hces = sca_unused_length // self.parameters["Length"]

            unused_part_of_hce = ((sca_unused_length %
                                   self.parameters["Length"]) /
                                  self.parameters["Length"])

            if self.hce_order < unused_hces:
                pr_geo = 0.0

            elif self.hce_order == unused_hces:
                pr_geo = 1 - ((sca_unused_length % self.parameters["Length"]) /
                                  self.parameters["Length"])
            else:
                pr_geo = 1.0

            if pr_geo > 1.0 or pr_geo < 0.0:
                print("ERROR pr_geo out of limits", pr_geo)

        return pr_geo


    def get_pr_shadows(self, solarpos):

        # Llamado "sombras" en Tesis. Pérdidas por sombras. ¿modelar sobre el SCA?
        # Sombras debidas a otros lazos

        shadowing = (1 -
                     np.sin(np.radians(abs(solarpos['elevation'][0]))) *
                     self.sca.loop.parameters['row_spacing'] /
                     self.sca.parameters['Aperture'])

        if shadowing < 0.0 or shadowing > 1.0:
            shadowing = 0.0

        if solarpos['elevation'][0] < 0:
            shadowing = 1.0

        pr_shadows = 1 - shadowing


        if  pr_shadows > 1 or  pr_shadows < 0:
            print("ERROR",  pr_shadows)

        return pr_shadows


    def get_hext(self, wspd):

        #  TO-DO:

        return 0.0

    def get_hint(self, t, p, fluid):


        #  Prandtl number
        prf = fluid.get_prandtl(t, p)

        kf = fluid.get_thermal_conductivity(t, p)

        mu = fluid.get_dynamic_viscosity(t, p)

        dri = self.parameters['Absorber tube inner diameter']

        #  Reynolds number for absorber tube inner diameter, dri
        redri = 4 * self.sca.loop.massflow / (mu * np.pi * dri)

        #  Nusselt num. Dittus-Boelter correlation. Eq. 4.14 Barbero2016
        nudb = 0.023 * redri**0.8 * prf** 0.4

        #  Internal transmission coefficient.
        hint = kf * nudb / dri

        return hint

    def get_emittance(self, tro, wspd):


        #  Eq. 5.2 Barbero
        eext = (self.parameters['Absorber emittance factor A0'] +
                self.parameters['Absorber emittance factor A1'] *
                (tro - 273.15))
        """
        Lineal Increase if wind speed lower than 4 m/s up to 1% at 4 m/s
        Lineal increase over 4 m/s up to 2% at 7 m/s
        """
        if wspd <4:
            eext = eext * (1 + 0.01 * wspd / 4)

        else:
            eext = eext * (1 + 0.01 * (0.3333 * wspd - 0.3333))

        return eext


    def get_absorptance(self):

        #
        alpha = self.parameters['Absorber absorptance']

        return alpha


    def get_transmittance(self):

        #dependiendo si hay vídrio y si hay vacío
        tau = self.parameters['Envelope transmittance']

        return tau


    def get_reflectance(self):

        return self.sca.parameters['Reflectance']


    def get_qlost_brackets(self, tf, text):

        #  Ec. 4.12

        n = self.parameters['Length'] / self.parameters['Brackets'] + \
            + (self.hce_order == 0)
        pb = 0.2032
        acsb = 1.613e-4
        kb = 48
        hb = 20
        tbase = tf - 10

        L = self.parameters['Length']

        return n * (np.sqrt(pb * kb * acsb * hb) * (tbase - text)) / L


class SCA(object):


    def __init__(self, loop, sca_order, settings):

        self.loop = loop
        self.sca_order = sca_order
        self.hces = []
        self.status = 'focused'
        self.tracking_angle = 0.0
        self.parameters = dict(settings)
        self.surface_tilt = 0.0
        self.surface_azimuth = 180.0


    def get_solar_fraction(self):

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
                             self.parameters['Availability'])
        else:
            solarfraction = 1.0

        if  solarfraction > 1 or  solarfraction < 0:
            print("ERROR",  solarfraction)

        return solarfraction


    def get_IAM(self, theta):

        F0 = self.parameters['IAM Coefficient F0']
        F1 = self.parameters['IAM Coefficient F1']
        F2 = self.parameters['IAM Coefficient F2']


        if (theta > 0 and theta < 80):
            kiam = (F0 + (F1 * np.radians(theta) + F2 * np.radians(theta)**2) /
                    np.cos(np.radians(theta)))

            if kiam > 1:
                kiam = 2 - kiam

            if kiam < 0:
                kiam = 0.0

        else:
            kiam = 0.0

        if  kiam > 1 or  kiam < 0:
            print("ERROR",  kiam, theta)

        return kiam


    def get_aoi(self, solarpos):

        sigmabeta = 0.0
        beta0 = 0.0

        if self.loop.parameters['Tracking Type'] == 1: # N-S single axis tracker
            if solarpos['azimuth'][0] > 0 and solarpos['azimuth'][0] <= 180:
                surface_azimuth = 90 # Surface facing east
            else:
                surface_azimuth = 270 # Surface facing west
        elif self.loop.parameters['Tracking Type'] == 2:  # E-W single axis tracker
            surface_azimuth = 180  # Surface facing the equator

        #  En esta fórmula asumo que el seguimiento del SCA es perfecto
        #  pero hay que ver la posibilidad de modelar cierto error o desfase
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


class __Loop__(object):


    def __init__(self, settings):

        self.scas = []
        self.parameters = settings

        self.tin = 0.0
        self.tout = 0.0
        self.pin = 0.0
        self.pout = 0.0
        self.massflow = 0.0
        self.req_massflow = 0.0  # Required massflow (to achieve setpoint tout)
        self.pr_req_massflow = 0.0
        self.pr_act_massflow = 0.0
        self.pr_opt = 0.0
        self.qabs = 0.0
        self.qlost_brackets = 0.0

        self.act_tin = 0.0
        self.act_tout = 0.0
        self.act_pin = 0.0
        self.act_pout = 0.0
        self.act_massflow = 0.0  # Actual massflow measured by the flowmeter

        self.wasted_power = 0.0
        self.tracking = True

    def initialize(self, source, values = None):

        if source == 'rated':
            self.massflow = self.parameters['rated_massflow']
            self.tin = self.parameters['rated_tin']
            self.pin = self.parameters['rated_pin']
            # self.tout = self.parameters['rated_tout']
            # self.pout = self.parameters['rated_pout']

        elif source == 'actual':
            self.massflow = self.act_massflow
            self.tin = self.act_tin
            self.pin = self.act_pin
            # self.tout = self.act_tout
            # self.pout = self.act_pout

        elif source == 'values' and values is not None:
            self.massflow = values['massflow']
            self.tin = values['tin']
            self.pin = values['pin']

        else:
            print('Select source [rated|actual|values]')
            sys.exit()


    def load_actual(self, subfield = None):

        if subfield == None:
            subfield = self.subfield

        self.act_massflow = subfield.act_massflow / len(subfield.loops)
        self.act_tin = subfield.act_tin
        self.act_tout = subfield.act_tout
        self.act_pin = subfield.act_pin
        self.act_pout = subfield.act_pout

    def set_loop_values_from_HCEs(self, type_of_massflow):

        pr_list = []
        qabs_list = []
        qlost_brackets_list = []
        qperd_list = []
        pr_opt_list = []

        for s in self.scas:
            for h in s.hces:
                pr_list.append(h.pr)
                qabs_list.append(h.qabs)
                qlost_brackets_list.append(h.qlost_brackets)
                qperd_list.append(h.qperd)
                pr_opt_list.append(h.pr_opt)

        if type_of_massflow == 'required':
            self.pr_req_massflow = np.mean(pr_list)
        elif type_of_massflow == 'actual':
            self.pr_act_massflow = np.mean(pr_list)

        self.qabs = np.sum(qabs_list)
        self.qlost_brackets = np.sum(qlost_brackets_list)
        self.qperd = np.sum(qperd_list)
        self.pr_opt = np.mean(pr_opt_list)


    def load_from_base_loop(self, base_loop):

        self.massflow = base_loop.massflow
        self.req_massflow = base_loop.req_massflow
        self.act_massflow = base_loop.act_massflow
        self.tin = base_loop.tin
        self.act_tin = base_loop.act_tin
        self.pin = base_loop.pin
        self.act_pin = base_loop.act_pin
        self.tout = base_loop.tout
        self.act_tout = base_loop.act_tout
        self.pout = base_loop.pout
        self.act_pout = base_loop.act_pout
        self.pr_req_massflow = base_loop.pr_req_massflow
        self.pr_act_massflow = base_loop.pr_act_massflow
        self.wasted_power = base_loop.wasted_power

    def calc_loop_pr_for_massflow(self, row, solarpos, htf, model):

        for s in self.scas:
            aoi = s.get_aoi(solarpos)
            for h in s.hces:
                h.set_pr_opt(aoi)
                h.set_qabs(aoi, solarpos, row)
                h.set_tin()
                h.set_pin()
                h.tout = h.tin
                model.calc_pr(h, htf, row)

        self.tout = self.scas[-1].hces[-1].tout
        self.pout = self.scas[-1].hces[-1].pout
        self.pr_req_massflow = 0.0
        self.set_loop_values_from_HCEs('actual')
        self.set_wasted_power(htf)

    def calc_loop_pr_for_tout(self, row, solarpos, htf, model):

        dri = self.scas[0].hces[0].parameters['Absorber tube inner diameter']
        min_reynolds = self.scas[0].hces[0].parameters['Min Reynolds']

        min_massflow = htf.get_massflow_from_Reynolds(dri, self.tin, self.pin,
                                                      min_reynolds)

        max_error = model.max_err_t  # % desviation tolerance
        search = True

        while search:

            self.calc_loop_pr_for_massflow(row, solarpos, htf, model)

            err = 100 * (abs(self.tout-self.parameters['rated_tout']) /
                   self.parameters['rated_tout'])

            if err > max_error:

                if self.tout >= self.parameters['rated_tout']:
                    self.massflow *= (1 + err / 100)
                    search = True
                elif (self.massflow > min_massflow and
                      self.massflow >
                      self.parameters['min_massflow']):
                    self.massflow *= (1 - err / 100)
                    search = True
                else:
                    self.massflow = max(
                        min_massflow,
                        self.parameters['min_massflow'])
                    self.calc_loop_pr_for_massflow(row, solarpos, htf, model)
                    search = False
            else:
                search = False

        self.tout = self.scas[-1].hces[-1].tout
        self.pout = self.scas[-1].hces[-1].pout
        self.req_massflow = self.massflow
        self.pr_act_massflow = 0.0
        self.set_loop_values_from_HCEs('required')


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

    def check_min_massflow(self, htf):

        dri = self.parameters['Absorber tube inner diameter']
        t = self.tin
        p = self.pin
        re = self.parameters['Min Reynolds']
        loop_min_massflow = htf.get_massflow_from_Reynolds(
                dri, t, p , re)

        if  self.massflow < loop_min_massflow:
            print("Too low massflow", self.massflow ,"<",
                  loop_min_massflow)

    def set_wasted_power(self, htf):

        if self.tout > self.parameters['tmax']:
            HH = htf.get_deltaH(self.tout, self.pout)
            HL = htf.get_deltaH(self.parameters['tmax'], self.pout)
            self.wasted_power = self.massflow * (HH - HL)
        else:
            self.wasted_power = 0.0

    def get_values(self, type = None):

        _values = {}

        if type is None:
            _values =  {'tin': self.tin, 'tout': self.tout,
                        'pin': self.pin, 'pout': self.pout,
                        'mf': self.massflow}
        elif type =='required':
            _values =  {'tin': self.tin, 'tout': self.tout,
                        'pin': self.pin, 'pout': self.pout,
                        'req_mf': self.req_massflow}
        elif type =='actual':
            _values =  {'actual_tin': self.actual_tin, 'tout': self.tout,
                        'pin': self.pin, 'pout': self.pout,
                        'mf': self.req_massflow}

        return _values

    def get_id(self):

        return 'LO.{0}.{1:000}'.format(self.subfield.name, self.loop_order)

class Loop(__Loop__):

    def __init__(self, subfield, loop_order, settings):

        self.subfield = subfield
        self.loop_order = loop_order

        super().__init__(settings)


    # def load_actual(self):

    #     self.act_massflow = self.subfield.act_massflow / len(self.subfield.loops)
    #     self.act_tin = self.subfield.act_tin
    #     self.act_pin = self.subfield.act_pin
    #     self.act_tout = self.subfield.act_tout
    #     self.act_pout = self.subfield.act_pout


class BaseLoop(__Loop__):

    def __init__(self, settings, sca_settings, hce_settings):

        super().__init__(settings)

        self.parameters_sca = sca_settings
        self.parameters_hce = hce_settings

        for s in range(settings['scas']):
            self.scas.append(SCA(self, s, sca_settings))
            for h in range(settings['hces']):
                self.scas[-1].hces.append(
                    HCE(self.scas[-1], h, hce_settings))


    def get_id(self, subfield = None):

        id = ''
        if subfield is not None:
            id = 'BL.'+subfield.name
        else:
            id = 'BL'

        return id

    def get_qlost_brackets(self, tf, text):

        #  Ec. 4.12

        L = self.parameters_hce['Length']
        n = (L / self.parameters_hce['Brackets'])
        pb = 0.2032
        acsb = 1.613e-4
        kb = 48
        hb = 20
        tbase = tf - 10

        return n * (np.sqrt(pb * kb * acsb * hb) * (tbase - text)) / L

    def get_pr_opt_peak(self):

        alpha = self.parameters_hce['Absorber absorptance']
        tau = self.parameters_hce['Envelope transmittance']
        rho = self.parameters_sca['Reflectance']
        gamma = self.get_solar_fraction()

        pr_opt_peak = alpha * tau * rho * gamma

        if pr_opt_peak > 1 or pr_opt_peak < 0:
            print("ERROR pr_opt_peak", pr_opt_peak)

        return pr_opt_peak


    def get_pr_geo(self, aoi):

        if aoi > 90:
            pr_geo = 0.0

        else:
            # Llamado "bordes" en Tesis. Pérdidas de los HCE de cabecera según aoi
            sca_unused_length = (self.parameters_sca["Focal Len"] *
                                 np.tan(np.radians(aoi)))

            pr_geo = 1 - sca_unused_length / (self.parameters['hces'] * \
                                              self.parameters_hce["Length"])

            if pr_geo > 1.0 or pr_geo < 0.0:
                print("ERROR pr_geo out of limits", pr_geo)

        return pr_geo


    def get_pr_shadows(self, solarpos):

        # Llamado "sombras" en Tesis. Pérdidas por sombras. ¿modelar sobre el SCA?
        # Sombras debidas a otros lazos

        shadowing = (1 -
                     np.sin(np.radians(abs(solarpos['elevation'][0]))) *
                     self.parameters['row_spacing'] /
                     self.parameters_sca['Aperture'])

        if shadowing < 0.0 or shadowing > 1.0:
            shadowing = 0.0

        if solarpos['elevation'][0] < 0:
            shadowing = 1.0

        pr_shadows = 1 - shadowing

        if pr_shadows > 1 or pr_shadows < 0:
            print("ERROR",  pr_shadows)

        return pr_shadows

    def get_solar_fraction(self):

        solarfraction = (self.parameters_sca['Reflectance'] *
                         self.parameters_sca['Geom.Accuracy'] *
                         self.parameters_sca['Track Twist'] *
                         self.parameters_sca['Cleanliness'] *
                         self.parameters_sca['Dust'] *
                         self.parameters_sca['Factor'] *
                         self.parameters_sca['Availability'])

        if solarfraction > 1 or solarfraction < 0:
            print("ERROR",  solarfraction)

        return solarfraction

    def get_IAM(self, theta):

        F0 = self.parameters_sca['IAM Coefficient F0']
        F1 = self.parameters_sca['IAM Coefficient F1']
        F2 = self.parameters_sca['IAM Coefficient F2']

        if (theta > 0 and theta < 80):
            kiam = (F0 + (F1 * np.radians(theta) + F2 * np.radians(theta)**2) /
                    np.cos(np.radians(theta)))
            if kiam > 1:
                kiam = 2 - kiam
            if kiam < 0:
                kiam = 0.0
        else:
            kiam = 0.0
        if kiam > 1 or kiam < 0:
            print("ERROR kiam={0:.4f} for angle of incidence {1:.4f}".format(
                kiam, theta))

        return kiam

    def get_aoi(self, solarpos):

        sigmabeta = 0.0
        beta0 = 0.0

        if self.parameters['Tracking Type'] == 1:  # N-S single axis tracker
            if solarpos['azimuth'][0] > 0 and solarpos['azimuth'][0] <= 180:
                surface_azimuth = 90  # Surface facing east
            else:
                surface_azimuth = 270  # Surface facing west
        elif self.parameters['Tracking Type'] == 2:  # E-W single axis tracker
            surface_azimuth = 180  # Surface facing the equator

        #  En esta fórmula asumo que el seguimiento del SCA es perfecto
        #  pero hay que ver la posibilidad de modelar cierto error o desfase
        beta0 = np.degrees(np.arctan(
            np.tan(np.radians(solarpos['zenith'][0])) *
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


class Subfield(object):
    '''
    Parabolic Trough Solar Field

    '''

    def __init__(self, solarfield, settings):

        self.solarfield = solarfield
        self.name = settings['name']
        self.parameters = settings
        self.loops = []

        self.tin = 0.0
        self.tout = 0.0
        self.pin = 0.0
        self.pout = 0.0
        self.massflow = 0.0
        self.req_massflow = 0.0
        self.pr_req_massflow = 0.0
        self.wasted_power = 0.0

        self.act_tin = 0.0
        self.act_tout = 0.0
        self.act_pin = 0.0
        self.act_pout = 0.0
        self.act_massflow = 0.0
        self.pr_act_massflow = 0.0

        self.rated_tin = self.solarfield.rated_tin
        self.rated_tout = self.solarfield.rated_tout
        self.rated_pin = self.solarfield.rated_pin
        self.rated_pout = self.solarfield.rated_pout
        self.rated_massflow = (self.solarfield.rated_massflow *
                               self.parameters['loops'] /
                               self.solarfield.total_loops)

    def set_massflow(self):

        mf = 0.0
        for l in self.loops:
            mf += l.massflow

        self.massflow = mf


    def set_req_massflow(self):

        req_mf = 0.0
        for l in self.loops:
            req_mf += l.req_massflow

        self.req_massflow = req_mf


    def set_wasted_power(self):

        wasted_power = 0.0
        for l in self.loops:
            wasted_power += l.wasted_power

        self.wasted_power = wasted_power

    def set_pr_req_massflow(self):

        loops_var = []

        for l in self.loops:

            loops_var.append(l.pr_req_massflow * l.req_massflow)

        self.pr_req_massflow = np.sum(loops_var) / self.req_massflow

    def set_pr_act_massflow(self):

        loops_var = []

        for l in self.loops:
            loops_var.append(l.pr_act_massflow * l.act_massflow)

        self.pr_act_massflow =  np.sum(loops_var) / self.act_massflow


    def set_pout(self):

        loops_var = []

        for l in self.loops:
            loops_var.append(l.pout * l.massflow)

        self.pout = np.sum(loops_var) / self.massflow

    def set_tout(self, htf):
        '''
        Calculates HTF output temperature throughout the solar field as a
        weighted average based on the enthalpy of the mass flow in each
        loop of the solar field
        '''
        # H = 0.0
        # H2 = 0.0
        dH = 0.0

        for l in self.loops:
            dH += l.massflow * htf.get_deltaH(l.tout, l.pout)

        dH /= self.massflow
        self.tout = htf.get_T(dH, self.pout)

    def loops_avg_out(self):

        tavg = 0.0
        cont = 0

        for l in self.loops:
            cont += 1
            tavg += l.tout

        return tavg / cont


    def initialize(self, source, values = None):

        if source == 'rated':
            self.massflow = self.rated_massflow
            self.tin = self.rated_tin
            self.pin = self.rated_pin
            self.tout = self.rated_tout
            self.pout = self.rated_pout

        elif source == 'actual':
            self.massflow = self.act_massflow
            self.tin = self.act_tin
            self.pin = self.act_pin
            self.tout = self.act_tout
            self.pout = self.act_pout

        elif source == 'values' and values is not None:
            self.massflow = values['massflow']
            self.tin = values['tin']
            self.pin = values['pin']

        else:
            print('Select source [rated|actual|values]')
            sys.exit()


    def load_actual(self, row):

        self.act_massflow = row[1][self.get_id() +'.act_mf']
        self.act_tin = row[1][self.get_id() +'.act_tin']
        self.act_pin = row[1][self.get_id() +'.act_pin']
        self.act_tout = row[1][self.get_id() +'.act_tout']
        self.act_pout = row[1][self.get_id() +'.act_pout']


    def get_id(self):

        return 'SB.' + self.name

class SolarField(object):
    '''
    Parabolic Trough Solar Field

    '''

    def __init__(self, subfield_settings, loop_settings, sca_settings, hce_settings):


        self.subfields = []
        self.total_loops = 0

        self.tin = 0.0
        self.tout = 0.0
        self.pin = 0.0
        self.pout = 0.0
        self.massflow = 0.0
        self.req_massflow = 0.0
        self.pr_req_massflow = 0.0
        self.wasted_power = 0.0

        self.act_tin = 0.0
        self.act_tout = 0.0
        self.act_pin = 0.0
        self.act_pout = 0.0
        self.act_massflow = 0.0
        self.pr_act_massflow = 0.0

        self.rated_tin = loop_settings['rated_tin']
        self.rated_tout = loop_settings['rated_tout']
        self.rated_pin = loop_settings['rated_pin']
        self.rated_pout = loop_settings['rated_pout']
        self.rated_massflow = (loop_settings['rated_massflow'] *
                               self.total_loops)


        for sf in subfield_settings:
            self.total_loops += sf['loops']
            self.subfields.append(Subfield(self, sf))
            for l in range(sf['loops']):
                self.subfields[-1].loops.append(
                    Loop(self.subfields[-1], l, loop_settings))
                for s in range(loop_settings['scas']):
                    self.subfields[-1].loops[-1].scas.append(
                        SCA(self.subfields[-1].loops[-1], s, sca_settings))
                    for h in range (loop_settings['hces']):
                        self.subfields[-1].loops[-1].scas[-1].hces.append(
                            HCE(self.subfields[-1].loops[-1].scas[-1], h,
                                hce_settings))

        # FUTURE WORK
        self.storage_available = False
        self.operation_mode = "subfield_heating"


    def initialize(self, source, values = None):

        if source == 'rated':
            self.massflow = self.rated_massflow
            self.tin = self.rated_tin
            self.pin = self.rated_pin
            self.tout = self.rated_tout
            self.pout = self.rated_pout

        elif source == 'actual':
            self.massflow = self.act_massflow
            self.tin = self.act_tin
            self.pin = self.act_pin
            self.tout = self.act_tout
            self.pout = self.act_pout

        elif source == 'values' and values is not None:
            self.massflow = values['massflow']
            self.tin = values['tin']
            self.pin = values['pin']

        else:
            print('Select source [rated|actual|values]')
            sys.exit()


    def load_actual(self, htf):

        massflow = 0.0
        H_tin = 0.0
        H_tout = 0.0
        list_pin = []
        list_pout = []

        for sf in self.subfields:
            massflow += sf.act_massflow
            H_tin += (sf.act_massflow *
                      htf.get_deltaH(sf.act_tin, sf.act_pin))
            H_tout += (sf.act_massflow *
                       htf.get_deltaH(sf.act_tout, sf.act_pout))
            list_pin.append(sf.act_pin * sf.act_massflow)
            list_pout.append(sf.act_pout * sf.act_massflow)

        H_tin /= massflow
        H_tout /= massflow

        self.act_massflow = massflow
        self.act_pin = np.sum(list_pin) / massflow
        self.act_pout = np.sum(list_pout) / massflow
        self.act_tin = htf.get_T(H_tin, self.act_pin)
        self.act_tout = htf.get_T(H_tout, self.act_pout)


    def set_massflow(self):

        mf = 0.0
        req_mf = 0.0

        for sf in self.subfields:
            mf += sf.massflow
            req_mf += sf.req_massflow

        self.massflow = mf
        self.req_massflow = req_mf


    def set_req_massflow(self):

        mf = 0.0

        for sf in self.subfields:
            mf += sf.req_massflow

        self.req_massflow = mf


    def set_wasted_power(self):

        wasted_power = 0.0
        for s in self.subfields:
            wasted_power += s.wasted_power

        self.wasted_power = wasted_power


    def set_pr_req_massflow(self):

        subfields_var = []
        for s in self.subfields:
            subfields_var.append(s.pr_req_massflow * s.req_massflow)

        self.pr_req_massflow = np.sum(subfields_var) / self.req_massflow


    def set_pr_act_massflow(self):

        subfields_var = []
        for s in self.subfields:
            subfields_var.append(s.pr_act_massflow * s.act_massflow)

        self.pr_act_massflow = np.sum(subfields_var) / self.act_massflow


    def set_pout(self):

        subfields_var = []

        for sf in self.subfields:
            subfields_var.append(sf.pout * sf.massflow)

        self.pout = np.sum(subfields_var) / self.massflow

    def set_tout(self, htf):
        '''
        Calculates HTF output temperature throughout the solar plant as a
        weighted average based on the enthalpy of the mass flow in each
        loop of the solar field
        '''
        dH = 0.0

        for sf in self.subfields:
            dH += sf.massflow * htf.get_deltaH(sf.tout, sf.pout)

        dH /= self.massflow
        self.tout = htf.get_T(dH, self.pout)

    def set_act_tout(self, htf):
        '''
        Calculates HTF output temperature throughout the solar plant as a
        weighted average based on the enthalpy of the mass flow in each
        loop of the solar field
        '''

        dH_actual = 0.0

        for sf in self.subfields:

            dH_actual += sf.act_massflow * htf.get_deltaH(sf.tout, sf.pout)
            dH_actual += sf.act_massflow * htf.get_deltaH(sf.act_tout, sf.act_pout)

        dH_actual /= self.act_massflow
        self.act_tout = htf.get_T(dH_actual, self.act_pout)


    def set_tin(self, htf):
        '''
        Calculates HTF output temperature throughout the solar plant as a
        weighted average based on the enthalpy of the mass flow in each
        loop of the solar field
        '''
        H = 0.0

        for sf in self.subfields:

            H += (htf.get_deltaH(sf.tin, sf.pin) * sf.massflow)

        H /= self.massflow
        self.tin = htf.get_T(H, self.rated_pin)


    def set_pin(self):

        subfields_var = []

        for sf in self.subfields:
            subfields_var.append(sf.pin * sf.massflow)

        self.pin = np.sum(subfields_var) / self.massflow


    def set_act_pin(self):

        subfields_var = []

        for sf in self.subfields:
            subfields_var.append(sf.act_pin * sf.act_massflow)

        self.act_pin = np.sum(subfields_var) / self.act_massflow


    def get_thermalpoweroutput(self, htf):

        HL = htf.get_deltaH(self.tin, self.pin)
        HH = htf.get_deltaH(self.tout, self.pout)

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


class SolarFieldSimulation(object):
    '''
    Definimos la clase simulacion para representar las diferentes
    pruebas que lancemos, variando el archivo TMY, la configuracion del
    site, la planta, el modo de operacion o el modelo empleado.
    '''

    def __init__(self, settings):

        self.ID =  settings['simulation']['ID']
        self.simulation = settings['simulation']['simulation']
        self.benchmark = settings['simulation']['benchmark']
        self.datatype = settings['simulation']['datatype']
        self.fastmode = settings['simulation']['fastmode']
        self.tracking = True
        self.solarfield = None
        self.powersystem = None
        self.htf = None
        self.coldfluid = None
        self.site = None
        self.datasource = None
        self.powercycle = None
        self.parameters = settings

        if settings['model']['name'] == 'Barbero4thOrder':
            self.model = ModelBarbero4thOrder(settings['model'])
        elif settings['model']['name'] == 'Barbero1stOrder':
            self.model = ModelBarbero1stOrder(settings['model'])
        elif settings['model']['name'] == 'BarberoSimplified':
            self.model = ModelBarberoSimplified(settings['model'])

        if self.datatype == 1:
            self.datasource = Weather(settings['simulation'])
        elif self.datatype == 2:
            self.datasource = FieldData(settings['simulation'],
                                        settings['tags'])

        if not hasattr(self.datasource, 'site'):
            self.site = Site(settings['site'])
        else:
            self.site = Site(self.datasource.site_to_dict())


        if settings['HTF']['source'] == "CoolProp":
            if settings['HTF']['CoolPropID'] not in Fluid._COOLPROP_FLUIDS:
                print("Not CoolPropID valid")
                sys.exit()
            else:
                self.htf = FluidCoolProp(settings['HTF'])

        else:
            self.htf = FluidTabular(settings['HTF'])

        self.solarfield = SolarField(settings['subfields'],
                                   settings['loop'],
                                   settings['SCA'],
                                   settings['HCE'])

        self.base_loop = BaseLoop(settings['loop'],
                                  settings['SCA'],
                                  settings['HCE'])


    def runSimulation(self):

        self.show_message()

        for row in self.datasource.dataframe.iterrows():

            solarpos = self.site.get_solarposition(row)

            self.gather_general_data(row, solarpos)

            if solarpos['zenith'][0] < 90:
                self.tracking = True
            else:
                self.tracking = False

            if self.simulation:
                self.simulate_solarfield(solarpos, row)
                self.gather_simulation_data(row)

            if self.benchmark and self.datatype == 2:  # 2: Field Data File available
                self.benchmarksolarfield(solarpos, row)
                self.gather_benchmark_data(row)

            self.plantperformance(row)
            str_data = ("{0} Ang. Zenith: {1:.2f} DNI: {2} W/m2 " +
                         "Qm: {3:.1f}kg/s Tin: {4:.1f}K Tout: {5:1f}K act_Tout: {6:1f}K")

            print(str_data.format(row[0], solarpos['zenith'][0],
                                   row[1]['DNI'], self.solarfield.massflow,
                                   self.solarfield.tin, self.solarfield.tout,
                                   self.solarfield.act_tout))

        self.save_results()


    def simulate_solarfield(self, solarpos, row):

        flag_0 = datetime.now()
        self.base_loop.initialize('rated')

        if self.datatype == 2:
            for s in self.solarfield.subfields:
                s.load_actual(row)
                s.initialize('actual')
            self.solarfield.load_actual(self.htf)
            self.base_loop.tin = self.solarfield.act_tin

        # Force minimum massflow at night
        if solarpos['zenith'][0] > 90:
            self.base_loop.massflow = self.base_loop.parameters['min_massflow']
            self.base_loop.req_massflow = self.base_loop.massflow
            self.base_loop.calc_loop_pr_for_massflow(
                row, solarpos, self.htf, self.model)
        else:
            self.base_loop.calc_loop_pr_for_tout(
                row, solarpos, self.htf, self.model)

        if  self.fastmode:

            for s in self.solarfield.subfields:
                s.initialize('rated')
                for l in s.loops:
                    l.load_from_base_loop(self.base_loop)

                s.set_massflow()
                s.set_req_massflow()
                s.set_pr_req_massflow()
                s.set_pout()
                s.set_tout(self.htf)
                s.set_wasted_power()

            flag_1 = datetime.now()
            delta_1 = flag_1 - flag_0

        else:
            for s in self.solarfield.subfields:
                for l in s.loops:
                    l.initialize('rated')
                    if solarpos['zenith'][0] > 90:
                        l.massflow = self.base_loop.parameters['min_massflow']
                        l.calc_loop_pr_for_massflow(
                            row, solarpos, self.htf, self.model)
                        l.req_massflow = l.massflow

                    else:
                        # Start with the previous loop massflow, for a better convergence
                        if l.loop_order > 0:
                            l.massflow = \
                                l.subfield.loops[l.loop_order-1].massflow
                        else:
                            l.massflow = self.base_loop.massflow

                        l.calc_loop_pr_for_tout(
                            row, solarpos, self.htf, self.model)

                    s.set_massflow()
                    s.set_req_massflow()
                    s.set_pr_req_massflow()
                    s.set_tout(self.htf)
                    s.set_pout()
                    s.set_wasted_power()

            flag_2 = datetime.now()
            delta_2 = flag_2 - flag_0

        self.solarfield.set_massflow()
        self.solarfield.set_req_massflow()
        self.solarfield.set_pr_req_massflow()
        self.solarfield.set_tin(self.htf)
        self.solarfield.set_pin()
        self.solarfield.set_tout(self.htf)
        self.solarfield.set_pout()


    def benchmarksolarfield(self, solarpos, row):

        flag_0 = datetime.now()

        if self.fastmode:
            for s in self.solarfield.subfields:
                s.load_actual(row)
                s.initialize('actual')
                self.base_loop.load_actual(s)
                self.base_loop.initialize('actual')
                self.base_loop.calc_loop_pr_for_massflow(
                    row, solarpos, self.htf, self.model)
                self.base_loop.set_loop_values_from_HCEs('actual')

                for l in s.loops:
                    l.load_from_base_loop(self.base_loop)

                s.set_massflow()
                s.set_pr_act_massflow()
                s.set_pout()
                s.set_tout(self.htf)
                s.set_wasted_power()

            flag_1 = datetime.now()
            delta_1 = flag_1 - flag_0

        else:
            for s in self.solarfield.subfields:
                s.load_actual(row)
                s.initialize('actual')

                for l in s.loops:
                    l.load_actual()
                    l.initialize('actual')
                    l.calc_loop_pr_for_massflow(
                        row, solarpos, self.htf, self.model)
                    l.set_loop_values_from_HCEs('actual')

                s.set_massflow()
                s.set_pr_act_massflow()
                s.set_pout()
                s.set_tout(self.htf)
                s.set_wasted_power()

                values = {self.base_loop.get_id(s) +'.tin': self.base_loop.tin,
                   self.base_loop.get_id(s) +'.tout': self.base_loop.tout,
                   self.base_loop.get_id(s) +'.pin': self.base_loop.pin,
                   self.base_loop.get_id(s) +'.pout': self.base_loop.pout,
                   self.base_loop.get_id(s) +'.act_mf': self.base_loop.act_massflow,
                   self.base_loop.get_id(s) +'.act_pr': self.base_loop.pr_act_massflow,
                   self.base_loop.get_id(s) +'.wasted_power': self.base_loop.wasted_power,
                   s.get_id() + '.tout': s.tout,
                   s.get_id() + '.wasted_power': s.wasted_power,
                   s.get_id() + '.act_pr': s.pr_act_massflow}

                self.store_values(row, values)

            flag_2 = datetime.now()
            delta_2 = flag_2 - flag_0

        self.solarfield.load_actual(self.htf)
        self.solarfield.initialize('actual')
        self.solarfield.set_massflow()
        self.solarfield.set_tout(self.htf)
        self.solarfield.set_pout()
        self.solarfield.set_pr_act_massflow()
        self.solarfield.set_wasted_power()


    def plantperformance(self, row):

        self.solarfield.set_operation_mode()

        if self.solarfield.operation_mode == "subfield_not_heating":
            massflow_to_HE = 0
        else:
            massflow_to_HE = self.solarfield.massflow

        # self.heatexchanger.set_htf_in(massflow_to_HE,
        #                                  self.solarfield.tout,
        #                                  self.solarfield.pout)

        # self.heatexchanger.set_coldfluid_in(
        #         self.powercycle.massflow, self.powercycle.tout, self.powercycle.pin) #PENDIENTE

        # self.heatexchanger.set_thermalpowertransfered()
        # self.powercycle.calc_pr_NCA(self.powercycle.tout) #Provisonal temperatura agua coondensador
        # self.generator.calc_pr()

    def store_values(self, row, values):

        for v in values:
            self.datasource.dataframe.at[row[0], v] = values[v]

    def gather_general_data(self, row, solarpos):

        #TO-DO: CALCULOS PARA AGREGAR AL DATAFRAME
        #  Solar postion data
        self.datasource.dataframe.at[row[0], 'elevation'] = \
            solarpos['elevation'][0]
        self.datasource.dataframe.at[row[0], 'zenith'] = \
            solarpos['zenith'][0]
        self.datasource.dataframe.at[row[0], 'azimuth'] = \
            solarpos['azimuth'][0]
        aoi = self.base_loop.get_aoi(solarpos)
        self.datasource.dataframe.at[row[0], 'aoi'] = aoi
        self.datasource.dataframe.at[row[0], 'IAM'] = \
            self.base_loop.get_IAM(aoi)
        self.datasource.dataframe.at[row[0], 'pr_geo'] = \
            self.base_loop.get_pr_geo(aoi)
        self.datasource.dataframe.at[row[0], 'pr_shadows'] = \
            self.base_loop.get_pr_shadows(solarpos)
        self.datasource.dataframe.at[row[0], 'pr_opt_peak'] = \
            self.base_loop.get_pr_opt_peak()
        self.datasource.dataframe.at[row[0], 'solar_fraction'] = \
            self.base_loop.get_solar_fraction()


    def gather_simulation_data(self, row):

        #  Solarfield data
        self.datasource.dataframe.at[row[0], 'REQ_MASSFLOW'] = \
            self.solarfield.req_massflow
        self.datasource.dataframe.at[row[0], 'TIN'] = \
            self.solarfield.tin
        self.datasource.dataframe.at[row[0], 'TOUT'] = \
            self.solarfield.tout
        self.datasource.dataframe.at[row[0], 'PIN'] = \
            self.solarfield.pin
        self.datasource.dataframe.at[row[0], 'POUT'] = \
            self.solarfield.pout
        self.datasource.dataframe.at[row[0], 'PR_REQ_MF'] = \
            self.solarfield.pr_req_massflow

        power_th = self.solarfield.req_massflow * \
            (self.htf.get_deltaH(self.solarfield.tout,
                            self.solarfield.pout) -
             self.htf.get_deltaH(self.solarfield.tin,
                            self.solarfield.pin))

        self.datasource.dataframe.at[row[0], 'THERMAL_POWER'] = power_th

        values = {
            self.base_loop.get_id() + '.mf': self.base_loop.req_massflow,
            self.base_loop.get_id() + '.tin': self.base_loop.tin,
            self.base_loop.get_id() + '.tout': self.base_loop.tout,
            self.base_loop.get_id() + '.pin': self.base_loop.pin,
            self.base_loop.get_id() + '.pout': self.base_loop.pout,
            self.base_loop.get_id() + '.req_pr': \
                self.base_loop.pr_req_massflow,
            self.base_loop.get_id() + '.pr_opt': \
                self.base_loop.pr_opt,
            self.base_loop.get_id() + '.qabs': self.base_loop.qabs,
            self.base_loop.get_id() + '.qlost_brackets': \
                self.base_loop.qlost_brackets,
            self.base_loop.get_id() + '.qperd': self.base_loop.qperd}

        self.store_values(row, values)

        df = self.datasource.dataframe

        for s in self.solarfield.subfields:
            #  Agretate data from subfields
            values = {
                s.get_id() + '.req_mf': s.req_massflow,
                s.get_id() + '.tin': s.tin,
                s.get_id() + '.tout': s.tout,
                s.get_id() + '.pin': s.pin,
                s.get_id() + '.pout': s.pout,
                s.get_id() + '.req_pr': s.pr_req_massflow}

            self.store_values(row, values)

            for l in s.loops:
                #  Loop data
                values = {
                    l.get_id() + '.req_mf':  l.req_massflow,
                    l.get_id() + '.tin':  l.tin,
                    l.get_id() + '.tout':  l.tout,
                    l.get_id() + '.pin':  l.pin,
                    l.get_id() + '.pout':  l.pout,
                    l.get_id() + '.pr_req_mf':  l.pr_req_massflow,
                    l.get_id() + '.pr_opt':  l.pr_opt}

                self.store_values(row, values)

    def gather_benchmark_data(self, row):

        #  Solarfield data
        self.datasource.dataframe.at[row[0], 'ACT_MASSFLOW'] = \
            self.solarfield.act_massflow
        self.datasource.dataframe.at[row[0], 'ACT_TIN'] = \
            self.solarfield.act_tin
        self.datasource.dataframe.at[row[0], 'ACT_TOUT'] = \
            self.solarfield.act_tout
        self.datasource.dataframe.at[row[0], 'TOUT'] = \
            self.solarfield.tout
        self.datasource.dataframe.at[row[0], 'ACT_PIN'] = \
            self.solarfield.act_pin
        self.datasource.dataframe.at[row[0], 'ACT_POUT'] = \
            self.solarfield.act_pout
        self.datasource.dataframe.at[row[0], 'POUT'] = \
            self.solarfield.pout
        self.datasource.dataframe.at[row[0], 'PR_ACT_MF'] = \
            self.solarfield.pr_act_massflow
        self.datasource.dataframe.at[row[0], 'WASTED_POWER'] = \
            self.solarfield.wasted_power

        act_power_th = self.solarfield.act_massflow * \
            (self.htf.get_deltaH(self.solarfield.act_tout,
                                 self.solarfield.act_pout) -
             self.htf.get_deltaH(self.solarfield.act_tin,
                                 self.solarfield.act_pin))

        available_power_th = self.solarfield.act_massflow * \
            (self.htf.get_deltaH(self.solarfield.tout,
                                 self.solarfield.pout) -
             self.htf.get_deltaH(self.solarfield.act_tin,
                                 self.solarfield.act_pin))

        self.datasource.dataframe.at[row[0], 'THERMAL_POWER'] = act_power_th
        self.datasource.dataframe.at[row[0], 'THERMAL_POWER_AVAILABLE'] = \
            available_power_th

        values = {
            self.base_loop.get_id() + '.mf': self.base_loop.act_massflow,
            self.base_loop.get_id() + '.tin': self.base_loop.act_tin,
            self.base_loop.get_id() + '.tout': self.base_loop.act_tout,
            self.base_loop.get_id() + '.pin': self.base_loop.act_pin,
            self.base_loop.get_id() + '.pout': self.base_loop.act_pout,
            self.base_loop.get_id() + '.act_pr':
                self.base_loop.pr_act_massflow,
            self.base_loop.get_id() + '.pr_opt':
                self.base_loop.pr_opt,
            self.base_loop.get_id() + '.qabs': self.base_loop.qabs,
            self.base_loop.get_id() + '.qlost_brackets':
                self.base_loop.qlost_brackets,
            self.base_loop.get_id() + '.qperd': self.base_loop.qperd}

        self.store_values(row, values)

        for s in self.solarfield.subfields:

            #  Agretate data from subfields
            values = {
                s.get_id() + '.act_mf': s.act_massflow,
                s.get_id() + '.tin': s.act_tin,
                s.get_id() + '.tout': s.act_tout,
                s.get_id() + '.pin': s.act_pin,
                s.get_id() + '.pout': s.act_pout,
                s.get_id() + '.act_pr': s.pr_act_massflow,
                s.get_id() + '.wasted_power': s.wasted_power}

            self.store_values(row, values)

            for l in s.loops:
                #  Loop data
                values = {
                    l.get_id() + '.act_mf':  l.act_massflow,
                    l.get_id() + '.tin':  l.act_tin,
                    l.get_id() + '.tout':  l.act_tout,
                    l.get_id() + '.pin':  l.act_pin,
                    l.get_id() + '.pout':  l.act_pout,
                    l.get_id() + '.pr_act_mf':  l.pr_act_massflow,
                    l.get_id() + '.pr_opt': l.pr_opt,
                    l.get_id() + '.wasted_power': l.wasted_power}

                self.store_values(row, values)

    def show_results(self):

        mf_keys = ['SF.req_mf', 'SF.act_mf']

        for c in self.datasource.dataframe.columns:
            if 'BaseLoop.req_mf' in c:
                mf_keys.append(c)

        for c in self.datasource.dataframe.columns:
            if 'BaseLoop.act_mf' in c:
                mf_keys.append(c)

        self.datasource.dataframe[mf_keys].plot(
                        figsize=(20,10), linewidth=5, fontsize=20)
        plt.xlabel('Date', fontsize=20)
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)

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


    def testgeo(self):


        for row in self.datasource.dataframe.iterrows():

            solarpos = self.site.get_solarposition(row)

            for s in self.base_loop.scas:
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


    def show_message(self):

        print("Running simulation for source data file: {0}".format(
            self.parameters['simulation']['filename']))
        print("Simulation: {0} ; Benchmark: {1} ; FastMode: {2}".format(
            self.parameters['simulation']['simulation'],
            self.parameters['simulation']['benchmark'],
            self.parameters['simulation']['fastmode']))

        print("Site: {0} @ Lat: {1:.2f}º, Long: {2:.2f}º, Alt: {3} m".format(
            self.site.name, self.site.latitude,
            self.site.longitude, self.site.altitude))

        print("Loops:", self.solarfield.total_loops,
              'SCA/loop:', self.parameters['loop']['scas'],
              'HCE/SCA:', self.parameters['loop']['hces'])
        print("SCA model:", self.parameters['SCA']['Name'])
        print("HCE model:", self.parameters['HCE']['Name'])
        if self.parameters['HTF']['source'] == 'table':
            print("HTF form table:", self.parameters['HTF']['name'])
        elif self.parameters['HTF']['source'] == 'CoolProp':
            print("HTF form CoolProp:", self.parameters['HTF']['CoolPropID'])
        print("-------------------------------------------------------------")


class LoopSimulation(object):
    '''
    Definimos la clase simulacion para representar las diferentes
    pruebas que lancemos, variando el archivo TMY, la configuracion del
    site, la planta, el modo de operacion o el modelo empleado.
    '''

    def __init__(self, settings):


        self.tracking = True
        self.htf = None
        self.coldfluid = None
        self.site = None
        self.datasource = None
        self.parameters = settings

        if settings['model']['name'] == 'Barbero4thOrder':
            self.model = ModelBarbero4thOrder(settings['model'])
        elif settings['model']['name'] == 'Barbero1stOrder':
            self.model = ModelBarbero1stOrder(settings['model'])
        elif settings['model']['name'] == 'BarberoSimplified':
            self.model = ModelBarberoSimplified(settings['model'])

        self.datasource = TableData(settings['simulation'])

        self.site = Site(settings['site'])

        if settings['HTF']['source'] == "CoolProp":
            if settings['HTF']['CoolPropID'] not in Fluid._COOLPROP_FLUIDS:
                print("Not CoolPropID valid")
                sys.exit()
            else:
                self.htf = FluidCoolProp(settings['HTF'])
        else:
            self.htf = FluidTabular(settings['HTF'])

        self.base_loop = BaseLoop(settings['loop'],
                                  settings['SCA'],
                                  settings['HCE'])

    def runSimulation(self):

        self.show_message()

        flag_0 = datetime.now()

        for row in self.datasource.dataframe.iterrows():

            solarpos = self.site.get_solarposition(row)

            if solarpos['zenith'][0] < 90:
                self.tracking = True
            else:
                self.tracking = False

            self.simulate_base_loop(solarpos, row)


            # self.gather_data(row, solarpos)
            str_data = ("{0} Ang. Zenith: {1:.2f} DNI: {2} W/m2 " +
                         "Qm: {3:.1f}kg/s Tin: {4:.1f}K Tout: {5:1f}K")

            print(str_data.format(row[0], solarpos['zenith'][0],
                                   row[1]['DNI'], self.base_loop.act_massflow,
                                   self.base_loop.tin, self.base_loop.tout))

        print(self.datasource.dataframe)

        self.datasource.dataframe[['Z','E','aoi']].plot(
            use_index = False,
            x = 'Z',
            y = 'aoi',
            figsize=(20,10), linewidth=5, fontsize=20)

        rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        ## for Palatino and other serif fonts use:
        ## rc('font',**{'family':'serif','serif':['Palatino']})
        rc('text', usetex=True)

        plt.ylabel('Performance, pr', fontsize=20)
        plt.xlabel('$q_{abs} [kW]$', fontsize=20)

        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)

        flag_1 = datetime.now()
        delta_t = flag_1 - flag_0
        print("Total runtime: ", delta_t.total_seconds())

        self.save_results()


    def simulate_base_loop(self, solarpos, row):

        values = {'tin': 573,
                  'pin': 1900000,
                  'massflow': 4}
        self.base_loop.initialize('values', values)
        HCE_var = ''
        SCA_var = ''

        for c in row[1].keys():
            if c in self.base_loop.parameters_sca.keys():
                SCA_var = c

        for c in row[1].keys():
            if c in self.base_loop.parameters_hce.keys():
                HCE_var = c

        for s in self.base_loop.scas:
            if SCA_var != '':
                s.parameters[SCA_var] = row[1][SCA_var]
            aoi = s.get_aoi(solarpos)
            for h in s.hces:
                if HCE_var != '':
                    h.parameters[HCE_var] = row[1][HCE_var]
                h.set_pr_opt(aoi)
                h.set_qabs(aoi, solarpos, row)
                h.set_tin()
                h.set_pin()
                h.tout = h.tin
                self.model.calc_pr(h, self.htf, row)

        self.base_loop.tout = self.base_loop.scas[-1].hces[-1].tout
        self.base_loop.pout = self.base_loop.scas[-1].hces[-1].pout
        self.base_loop.set_loop_values_from_HCEs('actual')
        print('pr', self.base_loop.pr_act_massflow ,
              'tout', self.base_loop.tout,
              'massflow', self.base_loop.massflow)

        if HCE_var +  SCA_var != '':
            self.datasource.dataframe.at[row[0], HCE_var + SCA_var] = row[1][HCE_var + SCA_var]
        self.datasource.dataframe.at[row[0], 'pr'] = self.base_loop.pr_act_massflow
        self.datasource.dataframe.at[row[0], 'tout'] = self.base_loop.tout
        self.datasource.dataframe.at[row[0], 'pout'] = self.base_loop.pout
        self.datasource.dataframe.at[row[0], 'Z'] = solarpos['zenith'][0]
        self.datasource.dataframe.at[row[0], 'E'] = solarpos['elevation'][0]
        self.datasource.dataframe.at[row[0], 'aoi'] = aoi



    def save_results(self):


        try:
            initialdir = "./simulations_outputs/"
            prefix = datetime.today().strftime("%Y%m%d %H%M%S")
            filename = "Loop Simulation"
            sufix = ".csv"

            path = initialdir + prefix + filename + sufix

            self.datasource.dataframe.to_csv(path, sep=';', decimal = ',')

        except Exception:
            raise
            print('Error saving results, unable to save file: %r', path)


    def show_message(self):

        print("Running simulation for source data file: {0}".format(
            self.parameters['simulation']['filename']))

        print("Site: {0} @ Lat: {1:.2f}º, Long: {2:.2f}º, Alt: {3} m".format(
            self.site.name, self.site.latitude,
            self.site.longitude, self.site.altitude))

        print('SCA/loop:', self.parameters['loop']['scas'],
              'HCE/SCA:', self.parameters['loop']['hces'])
        print("SCA model:", self.parameters['SCA']['Name'])
        print("HCE model:", self.parameters['HCE']['Name'])
        print("HTF:", self.parameters['HTF']['name'])
        print("---------------------------------------------------")


class Air(object):


    def get_cinematic_viscosity(self, t):

        # t: air temperature [K]
        return 8.678862e-11 * t**2 + 4.069284e-08 * t + -4.288741e-06


    def get_reynolds(self, t, L, v):

        nu = self.get_cinematic_viscosity(t)

        return v * L / nu

class Fluid:

    _T_REF = 285.856
    _COOLPROP_FLUIDS = ['Water', 'INCOMP::TVP1', 'INCOMP::S800']

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

    def get_cp(self, p, t):
        pass

    def get_density(self, p, t):
        pass

    def get_thermal_conductivity(self, p, t):
        pass

    def get_deltaH(self, p, t):
        pass

    def get_T(self, h, p):
        pass

    def get_dynamic_viscosity(self, t, p):
        pass

    def get_Reynolds(self, dri, t, p, massflow):

        if t > self.tmax:
            t= self.tmax

        return (4 * massflow /
                (np.pi * dri * self.get_dynamic_viscosity(t,p)))

    def get_massflow_from_Reynolds(self, dri, t, p, re):

        if t > self.tmax:
            t = self.tmax

        return re * np.pi * dri * self.get_dynamic_viscosity(t,p) / 4

    def get_prandtl(self, t, p):

        #  Specific heat capacity
        cp = self.get_cp(t, p)

        #  Fluid dynamic viscosity
        mu = self.get_dynamic_viscosity(t, p)

        #  Fluid density
        rho = self.get_density(t, p)

        #  Fluid thermal conductivity
        kf = self.get_thermal_conductivity(t, p)

        #  Fluid thermal diffusivity
        alpha = kf / (rho * cp)

        #  Prandtl number
        prf = mu / alpha

        return prf

class FluidCoolProp(Fluid):

    def __init__(self, settings = None):

        if settings['source'] == 'table':
            self.name = settings['name']
            self.tmax = settings['tmax']
            self.tmin = settings['tmin']

        elif settings['source'] == 'CoolProp':
            self.tmax = PropsSI('T_MAX', settings['CoolPropID'])
            self.tmin = PropsSI('T_MIN', settings['CoolPropID'])
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


class FluidTabular(Fluid):

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

        href =(h0 + h1 *self._T_REF + h2 *self._T_REF**2 + h3 *self._T_REF**3 +
               h4 *self._T_REF**4 + h5 *self._T_REF**5)

        h = (h0 + h1 * t + h2 * t**2 + h3 * t**3 + h4 * t**4 + h5 * t**5)
        return (h - href)

    def get_T(self, h, p):

        t0, t1, t2, t3, t4, t5 = tuple(self.t)

        return (t0 + t1 * h + t2 * h**2 + t3 * h**3 +
                t4 * h**4 + t5 * h**5)


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

        return {"name": 'nombre_site',
                "latitude": self.site['latitude'],
                "longitude": self.site['longitude'],
                "altitude": self.site['altitude']}
        # return {"name": self.site['City'],
        #         "latitude": self.site['latitude'],
        #         "longitude": self.site['longitude'],
        #         "altitude": self.site['altitude']}

class FieldData(object):

    def __init__(self, settings, tags = None):
        self.filename = settings['filename']
        self.filepath = settings['filepath']
        self.file = self.filepath + self.filename
        self.tags = tags
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
                        self.dataframe = pd.read_csv(path, sep=';',
                                                     decimal= ',',
                                                     dayfirst=True,
                                                     index_col=0)
                        self.file = path
                    elif strext == ".xls":

                        self.dataframe = pd.read_excel(path)
                        self.file = path
                    else:
                        print("unknow extension ", strext)
                        return
            else:
                strfilename, strext = os.path.splitext(path)

                if  strext == ".csv":
                    self.dataframe = pd.read_csv(path, sep=';',
                                                 decimal= ',',
                                                 dayfirst=True,
                                                 index_col=0)
                    self.file = path
                elif strext == ".xls":
                    self.dataframe = pd.read_excel(path)
                    self.file = path
                else:
                    print("unknow extension ", strext)
                    return
        except Exception:
            raise
            txMessageBox.showerror('Error loading FieldData File',
                                   'Unable to open file: %r', self.file)

        self.dataframe.index = pd.to_datetime(self.dataframe.index,
                                              infer_datetime_format=True)

    def change_units(self):

        for c in self.dataframe.columns:
            if ('.act_t' in c) or ('DryBulb' in c) or ('Dew' in c):
                self.dataframe[c] += 273.15 # From Celsius Degrees to K
            if '.act_p' in c:
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

class TableData(object):

    def __init__(self, settings):
        self.filename = settings['filename']
        self.filepath = settings['filepath']
        self.file = self.filepath + self.filename
        self.dataframe = None

        self.openDataFile(self.file)


    def openDataFile(self, path = None):

        '''
        Table ddata
        '''

        try:
            if path is None:
                root = Tk()
                root.withdraw()
                path = askopenfilename(initialdir = ".data_files/",
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
                        self.dataframe = pd.read_csv(path, sep=';',
                                                     decimal= ',',
                                                     dayfirst=True)
                        self.file = path
                    elif strext == ".xls":

                        self.dataframe = pd.read_excel(path)
                        self.file = path
                    else:
                        print("unknow extension ", strext)
                        return
            else:
                strfilename, strext = os.path.splitext(path)

                if  strext == ".csv":
                    self.dataframe = pd.read_csv(path, sep=';',
                                                 decimal= ',',
                                                 dayfirst=True)
                    self.file = path
                elif strext == ".xls":
                    self.dataframe = pd.read_excel(path)
                    self.file = path
                else:
                    print("unknow extension ", strext)
                    return
        except Exception:
            raise
            txMessageBox.showerror('Error loading FieldData File',
                                   'Unable to open file: %r', self.file)


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

    def get_rated_massflow(self, ratedpower, solarfield, htf, coldfluid = None):

        HH = htf.get_deltaH(solarfield.rated_tout, solarfield.rated_pout)
        HL = htf.get_deltaH(solarfield.rated_tin, solarfield.rated_pin)
        rated_massflow = self.get_thermalpower(self.ratedpower) / (HH - HL)

        return rated_massflow

class HeatExchanger(object):

    def __init__(self, settings, htf, coldfluid):

        self.pr = settings['pr']
        self.thermalpowertransfered = 0.0

        self.htfmassflow = 0.0
        self.htftin = 0.0
        self.htftout = 393 + 273.15
        self.htfpin = 1500000
        self.htfpout = 0.0
        self.htf = htf

        self.coldfluidmassflow = 0.0
        self.coldfluidtin = 0.0
        self.coldfluidtout = 0.0
        self.colfluidpin = 0.0
        self.colfluidtout = 0.0
        self.coldfluid = coldfluid

    def set_htf_in(self, htfmassflow, htftin, htfpin):

        self.htfmassflow = htfmassflow
        self.htftin = htftin
        self.htfpin = htfpin

    def set_coldfluid_in(self, coldfluidmassflow, coldfluidtin, coldfluidpin):

        self.coldfluidmassflow = coldfluidmassflow
        self.coldfluidtin = coldfluidtin
        self.coldfluidpin = coldfluidpin

    def set_thermalpowertransfered(self):

        #provisional
        pinch = 10.0

        HH = self.htf.get_deltaH(self.htftin, self.htfpin)
        HL = self.htf.get_deltaH(self.coldfluidtin + pinch, self.htfpout)

        self.thermalpowertransfered = self.pr * ((HH-HL)*self.htfmassflow)
#        self.thermalpowertransfered = (
#                self.pr *
#                self.htfmassflow *
#                self.htf.get_cp(self.htftin, self.htfpin) *
#                self.htftin - hftout)


    def htf_tout():
        return


class HeatStorage(object):

    def __init__(self, settings):

        self.name = settings['heatstorage']

    def set_fluid_tout(self):
        pass

    def htf_tout():
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


    def get_solarposition(self, row):

        solarpos = pvlib.solarposition.get_solarposition(
                row[0],
                self.latitude,
                self.longitude,
                self.altitude,
                pressure=row[1]['Pressure'],
                temperature=row[1]['DryBulb'])

        return solarpos


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








