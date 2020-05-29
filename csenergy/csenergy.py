# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 08:26:43 2019
@author: pacomunuera
"""

import numpy as np
import scipy as sc
from scipy import constants
import math as mt
from decimal import Decimal
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
from datetime import datetime, timedelta
import time
import os.path
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
import pytz
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

    def calc_pr(self, hce, htf, row, qabs = None):

        if qabs is None:
            qabs = hce.qabs

        tin = hce.tin

        #  If the hce is the first one in the loop tf = tin, else
        #  tf equals tin plus half the jump of temperature in the previous hce

        if hce.hce_order == 0:
            tf = hce.tin  # HTF bulk temperature
        else:
            tf = hce.tin + 0.5 * (hce.sca.hces[hce.hce_order - 1].tout -
                                  hce.sca.hces[hce.hce_order - 1].tin)

        tri = tf  #  Absorber tube inner surface temperature
        massflow = hce.sca.loop.massflow
        wspd = row[1]['Wspd']  #  Wind speed
        text = row[1]['DryBulb']  #  Dry bulb ambient temperature
        sigma = sc.constants.sigma  #  Stefan-Bolztmann constant
        dro = hce.parameters['Absorber tube outer diameter']
        dri = hce.parameters['Absorber tube inner diameter']
        dgo = hce.parameters['Glass envelope outer diameter']
        dgi = hce.parameters['Glass envelope inner diameter']
        L = (hce.parameters['Length'] * hce.parameters['Bellows ratio'] *
            hce.parameters['Shield shading'])
        A = hce.sca.parameters['Aperture']
        x = 1 #  Calculation grid fits hce longitude

        #  HCE wall thermal conductivity
        # krec = hce.get_krec(tf)

        #  Specific Capacity
        cp = htf.get_cp(tf, hce.pin)

        #  Internal transmission coefficient.
        # hint = hce.get_hint(tf, hce.pin, htf)

        #  Internal heat transfer coefficient. Eq. 3.22 Barbero2016
        # urec = 1 / ((1 / hint) + (dro * np.log(dro / dri)) / (2 * krec))
        urec = hce.get_urec(tf, hce.pin, htf)
        #  We suppose thermal performance, pr = 1, at first if  the hce is
        #  the first one in the loop or pr_j = pr_j-1 if there is a previous
        #  HCE in the loop.
        pr = hce.get_previous_pr()
        tro = tf + qabs * pr / urec

        #  HCE emittance
        eext = hce.get_eext(tro, wspd)
        #  External Convective Heat Transfer equivalent coefficient
        hext = hce.get_hext(wspd)

        #  Thermal power lost througth  bracktets
        qlost_brackets = hce.get_qlost_brackets(tro, text)

        #  Thermal power lost. Eq. 3.23 Barbero2016
        qlost = sigma * eext * (tro**4 - text**4) + hext * (tro - text) + \
            qlost_brackets


        #  Critical Thermal power loss. Eq. 3.50 Barbero2016
        qcrit = sigma * eext * (tf**4 - text**4) + hext * (tf - text)

        #  Critical Internal heat transfer coefficient, Eq. 3.51 Barbero2016
        ucrit = 4 * sigma * eext * tf**3 + hext

        #  Transmission Units Number, Ec. 3.30 Barbero2016
        NTU = urec * x * L * np.pi * dro / (massflow * cp)

        #  if qabs > qcrit:
        if qabs > 1.1 * qcrit:

            #  We use Barbero2016's simplified model aproximation
            #  Eq. 3.63 Barbero2016
            fcrit = 1 / (1 + (ucrit / urec))

            #  Eq. 3.71 Barbero2016
            pr = fcrit * (1 - qcrit / qabs)

            errtro = 10.
            errpr = 1.
            step = 0

            while ((errtro > self.max_err_tro or errpr > self.max_err_pr) and
                   step < 1000):

                step += 1

                #  Eq. 3.32 Barbero2016
                f0 = qabs / (urec * (tf - text))

                #  Eq. 3.34 Barbero2016
                f1 = ((4 * sigma * eext * text**3) + hext) / urec
                f2 = 6 * (text**2) * (sigma * eext / urec) * (qabs / urec)
                f3 = 4 * text * (sigma * eext / urec) * ((qabs / urec)**2)
                f4 = (sigma * eext / urec) * ((qabs / urec)**3)

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
                # krec = hce.get_krec(tf)

                #  Specific Capacity
                cp = htf.get_cp(tf, hce.pin)

                #  Internal transmission coefficient.
                # hint = hce.get_hint(tf, hce.pin, htf)

                #  Internal heat transfer coefficient. Eq. 3.22 Barbero2016
                # urec = 1 / ((1 / hint) + (dro * np.log(dro / dri)) / (2 * krec))
                urec = hce.get_urec(tf, hce.pin, htf)

                #  HCE emittance
                eext = hce.get_emittance(tro, wspd)

                #  External Convective Heat Transfer equivalent coefficient
                hext = hce.get_hext(wspd)

                #  We calculate tro again.

                tro2 = tf + qabs * pr / urec
                errtro = abs(tro2-tro)
                tro = tro2

                #  Increase qlost with the thermal power lost througth bracktets
                qlost_brackets = hce.get_qlost_brackets(tro, text)

                #  Thermal power loss. Eq. 3.23 Barbero2016
                qlost = sigma * eext * (tro**4 - text**4) + \
                    hext * (tro - text) + qlost_brackets

                #  Critical Thermal power loss. Eq. 3.50 Barbero2016
                qcrit = sigma * eext * (tf**4 - text**4) + hext * (tf - text)

                #  Critical Internal heat transfer coefficient, Eq. 3.51 Barbero2016
                ucrit = 4 * sigma * eext * tf**3 + hext

                #  Transmission Units Number, Ec. 3.30 Barbero2016
                NTU = urec * x * L * np.pi * dro / (massflow * cp)


            if step == 1000:
                print('No se alcanzó convergencia. HCE', hce.get_index())
                print(qabs, qcrit, urec, ucrit)
            # hce.pr = hce.pr * (1 - qlost_brackets / qabs)
            hce.qlost = qlost
            hce.qlost_brackets = qlost_brackets

        else:
            errtro = 10.0
            while (errtro > self.max_err_tro):

                hce.pr = 0.0
                #  tro = tf
                tro = 0.5 * (hce.tin + hce.tout)
                #  HCE emittance
                eext = hce.get_emittance(tro, wspd)
                #  External Convective Heat Transfer equivalent coefficient
                hext = hce.get_hext(wspd)

                #  Thermal power lost. Eq. 3.23 Barbero2016
                qlost = sigma * eext * (tro**4 - text**4) + \
                    hext * (tro - text) + hce.get_qlost_brackets(tro, text)

                #  Thermal power lost througth  bracktets
                qlost_brackets = hce.get_qlost_brackets(tro, text)

                hce.qlost = qlost
                hce.qlost_brackets =  qlost_brackets
                hce.set_tout(htf)
                hce.set_pout(htf)
                errtro = abs(tro - 0.5 * (hce.tin + hce.tout))

class ModelBarbero1stOrder(Model):

    def __init__(self, settings):

        self.parameters = settings
        self.max_err_t = self.parameters['max_err_t']

    def calc_pr(self, hce, htf, row, qabs = None):

        if qabs is None:
            qabs = hce.qabs

        tin = hce.tin

        #  If the hce is the first one in the loop tf = tin, else
        #  tf equals tin plus half the jump of temperature in the previous hce
        if hce.hce_order == 0:
            tf = hce.tin  # HTF bulk temperature
        else:
            tf = hce.tin + 0.5 * (hce.sca.hces[hce.hce_order - 1].tout -
                                  hce.sca.hces[hce.hce_order - 1].tin)
        tri = tf  #  Absorber tube inner surface temperature
        massflow = hce.sca.loop.massflow
        wspd = row[1]['Wspd']  #  Wind speed
        text = row[1]['DryBulb']  #  Dry bulb ambient temperature
        sigma = sc.constants.sigma  #  Stefan-Bolztmann constant
        dro = hce.parameters['Absorber tube outer diameter']
        dri = hce.parameters['Absorber tube inner diameter']
        dgo = hce.parameters['Glass envelope outer diameter']
        dgi = hce.parameters['Glass envelope inner diameter']
        L = (hce.parameters['Length'] * hce.parameters['Bellows ratio'] *
            hce.parameters['Shield shading'])
        A = hce.sca.parameters['Aperture']
        x = 1 #  Calculation grid fits hce longitude

        #  HCE wall thermal conductivity
        # krec = hce.get_krec(tf)

        #  Specific Capacity
        cp = htf.get_cp(tf, hce.pin)

        #  Internal transmission coefficient.
        # hint = hce.get_hint(tf, hce.pin, htf)

        #  Internal heat transfer coefficient. Eq. 3.22 Barbero2016
        # urec = 1 / ((1 / hint) + (dro * np.log(dro / dri)) / (2 * krec))
        urec = hce.get_urec(tf, hce.pin, htf)
        #  We suppose performance, pr = 1, at first
        pr = 1.0
        tro = tf + qabs * pr / urec

        #  HCE emittance
        eext = hce.get_eext(tro, wspd)
        #  External Convective Heat Transfer equivalent coefficient
        hext = hce.get_hext(wspd)

        #  Thermal power lost througth  bracktets
        qlost_brackets = hce.get_qlost_brackets(tf, text)

        #  Thermal power lost. Eq. 3.23 Barbero2016
        qlost = sigma * eext * (tro**4 - text**4) + hext * (tro - text) + \
            qlost_brackets

        #  Critical Thermal power loss. Eq. 3.50 Barbero2016
        qcrit = sigma * eext * (tf**4 - text**4) + hext * (tf - text) + \
            qlost_brackets

        #  Critical Internal heat transfer coefficient, Eq. 3.51 Barbero2016
        ucrit = 4 * sigma * eext * tf**3 + hext

        #  Ec. 3.63
        ## fcrit = (1 / ((4 * eext * tfe**3 / urec) + (hext / urec) + 1))
        fcrit = 1 / (1 + (ucrit / urec))

        #  Ec. 3.64
        Aext = np.pi * dro * x  # Pendiente de confirmar
        NTUperd = ucrit * Aext / (massflow * cp)

        if qabs > 1.1 * qcrit:
            hce.pr = ((1 - (qcrit / qabs)) *
                  (1 / (NTUperd * x)) *
                  (1 - np.exp(-NTUperd * fcrit * x)))
            # hce.pr = hce.pr * (1 - qlost_brackets / qabs)
        else:
            hce.pr = 0.0

        hce.qlost = qlost
        hce.qlost_brackets = qlost_brackets
        hce.set_tout(htf)
        hce.set_pout(htf)

class ModelBarberoSimplified(Model):

    def __init__(self, settings):

        self.parameters = settings
        self.max_err_t = self.parameters['max_err_t']

    def calc_pr(self, hce, htf, row, qabs=None):

        if qabs is None:
            qabs = hce.qabs

        tin = hce.tin
        #  If the hce is the first one in the loop tf = tin, else
        #  tf equals tin plus half the jump of temperature in the previous hce
        if hce.hce_order == 0:
            tf = hce.tin  # HTF bulk temperature
        else:
            tf = hce.tin + 0.5 * (hce.sca.hces[hce.hce_order - 1].tout -
                                  hce.sca.hces[hce.hce_order - 1].tin)


        tri = tf  #  Absorber tube inner surface temperature
        massflow = hce.sca.loop.massflow
        wspd = row[1]['Wspd']  #  Wind speed
        text = row[1]['DryBulb']  #  Dry bulb ambient temperature
        sigma = sc.constants.sigma  #  Stefan-Bolztmann constant
        dro = hce.parameters['Absorber tube outer diameter']
        dri = hce.parameters['Absorber tube inner diameter']
        dgo = hce.parameters['Glass envelope outer diameter']
        dgi = hce.parameters['Glass envelope inner diameter']
        L = (hce.parameters['Length'] * hce.parameters['Bellows ratio'] *
            hce.parameters['Shield shading'])
        A = hce.sca.parameters['Aperture']
        x = 1 #  Calculation grid fits hce longitude

        #  HCE wall thermal conductivity
        # krec = hce.get_krec(tf)

        #  Specific Capacity
        cp = htf.get_cp(tf, hce.pin)

        #  Internal transmission coefficient.
        # hint = hce.get_hint(tf, hce.pin, htf)

        #  Internal heat transfer coefficient. Eq. 3.22 Barbero2016
        # urec = 1 / ((1 / hint) + (dro * np.log(dro / dri)) / (2 * krec))
        urec = hce.get_urec(tf, hce.pin, htf)

        #  We suppose performance, pr = 1, at first
        pr = 1.0
        tro = tf + qabs * pr / urec

        #  HCE emittance
        eext = hce.get_eext(tro, wspd)
        #  External Convective Heat Transfer equivalent coefficient
        hext = hce.get_hext(wspd)

        #  Thermal power lost througth  bracktets
        qlost_brackets = hce.get_qlost_brackets(tf, text)

        #  Thermal power loss. Eq. 3.23 Barbero2016
        qlost = sigma * eext * (tro**4 - text**4) + hext * (tro - text) + \
            qlost_brackets

        #  Critical Thermal power loss. Eq. 3.50 Barbero2016
        qcrit = sigma * eext * (tf**4 - text**4) + hext * (tf - text)

        #  Critical Internal heat transfer coefficient, Eq. 3.51 Barbero2016
        ucrit = 4 * sigma * eext * tf**3 + hext

        #  Ec. 3.63
        ## fcrit = (1 / ((4 * eext * tfe**3 / urec) + (hext / urec) + 1))
        fcrit = 1 / (1 + (ucrit / urec))

        if qabs > 1.1 * qcrit:

            hce.pr = fcrit * (1 - (qcrit / qabs))

        else:
            hce.pr = 0

        hce.qlost = qlost
        hce.qlost_brackets = qlost_brackets
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
        self.parameters = settings  # Set of parameters of the HCE
        self.tin = 0.0  # Temperature of the HTF when enters in the HCE
        self.tout = 0.0  # Temperature of the HTF when goes out the HCE
        self.pin = 0.0  # Pressure of the HTF when enters in the HCE
        self.pout = 0.0  # Pressure of the HTF when goes out the HCE
        self.pr = 0.0  # Thermal performance of the HCE
        self.pr_opt = 0.0  # Optical performance of the HCE+SCA set
        self.qabs = 0.0  # Thermal which reach the aboserber tube
        self.qlost = 0.0  # Thermal power lost througth out the HCE
        self.qlost_brackets = 0.0  # Thermal power lost in brackets and arms

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

        if self.pr > 0:

            q = (self.qabs * np.pi * self.pr *
                 self.parameters['Length'] *
                 self.parameters['Bellows ratio'] *
                 self.parameters['Shield shading'] *
                 self.parameters['Absorber tube outer diameter'])

        else:
            q =  ( -1 * (self.qlost) * np.pi *
                  self.parameters['Length'] *
                  self.parameters['Bellows ratio'] *
                  self.parameters['Shield shading'] *
                  self.parameters['Absorber tube outer diameter'])

        self.tout = htf.get_T2(self.tin, q, self.sca.loop.massflow, self.pin)


    def set_pout(self, htf):

        # TO-DO CÁLCULLO PERDIDA DE CARGA:
        # Ec. Colebrook-White para el cálculo del factor de fricción de Darcy

        # re_turbulent = 2300

        # k = self.parameters['Inner surface roughness']
        # D = self.parameters['Absorber tube inner diameter']
        # re = htf.get_Reynolds(D, self.tin, self.pin, self.sca.loop.massflow)


        # if re < re_turbulent:
        #     darcy_factor = 64 / re

        # else:
        #     # a = (k /  D) / 3.7
        #     a = k / 3.7
        #     b = 2.51 / re
        #     x = 1

        #     fx = lambda x: x + 2 * np.log10(a + b * x )

        #     dfx = lambda x: 1 + (2 * b) / (np.log(10) * (a + b * x))

        #     root = sc.optimize.newton(fx,
        #                               x,
        #                               fprime=dfx,
        #                               maxiter=10000)

        #     darcy_factor = 1 / (root**2)

        # rho = htf.get_density(self.tin, self.pin)
        # v = 4 * self.sca.loop.massflow / (rho * np.pi * D**2)
        # g = sc.constants.g

        # # Ec. Darcy-Weisbach para el cálculo de la pérdida de carga
        # deltap_mcl = darcy_factor * (self.parameters['Length'] / D) * \
        #     (v**2 / (2 * g))
        # deltap = deltap_mcl * rho * g
        # self.pout = self.pin - deltap
        self.pout = self.pin

    def set_pr_opt(self, aoi):

        IAM = self.sca.get_IAM(aoi)
        pr_opt_peak = self.get_pr_opt_peak()
        self.pr_opt = pr_opt_peak * IAM * np.cos(np.radians(aoi))

    def set_qabs(self, aoi, solarpos, row):
        ''' Total solar power that reach de absorber tube per longitude unit'''

        dni = row[1]['DNI']
        cg = (self.sca.parameters['Aperture'] /
              (np.pi*self.parameters['Absorber tube outer diameter']))

        pr_shadows = self.get_pr_shadows(solarpos)
        pr_borders = self.get_pr_borders(aoi)
        #  Ec. 3.20 Barbero
        self.qabs = self.pr_opt * cg * dni * pr_borders * pr_shadows

    def get_krec(self, t):

        # From Choom S. Kim for A
        #  kt = 100*(0.1241+0.00003279*t)

        # Ec. 4.22 Conductividad para el acero 321H, ver otras opciones.
        kt = 0.0153 * (t - 273.15) + 14.77

        return kt

    def get_urec(self, t, p, htf):

        #  HCE wall thermal conductivity
        krec = self.get_krec(t)

        #  Specific Capacity
        cp = htf.get_cp(t, p)

        #  Internal transmission coefficient.
        hint = self.get_hint(t, p, htf)

        #  Internal heat transfer coefficient. Eq. 3.22 Barbero2016
        return (1 / ((1 / hint) + (
            self.parameters['Absorber tube outer diameter'] *
            np.log(self.parameters['Absorber tube outer diameter'] /
                   self.parameters['Absorber tube inner diameter'])) /
            (2 * krec)))

    def get_previous_pr(self):

        if self.hce_order > 0:
            previous_pr = self.sca.hces[self.hce_order-1].pr
        elif self.sca.sca_order > 0:
            previous_pr = self.sca.loop.scas[self.sca.sca_order-1].hces[-1].pr
        else:
            previous_pr = 1.0

        return previous_pr

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


    def get_pr_borders(self, aoi):

        if aoi > 90:
            pr_borders = 0.0

        else:
            # Llamado "bordes" en Tesis. Pérdidas de los HCE de cabecera según aoi
            sca_unused_length = (self.sca.parameters["Focal Len"] *
                                 np.tan(np.radians(aoi)))

            unused_hces = sca_unused_length // self.parameters["Length"]

            unused_part_of_hce = ((sca_unused_length %
                                   self.parameters["Length"]) /
                                  self.parameters["Length"])

            if self.hce_order < unused_hces:
                pr_borders = 0.0

            elif self.hce_order == unused_hces:
                pr_borders = 1 - ((sca_unused_length % self.parameters["Length"]) /
                                  self.parameters["Length"])
            else:
                pr_borders = 1.0

            if pr_borders > 1.0 or pr_borders < 0.0:
                print("ERROR pr_bordes out of limits", pr_borders)

        return pr_borders


    def get_pr_shadows(self, solarpos):

        if solarpos['elevation'][0] < 0:
            shading = 1

        else:
            shading = 1 - (np.sin(np.radians(abs(solarpos['elevation'][0]))) *
                     self.sca.loop.parameters['row_spacing'] /
                     self.sca.parameters['Aperture'])

        if shading < 0.0 or shading > 1.0:
            shading = 0.0


        shading = 1 - shading


        if  shading > 1 or  shading < 0:
            print("ERROR shading",  shading)

        return shading

    def get_pr_shadows2(self, solarpos):

        sigmabeta = 0.0
        beta0 = 0.0

        if self.sca.loop.parameters['Tracking Type'] == 1:  # N-S single axis tracker
            if solarpos['azimuth'][0] > 0 and solarpos['azimuth'][0] <= 180:
                surface_azimuth = 90  # Surface facing east
            else:
                surface_azimuth = 270  # Surface facing west
        elif self.sca.loop.parameters['Tracking Type'] == 2:  # E-W single axis tracker
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

        #  Surface tilt
        beta = beta0 + 180 * sigmabeta
        # aoi = pvlib.irradiance.aoi(beta,
        #                            surface_azimuth,
        #                            solarpos['zenith'][0],
        #                            solarpos['azimuth'][0])

        Hs = abs(self.sca.parameters['Aperture'] -
               self.sca.loop.parameters['row_spacing'] *
               np.cos(np.radians(90 - beta)))

        Ls = abs(len(self.sca.loop.scas) * self.sca.parameters['SCA Length'] -
                 abs(self.sca.loop.parameters['row_spacing'] *
                     np.tan(np.radians(surface_azimuth -
                                       solarpos['azimuth'][0]))))

        shading = (Ls * Hs / (len(self.sca.loop.scas) *
                                self.sca.parameters['SCA Length'] *
                                self.sca.parameters['Aperture']))

        return shading





    def get_hext(self, wspd):

        #  TO-DO:

        return 0.0

    def get_hint(self, t, p, fluid):

        #  Prandtl number
        # prf = fluid.get_prandtl(t, p)
        # kf = fluid.get_thermal_conductivity(t, p)
        # mu = fluid.get_dynamic_viscosity(t, p)
        # dri = self.parameters['Absorber tube inner diameter']
        # #  Reynolds number for absorber tube inner diameter, dri
        # redri = 4 * self.sca.loop.massflow / (mu * np.pi * dri)
        # #  Nusselt num. Dittus-Boelter correlation. Eq. 4.14 Barbero2016
        # nudb = 0.023 * redri**0.8 * prf** 0.4
        # #  Internal transmission coefficient.
        # hint = kf * nudb / dri

        #  Gnielinski correlation. Eq. 4.15 Barbero2016
        kf = fluid.get_thermal_conductivity(t, p)
        dri = self.parameters['Absorber tube inner diameter']

        #  Prandtl number
        prf = fluid.get_prandtl(t, p)

        #  Reynolds number for absorber tube inner diameter, dri
        redri = fluid.get_Reynolds(dri, t, p, self.sca.loop.massflow)

        # We supose inner wall temperature is equal to fluid temperature
        prri = prf

        #  Skin friction coefficient
        cf = np.power(1.58 * np.log(redri) - 3.28, -2)
        nug = ((0.5 * cf * prf * (redri - 1000)) /
                (1 + 12.7 * np.sqrt(0.5 * cf) * (np.power(prf, 2/3) - 1))) * \
                    np.power(prf / prri, 0.11)

        #  Internal transmission coefficient.
        hint = kf * nug / dri

        return hint


    def get_eext(self, tro, wspd):


        #  Eq. 5.2 Barbero. Parameters given in Pg. 245
        eext = (self.parameters['Absorber emittance factor A0'] +
                self.parameters['Absorber emittance factor A1'] *
                (tro - 273.15))
        """
        Lineal Increase if wind speed lower than 4 m/s up to 1% at 4 m/s
        Lineal increase over 4 m/s up to 2% at 7 m/s
        """
        if wspd < 4:
            eext = eext * (1 + 0.01 * wspd / 4)

        else:
            eext = eext * (1 + 0.01 * (0.3333 * wspd - 0.3333))

        return eext


    def get_absorptance(self):

        return  self.parameters['Absorber absorptance']


    def get_transmittance(self):

        return self.parameters['Envelope transmittance']


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

            #  Cleanliness two times because it affects mirror and envelope
            solarfraction = (self.parameters['Geom.Accuracy'] *
                             self.parameters['Track Twist'] *
                             self.parameters['Cleanliness'] *
                             self.parameters['Cleanliness'] *
                             self.parameters['Factor'] *
                             self.parameters['Availability'])
        else:
            solarfraction = 1.0

        if  solarfraction > 1 or  solarfraction < 0:
            print("ERROR",  solarfraction)

        return solarfraction


    def get_IAM(self, aoi):

        F0 = self.parameters['IAM Coefficient F0']
        F1 = self.parameters['IAM Coefficient F1']
        F2 = self.parameters['IAM Coefficient F2']

        if (aoi > 0 and aoi < 80):
            kiam = (F0 + (F1 * np.radians(aoi) + F2 * np.radians(aoi)**2) /
                    np.cos(np.radians(aoi)))

            if kiam > 1.0:
                kiam = 1.0
        else:
            kiam = 0.0

        if  kiam > 1.0 or  kiam < 0.0:
            print("ERROR",  kiam, aoi)


        # if (aoi > 0 and aoi < 80):
        #     aoi = np.radians(aoi)
        #     kiam = (1- 2.23073E-4 * aoi - 1.1E-4 * aoi**2 +
        #             3.18596E-6 * aoi**3 - 4.85509E-8 * aoi**4)
        # else:
        #     kiam = 0

        # print('aoi', aoi, 'kiam', kiam)
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
        self.qabs = 0.0
        self.qlost = 0.0
        self.qlost_brackets = 0.0
        self.wasted_power = 0.0
        self.pr = 0.0
        self.pr_opt = 0.0

        self.act_tin = 0.0
        self.act_tout = 0.0
        self.act_pin = 0.0
        self.act_pout = 0.0
        self.act_massflow = 0.0  # Actual massflow measured by the flowmeter

        self.wasted_power = 0.0
        self.tracking = True

    def initialize(self, type_of_source, source=None):

        if type_of_source == 'rated':
            self.massflow = self.parameters['rated_massflow']
            self.tin = self.parameters['rated_tin']
            self.pin = self.parameters['rated_pin']

        elif type_of_source == 'subfield' and source is not None:
            self.massflow = source.massflow / len(source.loops)
            self.tin = source.tin
            self.pin = source.pin

        elif type_of_source  == 'solarfield' and source is not None:
            self.massflow = source.massflow / source.total_loops
            self.tin = solarfield.tin
            self.pin = solarfield.pin

        elif type_of_source == 'values' and source is not None:
            self.massflow = source['massflow']
            self.tin = source['tin']
            self.pin = source['pin']

        else:
            print("ERROR initialize()")


    def load_actual(self, subfield = None):

        if subfield == None:
            subfield = self.subfield

        self.act_massflow = subfield.act_massflow / len(subfield.loops)
        self.act_tin = subfield.act_tin
        self.act_tout = subfield.act_tout
        self.act_pin = subfield.act_pin
        self.act_pout = subfield.act_pout

    def set_loop_values_from_HCEs(self):

        pr_list = []
        qabs_list = []
        qlost_brackets_list = []
        qlost_list = []
        pr_opt_list = []

        for s in self.scas:
            for h in s.hces:
                pr_list.append(h.pr)
                qabs_list.append(
                    h.qabs *
                    np.pi *
                    h.parameters['Length'] *
                    h.parameters['Bellows ratio'] *
                    h.parameters['Shield shading'] *
                    h.parameters['Absorber tube outer diameter'])
                qlost_brackets_list.append(
                    h.qlost_brackets *
                    np.pi *
                    h.parameters['Length'] *
                    h.parameters['Absorber tube outer diameter'])
                qlost_list.append(
                    h.qlost *
                    np.pi *
                    h.parameters['Length'] *
                    h.parameters['Absorber tube outer diameter'])
                pr_opt_list.append(h.pr_opt)
        self.pr = np.mean(pr_list)
        self.qabs = np.sum(qabs_list)
        self.qlost_brackets = np.sum(qlost_brackets_list)
        self.qlost = np.sum(qlost_list)
        self.pr_opt = np.mean(pr_opt_list)

    def load_from_base_loop(self, base_loop):

        self.massflow = base_loop.massflow
        self.tin = base_loop.tin
        self.pin = base_loop.pin
        self.tout = base_loop.tout
        self.pout = base_loop.pout
        self.pr = base_loop.pr
        self.wasted_power = base_loop.wasted_power
        self.pr_opt = base_loop.pr_opt
        self.qabs = base_loop.qabs
        self.qlost = base_loop.qlost
        self.qlost_brackets = base_loop.qlost_brackets

    def calc_loop_pr_for_massflow(self, row, solarpos, htf, model):

        for s in self.scas:
            aoi = s.get_aoi(solarpos)
            for h in s.hces:
                h.set_pr_opt(aoi)
                h.set_qabs(aoi, solarpos, row)
                h.set_tin()
                h.set_pin()
                model.calc_pr(h, htf, row)

        self.tout = self.scas[-1].hces[-1].tout
        self.pout = self.scas[-1].hces[-1].pout
        self.set_loop_values_from_HCEs()
        self.set_wasted_power(htf)

    def calc_loop_pr_for_tout(self, row, solarpos, htf, model):

        dri = self.scas[0].hces[0].parameters['Absorber tube inner diameter']
        min_reynolds = self.scas[0].hces[0].parameters['Min Reynolds']

        min_massflow = htf.get_massflow_from_Reynolds(dri, self.tin, self.pin,
                                                      min_reynolds)

        max_error = model.max_err_t  # % desviation tolerance
        search = True
        step = 0

        while search:

            self.calc_loop_pr_for_massflow(row, solarpos, htf, model)

            err = abs(self.tout-self.parameters['rated_tout'])

            if err > max_error and step <1000:
                step += 1
                if self.tout >= self.parameters['rated_tout']:
                    self.massflow *= (1 + err / self.parameters['rated_tout'])
                    search = True
                elif (self.massflow > min_massflow and
                      self.massflow >
                      self.parameters['min_massflow']):
                    self.massflow *= (1 - err / self.parameters['rated_tout'])
                    search = True
                else:
                    self.massflow = max(
                        min_massflow,
                        self.parameters['min_massflow'])
                    self.calc_loop_pr_for_massflow(row, solarpos, htf, model)
                    search = False
            else:
                search = False
                if step>=1000:
                    print("Massflow convergence failed")

        self.tout = self.scas[-1].hces[-1].tout
        self.pout = self.scas[-1].hces[-1].pout
        self.set_loop_values_from_HCEs()
        self.wasted_power = 0.0


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
                             'qlost': h.qlost})

        loop_df = pd.DataFrame(data)

        if parameters:

            loop_df[parameters].plot(
                figsize=(20,10), linewidth=5, fontsize=20)
            plt.xlabel('HCE order', fontsize=20)

        else:

            loop_df[['num', 'tin', 'tout', 'pin', 'pout',
                     'pr', 'qabs', 'qlost']].plot(
                figsize=(20,10), linewidth=5, fontsize=20)

            plt.xlabel('HCE order', fontsize=20)


        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        # pd.set_option('display.max_colwidth', -1)
        print(loop_df[['num', 'tin', 'tout', 'pin', 'pout',
                       'pr', 'qabs', 'qlost']])
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

        # if self.tout > self.parameters['tmax']:
        #     HH = htf.get_deltaH(self.tout, self.pout)
        #     HL = htf.get_deltaH(self.parameters['tmax'], self.pout)
        #     self.wasted_power = self.massflow * (HH - HL)
        # else:
        #     self.wasted_power = 0.0

        if self.tout > self.parameters['tmax']:
            self.wasted_power = htf.get_deltaH2(
                self.tin,
                self.tout,
                self.massflow,
                self.pin)
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

        return 'LP.{0}.{1:03d}'.format(self.subfield.name, self.loop_order)

class Loop(__Loop__):

    def __init__(self, subfield, loop_order, settings):

        self.subfield = subfield
        self.loop_order = loop_order

        super().__init__(settings)


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
            id = 'LB.'+subfield.name
        else:
            id = 'LB.000'

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


    def get_pr_borders(self, aoi):

        if aoi > 90:
            pr_borders = 0.0

        else:
            # Llamado "bordes" en Tesis. Pérdidas de los HCE de cabecera según aoi
            sca_unused_length = (self.parameters_sca["Focal Len"] *
                                 np.tan(np.radians(aoi)))

            pr_borders = 1 - sca_unused_length / (self.parameters['hces'] * \
                                              self.parameters_hce["Length"])

            if pr_borders > 1.0 or pr_borders < 0.0:
                print("ERROR pr_geo out of limits", pr_borders)

        return pr_borders


    def get_pr_shadows(self, solarpos):

        if solarpos['elevation'][0] < 0:
            shading = 1

        else:
            shading = 1 - (np.sin(np.radians(abs(solarpos['elevation'][0]))) *
                     self.parameters['row_spacing'] /
                     self.parameters_sca['Aperture'])

        if shading < 0.0 or shading > 1.0:
            shading = 0.0


        shading = 1 - shading


        if  shading > 1 or  shading < 0:
            print("ERROR shading",  shading)

        return shading

    def get_pr_shadows2(self, solarpos):

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

        #  Surface tilt
        beta = beta0 + 180 * sigmabeta
        # aoi = pvlib.irradiance.aoi(beta,
        #                            surface_azimuth,
        #                            solarpos['zenith'][0],
        #                            solarpos['azimuth'][0])

        Hs = abs(self.parameters_sca['Aperture'] -
               self.parameters['row_spacing'] *
               np.cos(np.radians(beta)))

        Ls = abs(len(self.scas) * self.parameters_sca['SCA Length'] -
                 abs(self.parameters['row_spacing'] *
                     np.tan(np.radians(surface_azimuth -
                                       solarpos['azimuth'][0]))))

        shading = (Ls * Hs / (len(self.scas) *
                                self.parameters_sca['SCA Length'] *
                                self.parameters_sca['Aperture']))

        return shading

    def get_solar_fraction(self):

        #  Cleanliness two times because it affects mirror and envelope
        solarfraction = (self.parameters_sca['Geom.Accuracy'] *
                         self.parameters_sca['Track Twist'] *
                         self.parameters_sca['Cleanliness'] *
                         self.parameters_sca['Cleanliness'] *
                         self.parameters_sca['Factor'] *
                         self.parameters_sca['Availability'])

        if solarfraction > 1 or solarfraction < 0:
            print("ERROR",  solarfraction)

        return solarfraction

    def get_IAM(self, aoi):

        F0 = self.parameters_sca['IAM Coefficient F0']
        F1 = self.parameters_sca['IAM Coefficient F1']
        F2 = self.parameters_sca['IAM Coefficient F2']

        if (aoi > 0 and aoi < 80):
            kiam = (F0 + (F1 * np.radians(aoi) + F2 * np.radians(aoi)**2) /
                    np.cos(np.radians(aoi)))

            if kiam > 1.0:
                kiam = 1.0
        else:
            kiam = 0.0

        if  kiam > 1.0 or  kiam < 0.0:
            print("ERROR",  kiam, aoi)


        # if (aoi > 0 and aoi < 80):
        #     aoi = np.radians(aoi)
        #     kiam = (1- 2.23073E-4 * aoi - 1.1E-4 * aoi**2 +
        #             3.18596E-6 * aoi**3 - 4.85509E-8 * aoi**4)
        # else:
        #     kiam = 0

        # print('aoi', aoi, 'kiam', kiam)
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

        #  Surface tilt
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
        self.qabs = 0.0
        self.qlost = 0.0
        self.qlost_brackets = 0.0
        self.wasted_power = 0.0
        self.pr = 0.0
        self.pr_opt = 0.0

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

    def set_subfield_values_from_loops(self, htf):

        massflow_list = []
        pr_list = []
        pr_opt_list = []
        wasted_power_list = []
        pout_list = []
        enthalpy_list = []
        qlost_list = []
        qabs_list = []
        qlost_brackets_list = []

        for l in self.loops:
            massflow_list.append(l.massflow)
            pr_list.append(l.pr * l.massflow)
            pr_opt_list.append(l.pr_opt * l.massflow)
            wasted_power_list.append(l.wasted_power)
            pout_list.append(l.pout * l.massflow)
            enthalpy_list.append(l.massflow * htf.get_deltaH(l.tout, l.pout))
            qlost_list.append(l.qlost)
            qabs_list.append(l.qabs)
            qlost_brackets_list.append(l.qlost_brackets)

        self.massflow = np.sum(massflow_list)
        self.pr = np.sum(pr_list) / self.massflow
        self.pr_opt = np.sum(pr_opt_list) / self.massflow
        self.wasted_power = np.sum(wasted_power_list) / 1000000
        self.pout = np.sum(pout_list) / self.massflow
        self.tout = htf.get_T(np.sum(enthalpy_list) /
                              self.massflow, self.pout)
        self.qlost = np.sum(qlost_list)
        self.qabs = np.sum(qabs_list)
        self.qlost_brackets = np.sum(qlost_brackets_list)

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

        self.wasted_power = wasted_power / 1000000

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

        self.act_massflow = row[1][self.get_id() +'.a.mf']
        self.act_tin = row[1][self.get_id() +'.a.tin']
        self.act_pin = row[1][self.get_id() +'.a.pin']
        self.act_tout = row[1][self.get_id() +'.a.tout']
        self.act_pout = row[1][self.get_id() +'.a.pout']


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
        self.qabs = 0.0
        self.qlost = 0.0
        self.qlost_brackets = 0.0
        self.wasted_power = 0.0
        self.pr = 0.0
        self.pr_opt = 0.0
        self.pwr = 0.0

        self.act_tin = 0.0
        self.act_tout = 0.0
        self.act_pin = 0.0
        self.act_pout = 0.0
        self.act_massflow = 0.0
        self.act_pwr = 0.0

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

    def set_solarfield_values_from_subfields(self, htf):

        massflow_list = []
        pr_list = []
        pr_opt_list = []
        wasted_power_list = []
        pout_list = []
        enthalpy_list = []
        qlost_list = []
        qabs_list = []
        qlost_brackets_list = []

        for s in self.subfields:
            massflow_list.append(s.massflow)
            pr_list.append(s.pr * s.massflow)
            pr_opt_list.append(s.pr_opt * s.massflow)
            wasted_power_list.append(s.wasted_power)
            pout_list.append(s.pout * s.massflow)
            enthalpy_list.append(s.massflow * htf.get_deltaH(s.tout, s.pout))
            qlost_list.append(s.qlost)
            qabs_list.append(s.qabs)
            qlost_brackets_list.append(s.qlost_brackets)

        self.massflow = np.sum(massflow_list)
        self.pr = np.sum(pr_list) / self.massflow
        self.pr_opt = np.sum(pr_opt_list) / self.massflow
        self.wasted_power = np.sum(wasted_power_list)
        self.pout = np.sum(pout_list) /self.massflow
        self.tout = htf.get_T(np.sum(enthalpy_list) /
                              self.massflow, self.pout)
        self.qlost = np.sum(qlost_list)
        self.qabs = np.sum(qabs_list)
        self.qlost_brackets = np.sum(qlost_brackets_list)

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


    def set_thermal_power(self, htf, datatype):

        self.pwr = self.massflow * (
            htf.get_deltaH(self.tout, self.pout) -
            htf.get_deltaH(self.tin, self.pin))

        self.pwr /= 1000000

        if datatype == 2:
            self.act_pwr = self.act_massflow * (
                htf.get_deltaH(self.act_tout, self.act_pout) -
                htf.get_deltaH(self.act_tin, self.act_pin))

            self.act_pwr /= 1000000

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
        self.first_date = pd.to_datetime(settings['simulation']['first_date'])
        self.last_date = pd.to_datetime(settings['simulation']['last_date'])
        self.report_df = pd.DataFrame()

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

            if self.datatype == 1:  # Because tmy format include TZ info
                naive_datetime = datetime.strptime(
                row[0].strftime('%Y/%m/%d %H:%M'), "%Y/%m/%d %H:%M")
            else:
                naive_datetime = row[0]

            if (naive_datetime < self.first_date or
                naive_datetime > self.last_date):
                pass

            else:

                solarpos = self.site.get_solarposition(row)

                self.gather_general_data(row, solarpos)


                if solarpos['zenith'][0] < 90:
                    self.tracking = True
                else:
                    self.tracking = False

                if self.simulation:
                    self.simulate_solarfield(solarpos, row)
                    self.solarfield.set_thermal_power(self.htf, self.datatype)
                    self.gather_simulation_data(row)

                    str_data = ("SIMULATION: {0} " +
                                "DNI: {1:3.0f} W/m2 Qm: {2:4.1f}kg/s " +
                                "Tin: {3:4.1f}° Tout: {4:4.1f}°C")

                    print(str_data.format(row[0],
                                          row[1]['DNI'],
                                          self.solarfield.massflow,
                                          self.solarfield.tin - 273.15,
                                          self.solarfield.tout - 273.15))


                if self.benchmark and self.datatype == 2:  # 2: Field Data File available
                    self.benchmark_solarfield(solarpos, row)
                    self.solarfield.set_thermal_power(self.htf, self.datatype)
                    self.gather_benchmark_data(row)

                    str_data = ("BENCHMARK: {0} " +
                                "DNI: {1:3.0f} W/m2 act_Qm: {2:4.1f}kg/s " +
                                "act_Tin: {3:4.1f}° act_Tout: {4:4.1f}° " +
                                "Tout: {5:4.1f}°")

                    print(str_data.format(row[0],
                                          row[1]['DNI'],
                                          self.solarfield.act_massflow,
                                          self.solarfield.act_tin - 273.15,
                                          self.solarfield.act_tout - 273.15,
                                          self.solarfield.tout - 273.15))

        self.save_results()

    def simulate_solarfield(self, solarpos, row):

        if self.datatype == 1:
            for s in self.solarfield.subfields:
                s.initialize('rated')
                for l in s.loops:
                    l.initialize('rated')
            self.solarfield.initialize('rated')
            self.base_loop.initialize('rated')
            if self.fastmode:
                if solarpos['zenith'][0] > 90:
                    self.base_loop.massflow = \
                        self.base_loop.parameters['min_massflow']
                    self.base_loop.calc_loop_pr_for_massflow(
                        row, solarpos, self.htf, self.model)
                else:
                    self.base_loop.calc_loop_pr_for_tout(
                        row, solarpos, self.htf, self.model)
                for s in self.solarfield.subfields:
                    for l in s.loops:
                        l.load_from_base_loop(self.base_loop)
            else:
                for s in self.solarfield.subfields:
                    for l in s.loops:
                        if solarpos['zenith'][0] > 90:
                            l.massflow = \
                                self.base_loop.parameters['min_massflow']
                            l.calc_loop_pr_for_massflow(
                                row, solarpos, self.htf, self.model)
                        else:
                            if l.loop_order > 1:
                                # For a faster convergence
                                l.massflow = \
                                l.subfield.loops[l.loop_order - 1].massflow
                            l.calc_loop_pr_for_tout(
                                row, solarpos, self.htf, self.model)

        elif self.datatype == 2:
            #  1st, we initialize subfields because actual data are given for
            #  subfields. 2nd, we initialize solarfield.
            for s in self.solarfield.subfields:
                s.load_actual(row)
                s.initialize('actual')
            self.solarfield.load_actual(self.htf)
            self.solarfield.initialize('actual')
            if self.fastmode:
                # Force minimum massflow at night
                for s in self.solarfield.subfields:
                    self.base_loop.initialize('subfield', s)
                    if solarpos['zenith'][0] > 90:
                        self.base_loop.massflow = \
                            self.base_loop.parameters['min_massflow']
                        self.base_loop.calc_loop_pr_for_massflow(
                            row, solarpos, self.htf, self.model)
                    else:
                        self.base_loop.calc_loop_pr_for_tout(
                            row, solarpos, self.htf, self.model)
                    for l in s.loops:
                        l.load_from_base_loop(self.base_loop)
            else:
                for s in self.solarfield.subfields:
                    for l in s.loops:
                        # l.load_actual(s)
                        l.initialize('subfield', s)
                        if solarpos['zenith'][0] > 90:
                            l.massflow = l.parameters['min_massflow']
                            l.calc_loop_pr_for_massflow(
                                row, solarpos, self.htf, self.model)
                        else:
                            if l.loop_order > 1:
                                # For a faster convergence
                                l.massflow = \
                                l.subfield.loops[l.loop_order - 1].massflow
                            l.calc_loop_pr_for_tout(
                                row, solarpos, self.htf, self.model)

        for s in self.solarfield.subfields:
            s.set_subfield_values_from_loops(self.htf)

        self.solarfield.set_solarfield_values_from_subfields(self.htf)

    def benchmark_solarfield(self, solarpos, row):

        for s in self.solarfield.subfields:
            s.load_actual(row)
            s.initialize('actual')
        self.solarfield.load_actual(self.htf)
        self.solarfield.initialize('actual')

        if self.fastmode:
            for s in self.solarfield.subfields:
                self.base_loop.initialize('subfield', s)
                self.base_loop.calc_loop_pr_for_massflow(
                    row, solarpos, self.htf, self.model)
                self.base_loop.set_loop_values_from_HCEs()

                for l in s.loops:
                    l.load_from_base_loop(self.base_loop)

                s.set_subfield_values_from_loops(self.htf)

        else:
            for s in self.solarfield.subfields:
                for l in s.loops:
                    # l.load_actual()
                    l.initialize('subfield', s)
                    l.calc_loop_pr_for_massflow(
                        row, solarpos, self.htf, self.model)
                    l.set_loop_values_from_HCEs()

                s.set_subfield_values_from_loops(self.htf)

        self.solarfield.set_solarfield_values_from_subfields(self.htf)


    def store_values(self, row, values):

        for v in values:
            self.datasource.dataframe.at[row[0], v] = values[v]

    def gather_general_data(self, row, solarpos):

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
        self.datasource.dataframe.at[row[0], 'pr_shadows'] = \
            self.base_loop.get_pr_shadows(solarpos)
        self.datasource.dataframe.at[row[0], 'pr_borders'] = \
            self.base_loop.get_pr_borders(aoi)
        self.datasource.dataframe.at[row[0], 'pr_opt_peak'] = \
            self.base_loop.get_pr_opt_peak()
        self.datasource.dataframe.at[row[0], 'solar_fraction'] = \
            self.base_loop.get_solar_fraction()


    def gather_simulation_data(self, row):

        #  Solarfield data
        self.datasource.dataframe.at[row[0], 'SF.x.mf'] = \
            self.solarfield.massflow
        self.datasource.dataframe.at[row[0], 'SF.x.tin'] = \
            self.solarfield.tin - 273.15
        self.datasource.dataframe.at[row[0], 'SF.x.tout'] = \
            self.solarfield.tout - 273.15
        self.datasource.dataframe.at[row[0], 'SF.x.pin'] = \
            self.solarfield.pin
        self.datasource.dataframe.at[row[0], 'SF.x.pout'] = \
            self.solarfield.pout
        self.datasource.dataframe.at[row[0], 'SF.x.prth'] = \
            self.solarfield.pr
        self.datasource.dataframe.at[row[0], 'SF.x.prop'] = \
            self.solarfield.pr_opt
        self.datasource.dataframe.at[row[0], 'SF.x.qabs'] = \
            self.solarfield.qabs
        self.datasource.dataframe.at[row[0], 'SF.x.qlst'] = \
            self.solarfield.qlost
        self.datasource.dataframe.at[row[0], 'SF.x.qlbk'] = \
            self.solarfield.qlost_brackets
        self.datasource.dataframe.at[row[0], 'SF.x.pwr'] = \
            self.solarfield.pwr

        if  self.datatype == 2:
            if row[1]['GrossPower']>0:
                self.datasource.dataframe.at[row[0], 'SF.x.globalpr'] = \
                    row[1]['GrossPower'] / self.solarfield.pwr
            else:
                self.datasource.dataframe.at[row[0], 'SF.x.globalpr'] = 0
        else:
            self.datasource.dataframe.at[row[0], 'SF.x.globalpr'] = 0

        if self.fastmode:

            values = {
                self.base_loop.get_id() + '.x.mf': self.base_loop.massflow,
                self.base_loop.get_id() + '.x.tin':
                    self.base_loop.tin - 273.15,
                self.base_loop.get_id() + '.x.tout':
                    self.base_loop.tout  - 273.15,
                self.base_loop.get_id() + '.x.pin': self.base_loop.pin,
                self.base_loop.get_id() + '.x.pout': self.base_loop.pout,
                self.base_loop.get_id() + '.x.prth': self.base_loop.pr,
                self.base_loop.get_id() + '.x.prop': self.base_loop.pr_opt,
                self.base_loop.get_id() + '.x.qabs': self.base_loop.qabs,
                self.base_loop.get_id() + '.x.qlst': self.base_loop.qlost,
                self.base_loop.get_id() + '.x.qlbk': \
                    self.base_loop.qlost_brackets}
            self.store_values(row, values)

        for s in self.solarfield.subfields:
            #  Agretate data from subfields
            values = {
                s.get_id() + '.x.mf': s.massflow,
                s.get_id() + '.x.tin': s.tin  - 273.15,
                s.get_id() + '.x.tout': s.tout - 273.15,
                s.get_id() + '.x.pin': s.pin,
                s.get_id() + '.x.pout': s.pout,
                s.get_id() + '.x.prth': s.pr,
                s.get_id() + '.x.prop': s.pr_opt,
                s.get_id() + '.x.qabs': s.qabs,
                s.get_id() + '.x.qlst': s.qlost,
                s.get_id() + '.x.qlbk': s.qlost_brackets}

            self.store_values(row, values)

            if not self.fastmode:
                for l in s.loops:
                    #  Loop data
                    values = {
                        l.get_id() + '.x.mf':  l.massflow,
                        l.get_id() + '.x.tin':  l.tin - 273.15,
                        l.get_id() + '.x.tout':  l.tout - 273.15,
                        l.get_id() + '.x.pin':  l.pin,
                        l.get_id() + '.x.pout':  l.pout,
                        l.get_id() + '.x.prth':  l.pr,
                        l.get_id() + '.x.prop':  l.pr_opt,
                        l.get_id() + '.x.qabs':  l.qabs,
                        l.get_id() + '.x.qlost':  l.qlost,
                        l.get_id() + '.x.qlbk':  l.qlost_brackets}

                    self.store_values(row, values)

    def gather_benchmark_data(self, row):

        #  Solarfield data
        self.datasource.dataframe.at[row[0], 'SF.a.mf'] = \
            self.solarfield.massflow
        self.datasource.dataframe.at[row[0], 'SF.a.tin'] = \
            self.solarfield.tin - 273.15
        self.datasource.dataframe.at[row[0], 'SF.a.tout'] = \
            self.solarfield.act_tout - 273.15
        self.datasource.dataframe.at[row[0], 'SF.a.pwr'] = \
            self.solarfield.act_pwr
        self.datasource.dataframe.at[row[0], 'SF.b.tout'] = \
            self.solarfield.tout - 273.15
        self.datasource.dataframe.at[row[0], 'SF.b.prth'] = \
            self.solarfield.pr
        self.datasource.dataframe.at[row[0], 'SF.b.prop'] = \
            self.solarfield.pr_opt
        self.datasource.dataframe.at[row[0], 'SF.b.pwr'] = \
            self.solarfield.pwr
        self.datasource.dataframe.at[row[0], 'SF.b.wpwr'] = \
            self.solarfield.wasted_power
        self.datasource.dataframe.at[row[0], 'SF.a.pin'] = \
            self.solarfield.pin
        self.datasource.dataframe.at[row[0], 'SF.a.pout'] = \
            self.solarfield.act_pout
        self.datasource.dataframe.at[row[0], 'SF.b.pout'] = \
            self.solarfield.pout
        self.datasource.dataframe.at[row[0], 'SF.b.qabs'] = \
            self.solarfield.qabs
        self.datasource.dataframe.at[row[0], 'SF.b.qlst'] = \
            self.solarfield.qlost
        self.datasource.dataframe.at[row[0], 'SF.b.qlbk'] = \
            self.solarfield.qlost_brackets

        if self.solarfield.qabs > 0:
            self.datasource.dataframe.at[row[0], 'SF.a.prth'] = \
                10e6 * self.solarfield.act_pwr / self.solarfield.qabs
        else:
            self.datasource.dataframe.at[row[0], 'SF.a.prth'] = 0

        if row[1]['GrossPower']>0:
            self.datasource.dataframe.at[row[0], 'SF.a.globalpr'] = \
                row[1]['GrossPower'] / self.solarfield.act_pwr
        else:
            self.datasource.dataframe.at[row[0], 'SF.a.globalpr'] = 0

        if row[1]['GrossPower']>0:
            self.datasource.dataframe.at[row[0], 'SF.b.globalpr'] = \
                row[1]['GrossPower'] / self.solarfield.pwr
        else:
            self.datasource.dataframe.at[row[0], 'SF.a.globalpr'] = 0

        for s in self.solarfield.subfields:

            if self.fastmode:
                values = {
                    self.base_loop.get_id(s) + '.a.mf':
                        self.base_loop.massflow,
                    self.base_loop.get_id(s) + '.a.tin':
                        self.base_loop.tin - 273.15,
                    self.base_loop.get_id(s) + '.a.tout':
                        self.base_loop.act_tout - 273.15,
                    self.base_loop.get_id(s) + '.b.tout':
                        self.base_loop.tout - 273.15,
                    self.base_loop.get_id(s) + '.a.pin': self.base_loop.pin,
                    self.base_loop.get_id(s) + '.b.pout': self.base_loop.pout,
                    self.base_loop.get_id(s) + '.b.prth': self.base_loop.pr,
                    self.base_loop.get_id(s) + '.b.prop':
                        self.base_loop.pr_opt,
                    self.base_loop.get_id(s) + '.b.qabs': self.base_loop.qabs,
                    self.base_loop.get_id(s) + '.b.qlst': self.base_loop.qlost,
                    self.base_loop.get_id(s) + '.b.qlbk':
                        self.base_loop.qlost_brackets,
                    self.base_loop.get_id(s) + '.b.wpwr':
                        self.base_loop.wasted_power}
                self.store_values(row, values)

            #  Agretate data from subfields
            values = {
                s.get_id() + '.a.mf': s.act_massflow,
                s.get_id() + '.a.tin': s.act_tin - 273.15,
                s.get_id() + '.a.tout': s.act_tout - 273.15,
                s.get_id() + '.b.tout': s.tout - 273.15,
                s.get_id() + '.a.pin': s.act_pin,
                s.get_id() + '.a.pout': s.act_pout,
                s.get_id() + '.b.pout': s.pout,
                s.get_id() + '.b.prth': s.pr,
                s.get_id() + '.b.prop': s.pr_opt,
                s.get_id() + '.b.qabs': s.qabs,
                s.get_id() + '.b.qlst': s.qlost,
                s.get_id() + '.b.qlbk': s.qlost_brackets,
                s.get_id() + '.b.wpwr': s.wasted_power}

            self.store_values(row, values)

            if not self.fastmode:
                for l in s.loops:
                    #  Loop data
                    values = {
                        l.get_id() + '.a.mf':  l.act_massflow,
                        l.get_id() + '.a.tin':  l.act_tin - 273.15,
                        l.get_id() + '.a.tout':  l.act_tout - 273.15,
                        l.get_id() + '.b.tout':  l.tout - 273.15,
                        l.get_id() + '.a.pin':  l.act_pin,
                        l.get_id() + '.a.pout':  l.act_pout,
                        l.get_id() + '.b.pout':  l.pout,
                        l.get_id() + '.b.prth':  l.pr,
                        l.get_id() + '.b.prop':  l.pr_opt,
                        l.get_id() + '.b.qabs':  l.qabs,
                        l.get_id() + '.b.qlost':  l.qlost,
                        l.get_id() + '.b.qlbk':  l.qlost_brackets,
                        l.get_id() + '.b.wpwr': l.wasted_power}

                    self.store_values(row, values)


    def show_report(self, keys= None):

        # from IPython.display import display, Math, Latex
        # display(Math(r'F(k) = \int_{-\infty}^{\infty} f(x) e^{2\pi i k} dx'))


        self.report_df[keys].plot(
                        figsize=(20,10), linewidth=5, fontsize=20)
        plt.xlabel('Date', fontsize=20)
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)

    def save_results(self):

        keys = ['DNI','elevation','zenith','azimuth','aoi','IAM',
                'pr_shadows','pr_borders','pr_opt_peak','solar_fraction']

        keys_graphics_power = ['DNI']
        keys_graphics_temp = ['DNI']

        keys_a = ['NetPower','AuxPower','GrossPower',
                  'SF.a.mf','SF.a.tin', 'SF.a.tout',
                  'SF.a.pwr','SF.a.prth','SF.a.globalpr']
        keys_graphics_power_a = ['NetPower','AuxPower','GrossPower',
                                 'SF.a.pwr']
        keys_graphics_temp_a = ['SF.a.mf','SF.a.tin', 'SF.a.tout']

        keys_x = ['SF.x.mf','SF.x.tin','SF.x.tout', 'SF.x.pwr',
                  'SF.x.prth','SF.x.globalpr']
        keys_graphics_power_x = ['SF.x.pwr']
        keys_graphics_temp_x = ['SF.x.mf', 'SF.x.tout']

        keys_b = ['SF.b.tout','SF.b.pwr', 'SF.b.wpwr',
                  'SF.b.prth','SF.b.globalpr']
        keys_graphics_power_b = ['SF.b.pwr']
        keys_graphics_temp_b = ['SF.b.tout']

        if self.datatype == 2:
            keys += keys_a
            keys_graphics_power += keys_graphics_power_a
            keys_graphics_temp += keys_graphics_temp_a

        if self.simulation == True:
            keys += keys_x
            keys_graphics_power += keys_graphics_power_x
            keys_graphics_temp += keys_graphics_temp_x

        if self.benchmark == True:
            keys += keys_b
            keys_graphics_power += keys_graphics_power_b
            keys_graphics_temp += keys_graphics_temp_b

        self.report_df = self.datasource.dataframe

        # self.report_df = self.datasource.dataframe[
        #     (self.datasource.dataframe.index >= self.first_date) &
        #     (self.datasource.dataframe.index <= self.last_date)]

        # self.report_df = self.report_df[keys]


        try:
            initialdir = "./simulations_outputs/"
            prefix = datetime.today().strftime("%Y%m%d %H%M%S ")
            filename_complete = str(self.ID) + "_COMPLETE"
            filename_report = str(self.ID) + "_REPORT"
            sufix = ".csv"

            path_complete = initialdir + prefix + filename_complete + sufix
            path_report = initialdir + prefix + filename_report + sufix

            self.datasource.dataframe.to_csv(
                path_complete, sep=';', decimal = ',')
            self.report_df.to_csv(path_report, sep=';', decimal = ',')

        except Exception:
            raise
            print('Error saving results, unable to save file: %r', path)

        # self.show_report(keys_graphics_power)
        # self.show_report(keys_graphics_temp)


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

        print("Running simulation for source data file: {0} from: \
              {1} to {2}".format(
            self.parameters['simulation']['filename'],
            self.parameters['simulation']['first_date'],
            self.parameters['simulation']['last_date']))
        print("Model: {0}".format(
            self.parameters['model']['name']))
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
        print("--------------------------------------------------------------")


class LoopSimulation(object):
    """
    Definimos la clase simulacion para representar las diferentes
    pruebas que lancemos, variando el archivo TMY, la configuracion del
    site, la planta, el modo de operacion o el modelo empleado.

    """

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

    _T_REF = 285.856  # K, T_REF= 12.706 K
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

    def get_T2(self, tin, q, mf=None, p=None):
        pass

    def get_dynamic_viscosity(self, t, p):
        pass

    def get_Reynolds(self, dri, t, p, massflow):

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

        # #  Prandtl number
        # prandtl = mu / alpha
        prandtl = cp * mu / kf

        return prandtl

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

        #p = 1600000
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

        CP.set_reference_state(self.coolpropID,'ASHRAE')
        temperature = PropsSI('T', 'H', h, 'P', p, self.coolpropID)
        CP.set_reference_state(self.coolpropID, 'DEF')

        return temperature

    def get_T2(self, t,  q, mf = None, p = None):

        if t > self.tmax:
            t = self.tmax

        CP.set_reference_state(self.coolpropID,'ASHRAE')
        hin = PropsSI('H', 'T', t,  'P', p, self.coolpropID)
        try:
            temperature = PropsSI('T', 'H', hin  + q/mf, 'P', p, self.coolpropID)
        except:
            temperature = self.tmax
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


    def get_density(self, t, p):

        # Dowtherm A.pdf, 2.2 Single Phase Liquid Properties. pg. 8.

        poly = np.polynomial.polynomial.Polynomial(self.rho)

        # if t > self.tmax:
        #     t= self.tmax

        return poly(t) * (p * 1.0e-4)**1.0e-3


    def get_dynamic_viscosity(self, t, p):

        poly = np.polynomial.polynomial.Polynomial(self.mu)

        if t > self.tmax:
            mu_  = poly(self.tmax)
        else:
            mu_ = poly(t)
        # mu_ = 0.00012
        # mu_ = poly(t)
        return mu_

        # return mu_

    def get_cp(self, t, p):

        # cp0, cp1, cp2, cp3, cp4, cp5 = tuple(self.cp)

        # return (cp0 + cp1 * t + cp2 * t**2 + cp3 * t**3 +
        #         cp4 * t**4 + cp5 * t**5)

        poly = np.polynomial.polynomial.Polynomial(self.cp)
        if t > self.tmax:
            t= self.tmax

        return poly(t)


    def get_thermal_conductivity(self, t, p):
        ''' Saturated Fluid conductivity at temperature t '''

        poly = np.polynomial.polynomial.Polynomial(self.kt)
        # if t > self.tmax:
        #     t= self.tmax

        return poly(t)

    def get_deltaH(self, t, p):

        poly = np.polynomial.polynomial.Polynomial(self.h)

        # if t > self.tmax:
        #     t= self.tmax

        return poly(t) - poly(self._T_REF)

    def get_deltaH2(self, t1, t2, mf=None, p=None):

        cp0, cp1, cp2, cp3, cp4, cp5, cp6, cp7, cp8 = tuple(self.cp)

        h = (mf * (
            (cp0 * t2 + cp1 * t2**2 / 2 + cp2 * t2**3 / 3 +
             cp3 * t2**4 / 4 + cp4 * t2**5 / 5 + cp5 * t2**6 / 6 +
             cp6 * t2**7 / 7 + cp7 * t2**8 / 8 + cp8 * t2**9 / 9)
            -
            (cp0 * t1 + cp1 * t1**2 / 2 + cp2 * t1**3 / 3 +
             cp3 * t1**4 / 4 + cp4 * t1**5 / 5 + cp5 * t1**6 / 6 +
             cp6 * t1**7 / 7 + cp7 * t1**8 / 8 + cp8 * t1**9 / 9)))

        return h

    def get_T(self, h, p):

        poly = np.polynomial.polynomial.Polynomial(self.t)

        return poly(h)

    def get_T2(self, tin, h, mf=None, p=None):


        #tout = tin
        cp0, cp1, cp2, cp3, cp4, cp5, cp6, cp7, cp8 = tuple(self.cp)

        a0 = (h/mf + cp0 * tin + cp1 * tin**2 / 2 + cp2 * tin**3 / 3 +
             cp3 * tin**4 / 4 + cp4 * tin**5 / 5 + cp5 * tin**6 / 6 +
             cp6 * tin**7 / 7 + cp7 * tin**8 / 8 + cp8 * tin**9 / 9)

        factors = [a0, -cp0, -cp1 / 2, -cp2 / 3, -cp3 / 4, -cp4 / 5, -cp5 / 6,
                   -cp6 / 7, -cp7 / 8, -cp8 / 9]

        poly = np.polynomial.polynomial.Polynomial(factors)
        roots = poly.roots()

        tout_bigger = []
        tout_smaller = []

        for r in roots:
            if r.imag == 0.0:
                if r.real >= tin:
                    tout_bigger.append(r.real)
                else:
                    tout_smaller.append(r.real)
        if h > 0:
            tout = min(tout_bigger)
        elif h<0:
            tout = max(tout_smaller)
        else:
            tout = tin

        # for r in roots:
        #     if r.imag == 0 and r.real >= 0:
        #         tout = r.real

        return tout


class Weather(object):

    def __init__(self, settings = None):

        self.dataframe = None
        self.site = None
        self.weatherdata = None

        if settings is not None:
            # self.filename = settings['filename']
            # self.filepath = settings['filepath']
            # self.file = self.filepath + self.filename
            self.openWeatherDataFile(settings['filepath'] +
                                     settings['filename'])
            # self.file
        else:
            self.openWeatherDataFile()

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

        needed_columns = ['DNI', 'DryBulb', 'DewPoint', 'Wspd', 'Pressure']
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
            # if self.datatype == 1:  # Because tmy format include TZ info
            #     naive_datetime = datetime.strptime(
            #     row[0].strftime('%Y/%m/%d %H:%M'), "%Y/%m/%d %H:%M")
            # else:
            #     naive_datetime = row[0]

            # if (naive_datetime < self.first_date or
            #     naive_datetime > self.last_date):
            #     self.datasource.dataframe.drop(row[0], axis=0)



    def __init__(self, settings, tags = None):
        self.filename = settings['filename']
        self.filepath = settings['filepath']
        self.file = self.filepath + self.filename
        self.first_date = pd.to_datetime(settings['first_date'])
        self.last_date = pd.to_datetime(settings['last_date'])
        self.tags = tags
        self.dataframe = None

        self.openFieldDataFile(self.file)
        self.rename_columns()
        self.change_units()


    def openFieldDataFile(self, path = None):

        '''
        fielddata
        '''

        rows_list =[]
        index_count = 1  # Skip fist row in skiprows: it has got columns names
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
                    df = pd.read_csv(
                        path, sep=';',
                        decimal= ',',
                        #dayfirst=True,
                        #index_col=0,
                        usecols=[0])

                    df = pd.to_datetime(
                        df['date'], format = "%d/%m/%Y %H:%M")

                    for row in df:
                        if (row < self.first_date or
                            row > self.last_date):
                            rows_list.append(index_count)
                        index_count += 1

                    self.dataframe = pd.read_csv(
                        path, sep=';',
                        decimal= ',',
                        dayfirst=True,
                        index_col=0,
                        skiprows=rows_list)

                    self.file = path

                else:
                    print("unknow extension ", strext)
                    return

        except Exception:
            raise
            txMessageBox.showerror('Error loading FieldData File',
                                   'Unable to open file: %r', self.file)

        self.dataframe.index = pd.to_datetime(self.dataframe.index,
                                              format= "%d/%m/%Y %H:%M")





        # try:
        #     if path is None:
        #         root = Tk()
        #         root.withdraw()
        #         path = askopenfilename(initialdir = ".fielddata_files/",
        #                            title = "choose your file",
        #                            filetypes = (("csv files","*.csv"),
        #                                         ("all files","*.*")))
        #         root.update()
        #         root.destroy()
        #         if path is None:
        #             return
        #         else:
        #             strfilename, strext = os.path.splitext(path)
        #             if  strext == ".csv":
        #                 self.dataframe = pd.read_csv(
        #                     path, sep=';',
        #                     decimal= ',',
        #                     dayfirst=True,
        #                     index_col=0)

        #                 self.file = path
        #             elif strext == ".xls":

        #                 self.dataframe = pd.read_excel(path)
        #                 self.file = path
        #             else:
        #                 print("unknow extension ", strext)
        #                 return
        #     else:
        #         strfilename, strext = os.path.splitext(path)

        #         if  strext == ".csv":
        #             self.dataframe = pd.read_csv(
        #                 path, sep=';',
        #                 decimal= ',',
        #                 dayfirst=True,
        #                 index_col=0)

        #             self.file = path
        #         elif strext == ".xls":
        #             self.dataframe = pd.read_excel(path)
        #             self.file = path
        #         else:
        #             print("unknow extension ", strext)
        #             return
        # except Exception:
        #     raise
        #     txMessageBox.showerror('Error loading FieldData File',
        #                            'Unable to open file: %r', self.file)

        # self.dataframe.index = pd.to_datetime(self.dataframe.index,
        #                                       format= "%d/%m/%Y %H:%M")

    def change_units(self):

        for c in self.dataframe.columns:
            if ('.a.t' in c) or ('DryBulb' in c) or ('Dew' in c):
                self.dataframe[c] += 273.15 # From Celsius Degrees to K
            if '.a.p' in c:
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
    Power Plant as a set composed by a solarfield, a HeatExchanger, a
    PowerCycle and a BOPSystem

    '''

    def __init__(self, settings):

        pass


class Site(object):
    def __init__(self, settings):

        self.name = settings['name']
        self.latitude = settings['latitude']
        self.longitude = settings['longitude']
        self.altitude = settings['altitude']


    def get_solarposition(self, row):

        solarpos = pvlib.solarposition.get_solarposition(
            row[0] + timedelta(hours=0.5),
            self.latitude,
            self.longitude,
            self.altitude,
            pressure=row[1]['Pressure'],
            temperature=row[1]['DryBulb'])

        return solarpos

    # def get_hour_angle(self, row, equiation_of_time):

    #     hour_angle = pvlib.solarposition.hour_angle(
    #         row[0] + timedelta(hours=0.5),
    #         self.longitude,
    #         equation_of_time)

    #     return hour_angle

    # def get_equation_of_time(self, row):



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


class Test(object):
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
        self.first_date = pd.to_datetime(settings['simulation']['first_date'])
        self.last_date = pd.to_datetime(settings['simulation']['last_date'])
        self.report_df = pd.DataFrame()

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
        elif self.datatype ==3:
            self.datasource = TableData(settings['simulation'])

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


    def run_test(self):

        print('running test')

        print(self.datasource)

        for row in self.datasource.dataframe.iterrows():

            self.base_loop.massflow = row[1]['mf']
            kf = self.htf.get_thermal_conductivity(
                row[1]['tin'], row[1]['pin'])
            dri = self.base_loop.parameters_hce['Absorber tube inner diameter']

            #  Prandtl number
            prf = self.htf.get_prandtl(row[1]['tin'], row[1]['pin'])
            mu = self.htf.get_dynamic_viscosity(row[1]['tin'], row[1]['pin'])

            #  Reynolds number for absorber tube inner diameter, dri
            redri = self.htf.get_Reynolds(dri, row[1]['tin'], row[1]['pin'],
                                          row[1]['mf'])
            cp = self.htf.get_cp(row[1]['tin'], row[1]['pin'])
            rho = self.htf.get_density(row[1]['tin'], row[1]['pin'])
            # We supose inner wall temperature is equal to fluid temperature
            prri = prf
            dro = self.base_loop.parameters_hce['Absorber tube outer diameter']


            #  HCE wall thermal conductivity
            krec = self.base_loop.scas[0].hces[0].get_krec(row[1]['tin'])
            hint = self.base_loop.scas[0].hces[0].get_hint(
                    row[1]['tin'], row[1]['pin'], self.htf)
            urec = 1 / ((1 / hint) + (dro * np.log(dro / dri)) / (2 * krec))
            #  Skin friction coefficient
            cf = np.power(1.58 * np.log(redri) - 3.28, -2)

            #  Gnielinski correlation. Eq. 4.15 Barbero2016
            nug = ((0.5 * cf * prf * (redri - 1000)) /
                   (1 + 12.7 * np.sqrt(0.5 * cf) * (np.power(prf, 2/3) - 1))) * \
                       np.power(prf / prri, 0.11)

            self.datasource.dataframe.at[row[0], 'h_int'] = \
                self.base_loop.scas[0].hces[0].get_hint(
                    row[1]['tin'], row[1]['pin'], self.htf)
            self.datasource.dataframe.at[row[0], 'krec'] = \
                self.base_loop.scas[0].hces[0].get_krec(
                    row[1]['tin'])
            self.datasource.dataframe.at[row[0], 'urec'] = urec
            self.datasource.dataframe.at[row[0], 'krec'] = krec
            self.datasource.dataframe.at[row[0], 'qloss'] = \
                self.base_loop.scas[0].hces[0].get_qlost_brackets(
                    row[1]['tin'], row[1]['DryBulb'])
            self.datasource.dataframe.at[row[0], 'mu'] = mu
            self.datasource.dataframe.at[row[0], 'cp'] = cp
            self.datasource.dataframe.at[row[0], 'rho'] = rho
            self.datasource.dataframe.at[row[0], 'nug'] = nug
            self.datasource.dataframe.at[row[0], 'cf'] = cf
            self.datasource.dataframe.at[row[0], 'redri'] = redri
            self.datasource.dataframe.at[row[0], 'prf'] = prf
            self.datasource.dataframe.at[row[0], 'kf'] = kf
            self.datasource.dataframe.at[row[0], 'pr_opt_peak'] = \
                self.base_loop.get_pr_opt_peak()
            self.datasource.dataframe.at[row[0], 'solar_fraction'] = \
                self.base_loop.get_solar_fraction()

        try:
            initialdir = "./simulations_outputs/"
            prefix = datetime.today().strftime("%Y%m%d %H%M%S ")
            filename_test = str(self.ID) + "_TEST"
            sufix = ".csv"

            path_complete = initialdir + prefix + filename_test + sufix

            self.datasource.dataframe.to_csv(
                path_complete, sep=';', decimal = ',')

        except Exception:
            raise
            print('Error saving results, unable to save file: %r', path)

# if __name__=='__main__':

#     cfg = sys.argv[1:]

#     SIM = cs.SolarFieldSimulation(cfg)

#     FLAG_00 = datetime.now()
#     SIM.runSimulation()
#     FLAG_01 = datetime.now()
#     DELTA_01 = FLAG_01 - FLAG_00
#     print("Total runtime: ", DELTA_01.total_seconds())

