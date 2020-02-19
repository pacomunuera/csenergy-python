import os, sys
import pvlib
from pvlib import pvsystem
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn
import inspect
import numpy as np
import pandas as pd
import io
import requests

lat = 24.133889
lon = 54.520833
altitude = 31
tz = 'Asia/Dubai'
surface_azimuth = 360
albedo = 0.35
#Varibles related to Single Axix Tracker Function
axis_tilt = 0
axis_azimuth = 180
max_angle = 90
gcr = 269/500

file = r"C:\Users\rf5896\Desktop\Daniel Sir\Jawairia\Solargis_TMY60_DATA_ORIGINAL_P50.csv"
df = pd.read_csv(file, skiprows = 56, header = None, sep = '\n')
df = df[0].str.split(';', expand=True)
header = df.iloc[0]

tmy = pd.DataFrame(df.values[1:], columns = header)

concat_Datetime = tmy['Day'].astype(str) + ' ' + tmy['Time'].astype(str)

tmy['timestamp'] = concat_Datetime
tmy['timestamp'] = pd.to_datetime(tmy['timestamp'], format = '%j %H:%M')


tmy.set_index(tmy['timestamp'],inplace=True)

tmy.drop(tmy.columns[[0,1,5,6,9,10,12,13]], axis = 1, inplace = True)
tmy.columns = ['ghi','dni','dhi','temp_air','wind_speed','pressure']
tmy = tmy.apply(pd.to_numeric)
tmy.columns = tmy.columns.str.replace(' ', '')

tmy['Month'] = tmy.index
tmy.Month = tmy.Month.dt.month	
monthly_tmy = tmy.groupby('Month').sum()/1000

loc = pvlib.location.Location(latitude=lat,longitude= lon,tz=tz, 
                                     altitude = altitude, name = 'Abu Dhabi')

solpos=pvlib.solarposition.get_solarposition(tmy.index, loc.latitude, loc.longitude,
                                               loc.altitude)

dni_extraterresial = pvlib.irradiance.get_extra_radiation(tmy.index)

airmass = pvlib.atmosphere.get_relative_airmass(zenith = solpos['apparent_zenith'])

single_axis_tracker = pvlib.tracking.singleaxis(apparent_zenith = solpos['apparent_zenith'], 
                                                apparent_azimuth = solpos['azimuth'], 
                                                axis_tilt = axis_tilt, 
                                                axis_azimuth = axis_azimuth, 
                                                max_angle = max_angle, backtrack = True, 
                                                gcr = gcr)


poa=pvlib.irradiance.get_total_irradiance(surface_tilt = single_axis_tracker['surface_tilt'], 
                                          surface_azimuth = 180,
                                          solar_zenith = solpos['apparent_zenith'],
                                          solar_azimuth = solpos['azimuth'] , dni = tmy['dni'], 
                                          ghi = tmy['ghi'], dhi = tmy['dhi'],
                                          dni_extra = dni_extraterresial, airmass = airmass,
                                          albedo = albedo, 
                                          model='perez',
                                          model_perez = 'allsitescomposite1990') 

poa['Month'] = poa.index
poa.Month = poa.Month.dt.month	
monthly_poa = poa.groupby('Month').sum()/1000



