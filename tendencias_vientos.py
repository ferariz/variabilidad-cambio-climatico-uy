#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 de tarde en Camelia 2019

Calculamos las tendencias lineales en los vientos en la región 15-45S, 90-30W.
Se hace en el marco del trabajo de los SLP WR y la variabilidad climática

@author: fernando
"""

####################################################################################################
##### Importamos las bibliotecas necesarias para la adquisición de datos y manejo de variables #####
####################################################################################################

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

import iris
import iris.coord_categorisation
import iris.quickplot as qplt

import cartopy
import cartopy.feature as cfeat
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


####################################################################################################
#### Adquirimos los datos de presión superficial diarios del reanalisis de NCEP CDAS1 ##############
#### La grilla tiene una resolución de  2.5∘ , que resulta en  N=144×73=10512 nodos ################
#### Dado que los extremos de latitud son exactamente en los polos, N=10226 ########################
#### El período temporal va desde el 1 ene 1948 al 16 dic 2018,  T′=25971  #########################
#### Vale aclarar que no vamos a tener en cuenta los 29 de febrero, entonces  T=27953 ##############
#### (años bisiestos: 1948, 1952, 1956, 1960, 1964, 1968, 1972, 1976, 1980, 1984, 1988, 1992, 1996,#
#### 2000, 2004, 2008, 2012, 2016).#################################################################

filename = '/home/meteo/datos/Uy/uwnd.1000mb.mon.mean.9030W1545S.nc'
filename2 = '/home/meteo/datos/Uy/vwnd.1000mb.mon.mean.9030W1545S.nc'

ncin = Dataset(filename, 'r')
u = ncin.variables['u'][:]
lons = ncin.variables['X'][:]
lats = ncin.variables['Y'][:]
ncin.close()

ncin = Dataset(filename2, 'r')
v = ncin.variables['v'][:]
ncin.close()

u = np.squeeze(u)
v = np.squeeze(v)


####################################################################################################
#####                          Selección del dominio espacial                                  #####
#####      En nuestro caso nos vamos a concentrar en la región  300E- 310E y  27S- 37S         #####
####################################################################################################

lat = lats[5:10]
lon = lons[11:17]
nlat = 5
nlon = 6
u = u[:,5:10,11:17]
v = v[:,5:10,11:17]

t_e = np.arange(8,u.shape[0],12)
mes = 'septiembre'
u_e = np.squeeze(u[t_e.astype(int),:,:])
v_e = np.squeeze(v[t_e.astype(int),:,:])

####################################################################################################
#####                          Descontamos la tendencia lineal                                 #####
####################################################################################################
t_medio = u_e.shape[0]/2
den = 0
# calculamos valores medios
for t in range(u_e.shape[0]):
    den += (t-t_medio)*(t-t_medio)
u_medio = np.zeros((u_e.shape[1],u_e.shape[2]))
v_medio = np.zeros((v_e.shape[1],v_e.shape[2]))
for i in range(u.shape[1]):
    print(i)
    for j in range(u_e.shape[2]):
        for t in range(u_e.shape[0]):
            u_medio[i,j] += u_e[t,i,j]/u_e.shape[0]
            v_medio[i,j] += v_e[t,i,j]/v_e.shape[0]
# calculamos coeficientes de regresión lineal  y(t) = a + b*t
u_stl = np.zeros(u_e.shape)
v_stl = np.zeros(u_e.shape)
uwind = np.zeros((u_e.shape[1],u_e.shape[2]))  # en uwind y vwind guardamos los coeficientes b que representan la tendencia
vwind = np.zeros((u_e.shape[1],u_e.shape[2]))
for i in range(u_e.shape[1]):
    print(i)
    for j in range(u_e.shape[2]):
        a = 0
        b = 0
        a1 = 0
        b1 = 0
        for t in range(u_e.shape[0]):
            b += (t-t_medio)*(u_e[t,i,j]-u_medio[i,j])
            b1 += (t-t_medio)*(v_e[t,i,j]-v_medio[i,j])
        b /= den
        b1 /= den
        a = u_medio[i,j]
        a -= b*t_medio
        a1 = v_medio[i,j]
        a1 -= b1*t_medio
        uwind[i,j] = b
        vwind[i,j] = b1
        # descontamos la tendencia lineal
        for t in range(u_e.shape[0]):
            u_stl[t,i,j] = u_e[t,i,j] - b*t
            v_stl[t,i,j] = v_e[t,i,j] - b1*t

####################################################################################################
#####                        Graficamos la tendencia lineal de cada mes                        #####
####################################################################################################


from iris.coords import DimCoord
from iris.cube import Cube

#def main():
    # Create a cube containing the wind speed
latitude = DimCoord(np.linspace(-27.5, -37.5, 5),
                    standard_name='latitude',
                    units='degrees')
longitude = DimCoord(np.linspace(-62.5, -50, 6),
                    standard_name='longitude',
                    units='degrees')
uw = Cube(uwind,
         dim_coords_and_dims=[(latitude, 0),
                              (longitude, 1)],units = 'm s-1')
new_cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
uw.coord(axis='x').coord_system = new_cs
uw.coord(axis='y').coord_system = new_cs
vw = Cube(vwind,
         dim_coords_and_dims=[(latitude, 0),
                              (longitude, 1)],units = 'm s-1')
vw.coord(axis='x').coord_system = new_cs
vw.coord(axis='y').coord_system = new_cs

ulon = uw.coord('longitude')
vlon = vw.coord('longitude')

windspeed = (uw ** 2 + vw ** 2) ** 0.5
windspeed.rename('windspeed')

x = ulon.points
y = uw.coord('latitude').points
u = uw.data
v = vw.data

# Normalise the data for uniform arrow size
u_norm = u / np.sqrt(u ** 2.0 + v ** 2.0)
v_norm = v / np.sqrt(u ** 2.0 + v ** 2.0)

### Toda la región
plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-56.5))

# Get the coordinate reference system used by the data
transform = ulon.coord_system.as_cartopy_projection()
    
# Plot the wind speed as a contour plot
qplt.contourf(windspeed, [0,0.003,0.006,0.009,0.012,0.015,0.018,0.021,0.024,0.027,0.03,0.033,0.036,0.039,0.042,0.045,0.048])

# Add arrows to show the wind vectors
plt.quiver(x, y, u_norm, v_norm, pivot='middle', transform = transform, scale = 12)
    
plt.title('Tendencia lineal de vientos en '+mes)
# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
#states_provinces = cfeat.NaturalEarthFeature(
#    category='cultural',
#    name='admin_1_states_provinces_shp',
    #name='admin_1_states_provinces_lines',
#    scale='10m',
#    facecolor='none', linewidth=0.5)
ax.coastlines()
#ax.gridlines()
#ax.add_feature(states_provinces, edgecolor='black')
ax.add_feature(cfeat.RIVERS)
ax.add_feature(cfeat.BORDERS)
ax.set_extent([297.5, 310, -27.5, -37.5], crs=ccrs.PlateCarree())
ax.set_xticks([300, 305, 310], crs=ccrs.PlateCarree())
ax.set_yticks([-30, -35], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True,
                                   number_format='.0f')
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
#plt.title('prec mm/día', fontsize=12)
plt.savefig('/home/meteo/investigacion/WR.slp/tend_wind_'+mes+'.png', dpi=150)
qplt.show()




###################################################################################3
###################################################################################3
###################################################################################3
###################################################################################3
###################################################################################3
###################################################################################3
###################################################################################3

