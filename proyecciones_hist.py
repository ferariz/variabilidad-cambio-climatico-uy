# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal.
"""

####################################################################################################
##### Importamos las bibliotecas necesarias para la adquisición de datos y manejo de variables #####
####################################################################################################
import netCDF4 as nc4
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
import cartopy.feature as cfeat
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import shapefile
from shapely.geometry import shape, Point

historico = ['/home/meteo/datos/CMIP5/CCSM4/sud_tas_Amon_CCSM4_1pctCO2_r1i1p1_194801-200512.nc',
             '/home/meteo/datos/CMIP6/EC_Earth3/sesa_tas_Amon_EC-Earth3_historical_r1i1p1f1_gr_194801-201412.nc',
             '/home/meteo/datos/CMIP6/EC_Earth3-Veg/sesa_tas_Amon_EC-Earth3-Veg_historical_r1i1p1f1_gr_194801-201412.nc',
             '/home/meteo/datos/CMIP6/CAMS-CSM1-0/sesa_tas_Amon_CAMS-CSM1-0_historical_r1i1p1f1_gn_194801-201412.nc',
             '/home/meteo/datos/CMIP6/CanESM5_CCCma/sud_tas_Amon_CanESM5_historical_r1i1p1f1_gn_195001-201412.nc',
             '/home/meteo/datos/CMIP6/BCC_CSM2/sud_tas_Amon_BCC-CSM2-MR_historical_r1i1p1f1_gn_195001-201412.nc',
             '/home/meteo/datos/CMIP6/MIROC6/sud_tas_Amon_MIROC6_historical_r1i1p1f1_gn_195001-201412.nc',
             '/home/meteo/datos/CMIP6/MRI_ESM2/sud_tas_Amon_MRI-ESM2-0_historical_r1i1p1f1_gn_195001-201412.nc']
observado = '/home/meteo/datos/Uy/CRU_tmp_sud_1948_2017.nc'
modelos = ['CCSM4','EC_Earth3','EC_Earth3-Veg','CAMS','CanESM5','BCC_CSM2','MIROC6','MRI_ESM2']

####################################################################################################
############      Acomodamos los datos para que tengan el mismo rango temporal     #################
####################################################################################################
tas2 = []
lon = []
lat = []
for i in range(len(modelos)):
    ncin = nc4.Dataset(historico[i], 'r')
    tas2.append(ncin.variables['tas'][:])
    lon.append(ncin.variables['lon'][:])
    lat.append(ncin.variables['lat'][:])
    ncin.close()

tas = []
alcance = min(tas2[0].shape[0],tas2[1].shape[0],tas2[2].shape[0],tas2[3].shape[0])
for i in range(len(modelos)):
    if i < 4:
        tas.append(tas2[i][24:int(alcance),:,:]) # los datos de estos modelos están desde 1948
    else:
        tas.append(tas2[i][0:int(alcance)-24,:,:]) # los datos de estos modelos están desde 1950
        

#cargamos las variables del observado tmp
ncin = nc4.Dataset(observado, 'r')
tmp = ncin.variables['tmp'][:]
lono = ncin.variables['lon'][:]
lato = ncin.variables['lat'][:]
ncin.close()
tmp=tmp[24:int(alcance),:,:]
tmp[tmp > 1000] = np.nan

####################################################################################################
###############     Definimos funciones para utilizar en el resto del código      ##################
####################################################################################################
def regrillado(obs,mod,lonobs,latobs,lonmod,latmod):
    '''función para regrillar las observaciones (resolución alta) a la grilla del modelo (resolución más baja) para una comparación apropiada'''
    #old grid dim
    loni=lonobs           # 'obs' son las temperaturas observadas CRU
    lati=latobs
    #new grid dim
    lonn=lonmod           # 'mod' son las temperaturas del modelo CMIP
    latn=latmod
    #create mesh
    X, Y = np.meshgrid(loni, lati)
    XI, YI = np.meshgrid(lonn,latn)
    X = X + 360          # esto es porque las observaciones tienen valores de longitud en (-180,180]
    #interp
    A2 = obs.transpose([1,2,0]).reshape(int(len(loni)*len(lati)),mod.shape[0])
    new_A=scipy.interpolate.griddata((X.flatten(),Y.flatten()),A2,(XI,YI),method='linear')
    new_A = new_A.transpose([2,0,1])
    return new_A    # devuelve el arreglo original regrillado según lonmod y latmod.

def climatologia(temp,latitud,longitud):
    '''función para calcular la climatología de temp(tiempo,lon,lat)'''
    new_temp = temp.reshape(int(temp.shape[0]/12),12,int(latitud.shape[0]),int(longitud.shape[0]))
    return np.mean(new_temp, axis = 0)

def graf_clim(clim,lati,long,model,mes):
    '''función para graficar las climatologías de cada mes'''
    plt.figure()
    meses = ['enero','febrero','marzo','abril','mayo','junio','julio','agosto','setiembre','octubre','noviembre','diciembre']
    clevs = [0,3,6,9,12,15,18,21,24,27,30,33]
    #clevs = [-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30]
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-56))
    fill = ax.contourf(np.real(long), np.real(lati), np.real(clim), np.real(clevs), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, extend='both')
    ax.coastlines()
    #ax.gridlines()
    ax.add_feature(cfeat.RIVERS)
    ax.add_feature(cfeat.BORDERS)
    ax.set_extent([min(long), max(long), max(lati), min(lati)], crs=ccrs.PlateCarree())
    ax.set_xticks([295, 300, 305, 310], crs=ccrs.PlateCarree())
    ax.set_yticks([-25, -30, -35], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True,
                                       number_format='.0f')
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    plt.title(model+' climatología '+ meses[mes])
    plt.colorbar(fill, orientation='horizontal')
    #plt.title('prec mm/día', fontsize=12)
    plt.savefig('/home/meteo/investigacion/WR.slp/proyecciones/'+model+'/temp/clim_'+model+'%d.png' % (mes), dpi=300)

def puntos_adentro_Uy(x,y,polygon):
    '''función para evaluar si una coordenada espacial está dentro del territorio uruguayo'''
#    r = shapefile.Reader('/home/meteo/investigacion/WR.slp/URY_adm/URY_adm0.shp')
    # get the shapes
#    shapes = r.shapes()
    # build a shapely polygon from your shape
#    polygon = shape(shapes[0])    
    # build a shapely point from your geopoint
    point = Point(x, y)
        # the contains function does exactly what you want
    return polygon.contains(point)

def dif_clim_tri(obs,mod,latitud,longitud,model,trimestre):
    '''función para calcular la diferencia de climatologías trimestrales'''
    # graficamos la diferencia de las climatologias trimestrales
    dif_clim = mod - obs
    plt.figure()
    clevs = [-6,-5,-4,-3,-2,-1, 0,1,2,3,4,5,6]
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-56))
    fill = ax.contourf(np.real(longitud), np.real(latitud), np.real(dif_clim), np.real(clevs), transform=ccrs.PlateCarree(), cmap=plt.cm.coolwarm, extend='both')
    ax.coastlines()
    #ax.gridlines()
    ax.add_feature(cfeat.RIVERS)
    ax.add_feature(cfeat.BORDERS)
    ax.set_extent([min(longitud), max(longitud), max(latitud), min(latitud)], crs=ccrs.PlateCarree())
    ax.set_xticks([295, 300, 305, 310], crs=ccrs.PlateCarree())
    ax.set_yticks([-25, -30, -35], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True,
                                       number_format='.0f')
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    plt.title(trimestre +' ('+model+' - CRU)')
    plt.colorbar(fill, orientation='horizontal')
    plt.savefig('/home/meteo/investigacion/WR.slp/proyecciones/'+model+'/temp/dif_clim_'+trimestre+'.png', dpi=300)
    
####################################################################################################
#######   Calculamos las climatologías mensuales de las observaciones y de los modelos   ###########
####################################################################################################
# calculamos la climatología mensual de los datos observados (CRU)
clim_cru = np.zeros((12,len(lato),len(lono)))
clim_cru = climatologia(tmp,lato,lono)
for m in range(12):
    graf_clim(clim_cru[m,:,:],lato,lono,'CRU',m)

# calculamos la climatología mensual para los modelos (CMIP5/6)
clim_modelos = []
for i in range(len(modelos)):
    clim_modelos.append(climatologia(tas[i] - 273.15,lat[i],lon[i]))
    for m in range(12):
        graf_clim(clim_modelos[i][m,:,:],lat[i],lon[i],modelos[i],m)

####################################################################################################
#######   Comparamos las climatologías trimestrales entre los modelos y las observaciones    #######
####################################################################################################
# Para una comparación más apropiada, regrillamos las observaciones (mejor resolución) a la grilla de cada modelo.
for i in range(len(modelos)):
    tmpr = regrillado(tmp,tas[i],lono,lato,lon[i],lat[i])   # REGRILLADO
    tmpr2 = tmpr.reshape(int(tmpr.shape[0]/12),12,int(lat[i].shape[0]),int(lon[i].shape[0]))

    tas3 = np.zeros((int(tas[i].shape[0]/12),12,int(lat[i].shape[0]),int(lon[i].shape[0])))
    tas3 = tas[i].reshape(int(tas[i].shape[0]/12),12,int(lat[i].shape[0]),int(lon[i].shape[0])) - 273.15

    # climatologia trimestrales de datos CRU regrillados
    tmpr_mam = np.mean(np.squeeze(np.mean(tmpr2[:,2:5,:,:], axis=1)),axis=0) 
    tmpr_jja = np.mean(np.squeeze(np.mean(tmpr2[:,5:8,:,:], axis=1)),axis=0) 
    tmpr_son = np.mean(np.squeeze(np.mean(tmpr2[:,8:11,:,:], axis=1)),axis=0) 
    tmpr_def = np.mean(np.squeeze(np.mean(tmpr2[:,[0,1,11],:,:], axis=1)),axis=0) 

    # climatologia trimestrales de datos CMIP5
    tas_mam = np.mean(np.squeeze(np.mean(tas3[:,2:5,:,:], axis=1)),axis=0)
    tas_jja = np.mean(np.squeeze(np.mean(tas3[:,5:8,:,:], axis=1)),axis=0)
    tas_son = np.mean(np.squeeze(np.mean(tas3[:,8:11,:,:], axis=1)),axis=0)
    tas_def = np.mean(np.squeeze(np.mean(tas3[:,[0,1,11],:,:], axis=1)),axis=0)

    # calculamos las diferencias y graficamos
    dif_clim_tri(tmpr_mam,tas_mam,lat[i],lon[i],modelos[i],'MAM')
    dif_clim_tri(tmpr_jja,tas_jja,lat[i],lon[i],modelos[i],'JJA')
    dif_clim_tri(tmpr_son,tas_son,lat[i],lon[i],modelos[i],'SON')
    dif_clim_tri(tmpr_def,tas_def,lat[i],lon[i],modelos[i],'DEF')

####################################################################################################
#####   Climatología de las obs y los modelos considerando un promedio espacial dentro de UY   #####
####################################################################################################
r = shapefile.Reader('/home/meteo/investigacion/WR.slp/URY_adm/URY_adm0.shp')   # carga el shapefile de Uruguay
shapes = r.shapes()                                                             # get the shapes
polygon = shape(shapes[0])                                                      # build a shapely polygon from your shape
clim_cru_uy = np.zeros((12))    # observaciones
for m in range(12):
    cont = 0
    for i in range(len(lato)):
        for j in range(len(lono)):
            if puntos_adentro_Uy(lono[j],lato[i],polygon):   # función definida arriba, evalúa si lono[j],lato[i] está dentro de Uy
                clim_cru_uy[m] += clim_cru[m,i,j]
                cont += 1
    print(m,cont)
    clim_cru_uy[m] /= cont
    
clim_mod_uy = np.zeros((8,12))    # modelos
for k in range(len(modelos)):
    for m in range(12):
        cont = 0
        for i in range(len(lat[k])):
            for j in range(len(lon[k])):
                if puntos_adentro_Uy(lon[k][j]-360,lat[k][i],polygon):
                    clim_mod_uy[k,m] += clim_modelos[k][m,i,j]
                    cont += 1
        clim_mod_uy[k,m] /= cont
    print(k,cont)
# graficamos las climatologías
plt.plot(clim_cru_uy, label = 'CRU', color= 'black')
plt.plot(clim_mod_uy[0,:], label = modelos[0])
plt.plot(clim_mod_uy[1,:], label = modelos[1])
plt.plot(clim_mod_uy[2,:], label = modelos[2])
plt.plot(clim_mod_uy[3,:], label = modelos[3])
plt.plot(clim_mod_uy[4,:], label = modelos[4])
plt.plot(clim_mod_uy[5,:], label = modelos[5])
plt.plot(clim_mod_uy[6,:], label = modelos[6])
plt.plot(clim_mod_uy[7,:], label = modelos[7])
plt.title('Climatología temperatura media Uruguay (1950 - 2005)')
plt.ylabel('$^\circ $C')
plt.legend(fancybox=True, framealpha=0.4)
plt.ylim((5,35))
plt.xticks(ticks=np.arange(12), labels=('Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Set', 'Oct', 'Nov', 'Dic'))
plt.savefig('/home/meteo/investigacion/WR.slp/proyecciones/clim_tas_hist.png', dpi=150)

####################################################################################################
########            Calculamos medidas para evaluar la habilidad de los modelos          ###########
####################################################################################################

# Calculamos la media de la climatología para las observaciones y para los modelos
media_uy_obs = np.mean(clim_cru_uy)
media_uy_mod = np.mean(clim_mod_uy, axis = 1)
sesgo_anual = np.zeros((8)) # Xm - Xobs anual
for i in range(8):
    sesgo_anual[i] = media_uy_mod[i] - media_uy_obs

# RMSE
rmse = np.zeros((8))
for i in range(8):
    for t in range(12):
        rmse[i] += (clim_mod_uy[i,t] - clim_cru_uy[t])**2
    rmse[i] = np.sqrt(rmse[i]/12)
rmse[rmse > 100] = np.nan

# calculo indice VI: es el cociente entre la desviación estándar de cada modelo y la desviación estándar de las observaciones
std_cru = np.std(tmp, axis = 0)
std_mod = []
for i in range(len(modelos)):
    std_mod.append(np.std(tas[i],axis = 0))

std_cru_uy = 0    # observaciones
cont = 0
for i in range(len(lato)):
    for j in range(len(lono)):
        if puntos_adentro_Uy(lono[j],lato[i],polygon):
            std_cru_uy += std_cru[i,j]
            cont += 1
std_cru_uy /= cont
    
std_mod_uy = np.zeros((8))    # modelos
for k in range(len(modelos)):
    cont = 0
    for i in range(len(lat[k])):
        for j in range(len(lon[k])):
            if puntos_adentro_Uy(lon[k][j] - 360,lat[k][i],polygon):
                std_mod_uy[k] += std_mod[k][i,j]
                cont += 1
    std_mod_uy[k] /= cont

VI = np.zeros((8))
VI = std_mod_uy/std_cru_uy
