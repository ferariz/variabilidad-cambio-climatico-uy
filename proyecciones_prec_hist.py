# -*- coding: utf-8 -*-
"""
Editor de Spyder
Programa para evaluar los modelos del CMIP6 en el período histórico.
En este caso vemos la precipitación en SESA (Southeastern South America) y 
Uruguay en particular.
#"""

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

historico = ['/home/meteo/datos/CMIP5/CCSM4/sesa_pr_Amon_CCSM4_1pctCO2_r1i1p1_194801-200512.nc',
             '/home/meteo/datos/CMIP6/EC_Earth3/sesa_pr_Amon_EC-Earth3_historical_r1i1p1f1_gr_194801-201412.nc',
             '/home/meteo/datos/CMIP6/EC_Earth3-Veg/sesa_pr_Amon_EC-Earth3-Veg_historical_r1i1p1f1_gr_194801-201412.nc',
             '/home/meteo/datos/CMIP6/CAMS-CSM1-0/sesa_pr_Amon_CAMS-CSM1-0_historical_r1i1p1f1_gn_194801-201412.nc',
             '/home/meteo/datos/CMIP6/CanESM5_CCCma/sud_pr_Amon_CanESM5_historical_r1i1p1f1_gn_195001-201412.nc',
             '/home/meteo/datos/CMIP6/BCC_CSM2/sud_pr_Amon_BCC_CSM2_MR_historical_r1i1p1f1_gn_195001_201412.nc',
             '/home/meteo/datos/CMIP6/MIROC6/sud_pr_Amon_MIROC6_historical_r1i1p1f1_gn_195001-201412.nc',
             '/home/meteo/datos/CMIP6/MRI_ESM2/sud_pr_Amon_MRI-ESM2-0_historical_r1i1p1f1_gn_195001-201412.nc']
observado = '/home/meteo/datos/Uy/CRU_pp_sud_1948_2017.nc'
modelos = ['CCSM4','EC_Earth3','EC_Earth3-Veg','CAMS','CanESM5','BCC_CSM2','MIROC6','MRI_ESM2']

####################################################################################################
######### Acomodamos los datos para que tengan el mismo rango temporal y misma unidad ##############
####################################################################################################
pr2 = []
lon = []
lat = []
for i in range(len(modelos)):
    ncin = nc4.Dataset(historico[i], 'r')
    pr2.append(ncin.variables['pr'][:])
    lon.append(ncin.variables['lon'][:])
    lat.append(ncin.variables['lat'][:])
    ncin.close()

pr = []
alcance = min(pr2[0].shape[0],pr2[1].shape[0],pr2[2].shape[0],pr2[3].shape[0])
# Convertir los datos de [kg.m-2.s-1] a [mm], si los datos fueran diarios es solo multiplicar por 86400s y pasa a mm/day
# como los datos vienen de medias mensuales de 30 días (ver info de netcdf) pra que quede en mm/dia hacemos 60*60*24*30=2592000
for i in range(len(modelos)):
    if i < 4:
        pr.append(pr2[i][24:int(alcance),:,:]*60*60*24) # los datos de estos modelos están desde 1948
    else:
        pr.append(pr2[i][0:int(alcance)-24,:,:]*60*60*24) # los datos de estos modelos están desde 1950

years = int((alcance-24)/12)
prr = []
for i in range(len(modelos)):
    prr.append(np.zeros((pr[i].shape)))
    for y in range(years):
        aux = y * 12
        for m in range(12):
            if (m == 0 or m == 2 or m == 4 or m == 6 or m == 7 or m == 9 or m == 11):
                prr[i][aux+m,:,:] = pr[i][aux+m,:,:] * 31
            elif (m == 3 or m == 5 or m == 8 or m == 10):
                prr[i][aux+m,:,:] = pr[i][aux+m,:,:] * 30
            else:
                prr[i][aux+m,:,:] = pr[i][aux+m,:,:] * 28
                
#cargamos las variables del observado 'pre' datos CRU
ncin = nc4.Dataset(observado, 'r')
pre = ncin.variables['pre'][:]
lono = ncin.variables['lon'][:]
lato = ncin.variables['lat'][:]
ncin.close()

pre=pre[24:alcance,:,:]
pre[pre > 10000] = np.nan

####################################################################################################
###############     Definimos funciones para utilizar en el resto del código      ##################
####################################################################################################
def regrillado(obs,mod,lonobs,latobs,lonmod,latmod):
    '''función para regrillar las observaciones (resolución alta) a la grilla del modelo (resolución más baja) para una comparación apropiada'''
    #old grid dim
    loni=lonobs           # 'obs' son las precipitaciones acumuladas observadas CRU
    lati=latobs
    #new grid dim
    lonn=lonmod           # 'mod' son las precipitaciones acumuladas del modelo CMIP
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

def climatologia(prec,latitud,longitud):
    '''función para calcular la climatología de prec(tiempo,lon,lat)'''
    new_prec = prec.reshape(int(prec.shape[0]/12),12,int(latitud.shape[0]),int(longitud.shape[0]))
    return np.mean(new_prec, axis = 0)

def graf_clim(clim,lati,long,model,mes):
    '''función para graficar las climatologías de cada mes'''
    plt.figure()
    meses = ['enero','febrero','marzo','abril','mayo','junio','julio','agosto','setiembre','octubre','noviembre','diciembre']
    clevs = [0,20,40,60,80,100,120,140,160,180,200]
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-56))
    fill = ax.contourf(np.real(long), np.real(lati), np.real(clim), np.real(clevs), transform=ccrs.PlateCarree(), cmap=plt.cm.PuBuGn, extend='both')
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
    plt.savefig('/home/meteo/investigacion/WR.slp/proyecciones/'+model+'/prec/clim_'+model+'%d.png' % (mes), dpi=300)

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
    clevs = [-300, -250, -200, -150, -100, -50, 50, 100, 150, 200, 250, 300]
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-56))
    fill = ax.contourf(np.real(longitud), np.real(latitud), np.real(dif_clim), np.real(clevs), transform=ccrs.PlateCarree(), cmap=plt.cm.BrBG, extend='both')
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
    plt.savefig('/home/meteo/investigacion/WR.slp/proyecciones/'+model+'/prec/dif_clim_'+trimestre+'.png', dpi=300)
    
####################################################################################################
#######   Calculamos las climatologías mensuales de las observaciones y de los modelos   ###########
####################################################################################################
# calculamos la climatología mensual de los datos observados (CRU)
clim_cru = np.zeros((12,len(lato),len(lono)))
clim_cru = climatologia(pre,lato,lono)
for m in range(12):
    graf_clim(clim_cru[m,:,:],lato,lono,'CRU',m)

# calculamos la climatología mensual para los modelos (CMIP5/6)
clim_modelos = []
for i in range(len(modelos)):
    clim_modelos.append(climatologia(prr[i],lat[i],lon[i]))
    for m in range(12):
        graf_clim(clim_modelos[i][m,:,:],lat[i],lon[i],modelos[i],m)

####################################################################################################
#######   Comparamos las climatologías trimestrales entre los modelos y las observaciones    #######
####################################################################################################
# Para una comparación más apropiada, regrillamos las observaciones (mejor resolución) a la grilla de cada modelo.
for i in range(len(modelos)):
    prer = regrillado(pre,prr[i],lono,lato,lon[i],lat[i])   # REGRILLADO
    prer2 = prer.reshape(int(prer.shape[0]/12),12,int(lat[i].shape[0]),int(lon[i].shape[0]))

    prr3 = np.zeros((int(prr[i].shape[0]/12),12,int(lat[i].shape[0]),int(lon[i].shape[0])))
    prr3 = prr[i].reshape(int(prr[i].shape[0]/12),12,int(lat[i].shape[0]),int(lon[i].shape[0]))

    # climatologia trimestrales de datos CRU regrillados
    prer_mam = np.mean(np.squeeze(np.sum(prer2[:,2:5,:,:], axis=1)),axis=0) 
    prer_jja = np.mean(np.squeeze(np.sum(prer2[:,5:8,:,:], axis=1)),axis=0) 
    prer_son = np.mean(np.squeeze(np.sum(prer2[:,8:11,:,:], axis=1)),axis=0) 
    prer_def = np.mean(np.squeeze(np.sum(prer2[:,[0,1,11],:,:], axis=1)),axis=0) 

    # climatologia trimestrales de datos CMIP5
    prr_mam = np.mean(np.squeeze(np.sum(prr3[:,2:5,:,:], axis=1)),axis=0)
    prr_jja = np.mean(np.squeeze(np.sum(prr3[:,5:8,:,:], axis=1)),axis=0)
    prr_son = np.mean(np.squeeze(np.sum(prr3[:,8:11,:,:], axis=1)),axis=0)
    prr_def = np.mean(np.squeeze(np.sum(prr3[:,[0,1,11],:,:], axis=1)),axis=0)

    # calculamos las diferencias y graficamos
    dif_clim_tri(prer_mam,prr_mam,lat[i],lon[i],modelos[i],'MAM')
    dif_clim_tri(prer_jja,prr_jja,lat[i],lon[i],modelos[i],'JJA')
    dif_clim_tri(prer_son,prr_son,lat[i],lon[i],modelos[i],'SON')
    dif_clim_tri(prer_def,prr_def,lat[i],lon[i],modelos[i],'DEF')

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
plt.figure()
plt.plot(clim_cru_uy, label = 'CRU', color= 'black')
plt.plot(clim_mod_uy[0,:], label = modelos[0])
plt.plot(clim_mod_uy[1,:], label = modelos[1])
plt.plot(clim_mod_uy[2,:], label = modelos[2])
plt.plot(clim_mod_uy[3,:], label = modelos[3])
plt.plot(clim_mod_uy[4,:], label = modelos[4])
plt.plot(clim_mod_uy[5,:], label = modelos[5])
plt.plot(clim_mod_uy[6,:], label = modelos[6])
plt.plot(clim_mod_uy[7,:], label = modelos[7])
plt.title('Climatología prec acumulada Uruguay (1950 - 2005)')
plt.ylabel('mm/mes')
plt.legend(fancybox=True, framealpha=0.4)
plt.ylim((30,180))
plt.xticks(ticks=np.arange(12), labels=('Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Set', 'Oct', 'Nov', 'Dic'))
plt.savefig('/home/meteo/investigacion/WR.slp/proyecciones/clim_pre_hist.png', dpi=150)

####################################################################################################
########            Calculamos medidas para evaluar la habilidad de los modelos          ###########
####################################################################################################

# Calculamos la media de la climatología para las observaciones y para los modelos
media_uy_obs = np.mean(clim_cru_uy)
media_uy_mod = np.mean(clim_mod_uy, axis = 1)
sesgo_anual = np.zeros((8)) # Xm - Xobs anual
for i in range(8):
    sesgo_anual[i] = media_uy_mod[i] - media_uy_obs

# sesgo en porcentaje, dividiendo el valor absoluto por la lluvia media observada
sesgo_porcentaje = np.zeros((8))
for i in range(8):
    sesgo_porcentaje[i] = sesgo_anual[i]*100/media_uy_obs

# RMSE
rmse = np.zeros((8))
for i in range(8):
    for t in range(12):
        rmse[i] += (clim_mod_uy[i,t] - clim_cru_uy[t])**2
    rmse[i] = np.sqrt(rmse[i]/12)
rmse[rmse > 10000] = np.nan

# calculo indice VI: es el cociente entre la desviación estándar de cada modelo y la desviación estándar de las observaciones
std_cru = np.std(pre, axis = 0)
std_mod = []
for i in range(len(modelos)):
    std_mod.append(np.std(prr[i],axis = 0))

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
