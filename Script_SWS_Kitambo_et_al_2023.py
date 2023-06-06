
#*********************************************************************************
#*******************    Surface Water Storage calculation    *********************
#*********************************************************************************
#*	 Laboratoire d’Etudes en Géophysique et Océanographie Spatiales (LEGOS)      *
#*	            Université Toulouse III - Paul Sabatier (UPS3)                   *
#*                  	           June 2023                                     *
#*                  	          Version 1.0                                    *
#*                              Benjamin Kitambo                                 *
#*********************************************************************************
#*********************************************************************************
#
#  SCRIPT Script_SWS_Kitambo_et_al_2023 (Script_SWS_Kitambo_et_al_2023.py) is the
#  only routine for this script.
#
#---------------------------------------------------------------------------------
#  Opening:
#
#    This routine was developed by Benjamin Kitambo and Fabrice Papa during the
#    research thesis of Benjamin Kitambo in early 2022, in Toulouse-France.
#
#    The reference of the publication in which the script was used is given below:
#
#   Kitambo, B. M., Papa, F., Paris, A., Tshimanga, R. M., Frappart, F., Calmant, S.,
#   Elmi, O., Fleischmann, A. S., Becker, M., Tourian, M. J., Jucá Oliveira, R. A.,
#   and Wongchuig, S.: A long-term monthly surface water storage dataset for the Congo
#   basin from 1992 to 2015, Earth Syst. Sci. Data Discuss. [preprint],
#   https://doi.org/10.5194/essd-2022-376, in review, 2022.
#---------------------------------------------------------------------------------
#
#---------------------------------------------------------------------------------
#  Script_SWS_Kitambo_et_al_2023.py
#---------------------------------------------------------------------------------
#  Discussion:
# 
#    This script consists of code to run the reproducibility of the hypsometric curve
#    approach dataset of surface water storage. The input data includes the digital
#    elevation model, latitude-longitude-width of each Global Inundation Exent from
#    Multi-satellite (GIEMS-2) pixel dataset, and the surface water extent from
#    GIEMS-2 dataset itself.
#
#  Usage:
#
#    uses libraries
#
#    Reads (Input data)
#
#    * Mosaic_f_AF.tif, stands for the digital elevation model.
#    * cell_lat_lon_width_Africa.csv, stands for the latitude-longitude-width of each
#      GIEMS-2 pixel dataset.
#    * Africa2.csv, stands for GIEMS-2 dataset.
#
#    Creates (Output)
#
#    * hypsoCurve, stands for the hypsometric curve file in csv format.
#    * df_hypso_corr, stands for the corrected hypsometric curve file in csv format.
#    * df_WaterStorage, stands for the hypsometric curve providing the surface water
#      extent area-storage relationship file in csv format.
#    * sws, stands for the monthly surface water storage variations in csv format.
#
#---------------------------------------------------------------------------------
#  Licensing:
#
#    This code is distributed under the GNU/GLP Licence
#
#  Version/Modified: 
#
#    06 June 2023
#    By: Benjamin Kitambo
#
#  Authors:
#
#    Original Python version by:
#    * Benjamin Kitambo
#    * Fabrice Papa
#    Present Python version by:
#    * Benjamin Kitambo
#
#  Main Reference:
#
#   Kitambo, B. M., Papa, F., Paris, A., Tshimanga, R. M., Frappart, F., Calmant, S.,
#   Elmi, O., Fleischmann, A. S., Becker, M., Tourian, M. J., Jucá Oliveira, R. A.,
#   and Wongchuig, S.: A long-term monthly surface water storage dataset for the Congo
#   basin from 1992 to 2015, Earth Syst. Sci. Data Discuss. [preprint],
#   https://doi.org/10.5194/essd-2022-376, in review, 2022.
#---------------------------------------------------------------------------------
#  Variables and Parameters:
#  *Variables declarations and procedures are all commented below.
#---------------------------------------------------------------------------------
# End of header
#---------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import geopandas
import numpy as np
import pandas as pd
import os
from osgeo import gdal
from geopandas.tools import sjoin
import sys
import csv
from scipy import interpolate
from itertools import accumulate

PATH_IN = r"D:\thesis_Kitambo\global_dataset"
PATH_IN1 = r"D:\FABDEM_AF\Output\Mosaic_f_AF.tif"
PATH_IN3 = r"D:\thesis_Kitambo\global_dataset\Africa"
env_path = r"C:\Users\Kitambo\PycharmProjects\CRB_SV_GIEMS\\"

##################Computation of the Hypsometric curve
list_nbr_cell=[]
list_lat=[]
list_lon=[]
list_eleva=[]
list_extract_tiles=[]

inputDEM = gdal.Open(PATH_IN1)
cell_width = pd.read_csv(os.path.join(PATH_IN, "cell_lat_lon_width_Africa.csv"),sep=",",index_col=False, header=None)

for (col_name, col_data) in cell_width.T.iteritems():
    list_data = np.asarray(col_data)
    df = pd.DataFrame()
    df = df.assign(array=list_data)

    ## GIEMS node (clipping extent)
    cell_dem = str(int(df['array'][0]))
    y1_4 = str('{:.2f}'.format(df['array'][1] + 0.125))
    y2_3 = str('{:.2f}'.format(df['array'][1] - 0.125))
    x1_2 = str('{:.3f}'.format(df['array'][2] - (df['array'][3] / 2)))
    x3_4 = str('{:.3f}'.format(df['array'][2] + (df['array'][3] / 2)))
    xmin, ymin, xmax, ymax = x1_2, y2_3, x3_4, y1_4

    ##Clip according to cell GIEMS extent
    clip = gdal.Warp(cell_dem + "_" + "_clip.tif", inputDEM, outputBounds = (xmin, ymin, xmax, ymax),
                     cropToCutline= True, format="GTiff", dstNodata = np.nan)
    
    ##Resample to 90m x 90m
    xRes= 0.0008333333333333328447
    yRes= -0.0008333333333333330616
    resa =  gdal.Warp(cell_dem + "_" + "_resa.tif", clip, xRes = xRes, yRes= yRes)

    ##Build hypsometric curve
    dem = np.array(resa.GetRasterBand(1).ReadAsArray()) ##Convert tiff in array  ##list_extract_tiles[0]
    order= np.sort(dem, axis=None)  ##Sort the flattened array
    dt = pd.DataFrame(order) ##Convert np.array into dataframe
    dt.columns = ['eleva']
    dt = dt.loc[dt['eleva'] >= 0,:]
    rem_nan= dt.dropna() ##Remove nan in the dataframe
    order_perce= np.array_split(rem_nan.to_numpy(), 101) ##Split sorted array altitude in percentage (100) of the length
    split_avg = [np.median(arr) for arr in order_perce] ##Get the average of each sub-array using median

    ##Built of dataframe
    list_nbr_cell.append(cell_dem)
    list_lat.append('{:.3f}'.format(df['array'][1]))
    list_lon.append('{:.3f}'.format(df['array'][2]))
    list_eleva.append(split_avg)

    ##Clear memory
    del resa
    del clip
    list_extract_tiles.clear()
    ##Delete created file in environment work directory
    for items in os.listdir(env_path):
        if items.endswith("tif"):
            ##Construct full file path
            itpath = env_path + items
            if os.path.isfile(itpath):
                os.remove(itpath)
dfina = pd.DataFrame(list_eleva)
dfina = dfina.assign(cell_nbr=list_nbr_cell, lat=list_lat, lon=list_lon)
cols = list(dfina.columns)
cols = cols[101:104] + cols[:101]
hypsoCurve = dfina[cols]


##################Correction of the hypsometric curve
xline= [ i for i in range(0, 101, 1)] # series of percentage value
cell_width.columns = ['cell_nbr', 'lat', 'lon', 'width']

lst_lat = []
lst_lon = []
lst_cell_nbr = []
lst_AM = []
lst_AMC = []
list_eleva = []
lst_const = []
lst_repl = []

for cell_nbr in hypsoCurve.iloc[:,0:1].values.tolist():
    lst_cell_nbr.append(cell_nbr[0])
    sub_dt = cell_width.loc[cell_width['cell_nbr'] == cell_nbr[0]]
    lat = sub_dt.iloc[0]['lat']
    lon = sub_dt.iloc[0]['lon']
    lst_lat.append(lat)
    lst_lon.append(lon)

    ## Dataframe of altitude
    hypso_speci = hypsoCurve[hypsoCurve.iloc[:,0]== cell_nbr[0]]
    hypso_speci_T= hypso_speci.T.iloc[3:]
    hypso_speci_T['x'] = xline
    hypso_speci_T.columns= ['alt', 'x'] ##Change alt in distor_alt for ALOS, ASTER, TANDEM-X
    hypso_speci_T = hypso_speci_T.reset_index() ####Reorder index values.

    ####Smooth ALOS, ASTER and TANDEM-X by Savitzky-Golay filter
    # yhat = savgol_filter(hypso_speci_T['distor_alt'], 21, 3)  #window size 21 (21 good of not distorting the signal tendency), polynomial order 3
    # hypso_speci_T['alt'] = yhat.round(2)

    lst_ele = []
    lst_pred = []
    hypso_speci_np = hypso_speci_T['alt'].to_numpy()
    area_perce = hypso_speci_T['x'].to_numpy()
    hypso_speci_np2 = hypso_speci_np.copy()

    x = np.arange(0,hypso_speci_np.size)
    for i in area_perce:
        alt = hypso_speci_np[i : i + 5]
        alt2 = hypso_speci_np[i - 2: i]
        std = round(np.std(alt),2)
        std2 = round(np.std(alt2),2)

        tresholdSTD = 0.7  ##Value of STD to be changed according to the amplitude from altimetry-derived surface water level
        if std > tresholdSTD:
            if std2 < tresholdSTD:
                f = interpolate.interp1d(x[i-2:i], hypso_speci_np2[i-2:i], fill_value="extrapolate")
                hypso_speci_np2[i] = f(x[i])
            else:
                hypso_speci_np2[i] = np.nan
    hypso_speci_np2[-1] = np.nan ##Set the last value to nan
    df = pd.DataFrame(hypso_speci_np2)
    df.columns = ["alt"]
    count_nan = int(((df.isna().sum()[0]) * 100) / len(df["alt"].tolist()))
    if count_nan == 99:
        df['alt'] = np.nan
    df2 = df.dropna(axis= 0)
    index = df2.index
    valist = np.asarray(df2["alt"].tolist())
    valist2 = valist.copy()
    for j in np.arange(1,valist.size):  ##Cross check for higher denivelation
        dif = valist[j] - valist[j - 1]
        if dif > tresholdSTD:
            const = dif ##Constant representing the denivelation to remove
            idx = j
            valist2[idx:] = valist2[idx:] - const
        else:
            const = 0
            idx = 0
    ####Get the index of nan values as list of tuples
    nanlist = list(map(tuple, np.where(np.isnan(hypso_speci_np2))))
    nanlist = np.asarray(nanlist)[0]
    # print(nanlist)
    if all(pd.isna(df['alt'])) == False: ###condition for entire nan dataframe
        fct = interpolate.interp1d(np.asarray(index), valist2, fill_value="extrapolate")
        ynew = fct(nanlist)
    else:
        ynew = np.empty((101, 1,)) * np.nan ##create range of nan

    df2 = df2.replace(df2["alt"].tolist(), list(valist2))
    dfval = pd.DataFrame(data=ynew, index=nanlist)
    dfcorr = pd.concat([df2, dfval], axis=1)
    dfcorr['inter'] = dfcorr.sum(axis=1)
    dfcorr = dfcorr.replace(0,np.nan)
    dfcorr['sur'] = xline
    list_eleva.append(dfcorr['inter'].tolist())

##Create dataframe of hypsometric curve corrected dataset
df = pd.DataFrame(list_eleva)
df = df.round(2)
df = df.assign(cell_nbr=lst_cell_nbr, lat=lst_lat, lon=lst_lon)
cols = list(df.columns)
cols = cols[101:104] + cols[:101]
df_hypso_corr = df[cols]

###########################Remove all artifacts (NaN values).
df_hypso_corr['mean'] = df_hypso_corr.iloc[:, 3:288].mean(axis = 1)
df_hypso_corr['mean'].replace(np.nan,0 , inplace = True)
df_hypso_corr.loc[(df_hypso_corr['mean'].astype(float) == 0), :] = np.nan ### Pixel with NaN value
df_hypso_corr.dropna(axis= 0, how='all', inplace = True)
df_hypso_corr = df_hypso_corr.loc[df['mean'] >= 0,:] ###Negative value
df_hypso_corr2 = df_hypso_corr.iloc[:,:-1]


##################Computation of the area – surface water storage relationship
giems = pd.read_csv(os.path.join(PATH_IN3, "Africa2.csv"),sep=",", header=0, index_col=None) # GIEMS dataset
date=pd.date_range('1992', freq='MS', periods=288)

lst_cell_nbr = []
lst_v_i = []
lst_lat = []
lst_lon = []
lst_AM = []
lst_vol = []
lst_vol2 = []

for cell_nbr in df_hypso_corr2.iloc[:,0].tolist():
    lst_cell_nbr.append(cell_nbr)
    ss_giems = giems.loc[giems['cell_x'] == cell_nbr]
    lat = ss_giems.iloc[0]['lat_x']
    lon = ss_giems.iloc[0]['lon_x']
    lst_lat.append(lat)
    lst_lon.append(lon)
    ss_giems =ss_giems.iloc[:, :-3]

    hypso_spec = df_hypso_corr2[df_hypso_corr2.iloc[:,0]== cell_nbr]
    hypso_spec_T= hypso_spec.T.iloc[3:]
    hypso_spec_T['sur'] = xline
    hypso_spec_T.columns= ["alt", "sur"]

    ##Volume
    for i in xline:
        alt_max = hypso_spec_T.loc[hypso_spec_T['sur'] == xline[i]]['alt'].values[0]
        alt_min = hypso_spec_T.loc[hypso_spec_T['sur'] == xline[i - 1]]['alt'].values[0]
        ##Formula representing a volume of a trunk or regular truncated pyramid
        vol = (((alt_max - alt_min) / 1000) * ((773 / 100) * (i + (i - 1) + np.sqrt(i * (i - 1))))) / 3
        ##Value 773 represents the maximum surface extent of GIEMS-2 dataset
        lst_vol.append(vol)
    cum = lst_vol[-101:] ## Consider only the last list of the corresponding cell
    cum[0] = 0.0  ## Replace the negative 0 by positive 0
    cum = list(accumulate(cum))  ## Cumulate the volume for each percentage of surface inondated
    lst_vol2.append(cum)

### Calculation of the average volume amplitude
df = pd.DataFrame()
df = df.assign(vol = cum)
df['sur'] = xline

df = pd.DataFrame(lst_vol2)
df = df.assign(cell_nbr=lst_cell_nbr, lat=lst_lat, lon=lst_lon)
cols = list(df.columns)
cols = cols[101:104] + cols[0:101]
df_WaterStorage = df[cols]

##################Computation of the time series of surface water storage variations
lst_cell_nbr = []
lst_lat = []
lst_lon = []
lst_stor = []

for cell_nbr in df_WaterStorage['cell_nbr'].tolist():
    lst_cell_nbr.append(cell_nbr)
    vol_S = df_WaterStorage.loc[df_WaterStorage['cell_nbr']== cell_nbr]
    vol_T= vol_S.T.iloc[3:104] ### for 3:104 because ignoring coordinate and amplitude columns
    vol_T['sur'] = xline
    vol_T.columns= ["vol", "sur"]

    ss_giems = giems.loc[giems['cell_x'] == cell_nbr]
    lat = ss_giems.iloc[0]['lat_x']
    lon = ss_giems.iloc[0]['lon_x']
    lst_lat.append(lat)
    lst_lon.append(lon)
    ss_giems =ss_giems.iloc[:, :-3]
    min_ss_giems = ss_giems.min(axis = 1)

    for sur in ss_giems.T.iloc[:,0].tolist():
        area = int(round((sur * 100)/ 773))
        area_min = int(round((min_ss_giems.iloc[0] * 100)/ 773))
        storage = vol_T.loc[vol_T['sur'] == area]['vol'].values[0]
        stora_min = vol_T.loc[vol_T['sur'] == area_min]['vol'].values[0]
        surf_stora = storage - stora_min
        lst_stor.append(surf_stora)
np = np.array_split(lst_stor, len(lst_cell_nbr))
df = pd.DataFrame(np)

df = df.assign(cell_nbr=lst_cell_nbr, lat=lst_lat, lon=lst_lon)
cols = list(df.columns)
cols = cols[288:292] + cols[0:288]
sws = df[cols]
##########################################################



