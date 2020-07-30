#Import necesary modules
#Use shorter names (np, pd, plt) instead of full (numpy, pandas, matplotlib.pylot) for convenience
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd
import h5py
import xarray as xr
import numpy as np
import pdb
import numpy.ma as ma

def getSnowandConverttoThickness(dF, snowDepthVar='snowDepth', 
                                 snowDensityVar='snowDensity',
                                 outVar='iceThickness'):
    """ Grid using nearest neighbour the NESOSIM snow depths to the 
    high-res ICESat-1 freeboard locations
    """
    
    # Convert freeboard to thickness
    # Need to copy arrays or it will overwrite the pandas column!
    freeboardT=np.copy(dF['freeboard'].values)
    snowDepthT=np.copy(dF[snowDepthVar].values)
    snowDensityT=np.copy(dF[snowDensityVar].values)
    ice_thickness = freeboard_to_thickness(freeboardT, snowDepthT, snowDensityT)
    #print(ice_thickness)
    dF[outVar] = pd.Series(np.array(ice_thickness), index=dF.index)
   
    return dF

def freeboard_to_thickness(freeboardT, snow_depthT, snow_densityT):
    """
    Hydrostatic equilibrium equation to calculate sea ice thickness 
    from freeboard and snow depth/density data

    Args:
        freeboardT (var): ice freeboard
        snow_depthT (var): snow depth
        snow_densityT (var): final snow density

    Returns:
        ice_thicknessT (var): ice thickness dereived using hydrostatic equilibrium

    """

    # Define density values
    rho_w=1024.
    rho_i=925.
    #rho_s=300.

    # set snow to freeboard where it's bigger than freeboard.
    snow_depthT[snow_depthT>freeboardT]=freeboardT[snow_depthT>freeboardT]

    ice_thicknessT = (rho_w/(rho_w-rho_i))*freeboardT - ((rho_w-snow_densityT)/(rho_w-rho_i))*snow_depthT

    return ice_thicknessT

def getWarrenData(dF, outSnowVar, outDensityVar='None'):
    """
    Assign Warren1999 snow dept/density climatology to dataframe

    Added 

    Args:
        dF (data frame): Pandas dataframe
        outSnowVar (string): name of Warren snow depth variable
        outDensityVar (string): name of Warren snow density variable
        

    Returns:
        dF (data frame): Pandas dataframe updated to include colocated Warren snow depth and density
        
    """
    
    # Generate empty lists
    snowDepthW99s=ma.masked_all(np.size(dF['freeboard'].values))
    if (outDensityVar!='None'):
        snowDensityW99s=ma.masked_all(np.size(dF['freeboard'].values))

    # Loop over all freeboard values (rows)
    for x in range(np.size(dF['freeboard'].values)):
        #print(x, dF['lon'].iloc[x], dF['lat'].iloc[x], dF['month'].iloc[x]-1)
        # SUbtract 1 from month as warren index in fucntion starts at 0
        snowDepthDayW99T, snowDensityW99T=WarrenClimatology(dF['lon'].iloc[x], dF['lat'].iloc[x], dF['datetime'].iloc[x].month-1)
        

        # Append values to list
        snowDepthW99s[x]=snowDepthDayW99T
        if (outDensityVar!='None'):
            snowDensityW99s[x]=snowDensityW99T

    # Assign list to dataframe as a series
    dF[outSnowVar] = pd.Series(snowDepthW99s, index=dF.index)
    if (outDensityVar!='None'):
        dF[outDensityVar] = pd.Series(snowDensityW99s, index=dF.index)
    

    return dF

def WarrenClimatology(lonT, latT, monthT):
    """
    Get Warren1999 snow depth climatology

    Args:
        lonT (var): longitude
        latT (var): latitude
        monthT (var): month with the index starting at 0
        
    Returns:
        Hs (var): Snow depth (m)
        rho_s (var): Snow density (kg/m^3)
        
    """

    H_0 = [28.01, 30.28, 33.89, 36.8, 36.93, 36.59, 11.02, 4.64, 15.81, 22.66, 25.57, 26.67]
    a = [.127, .1056, .5486, .4046, .0214, .7021, .3008, .31, .2119, .3594, .1496, -0.1876]
    b = [-1.1833, -0.5908, -0.1996, -0.4005, -1.1795, -1.4819, -1.2591, -0.635, -1.0292, -1.3483, -1.4643, -1.4229]
    c = [-0.1164, -0.0263, 0.0280, 0.0256, -0.1076, -0.1195, -0.0811, -0.0655, -0.0868, -0.1063, -0.1409, -0.1413]
    d = [-0.0051, -0.0049, 0.0216, 0.0024, -0.0244, -0.0009, -0.0043, 0.0059, -0.0177, 0.0051, -0.0079, -0.0316]
    e = [0.0243, 0.0044, -0.0176, -0.0641, -0.0142, -0.0603, -0.0959, -0.0005, -0.0723, -0.0577, -0.0258, -0.0029]

    # Convert lat and lon into degrees of arc, +x axis along 0 degrees longitude and +y axis along 90E longitude
    x = (90.0 - latT)*np.cos(lonT * np.pi/180.0)  
    y = (90.0 - latT)*np.sin(lonT*np.pi/180.0) 

    Hs = H_0[monthT] + a[monthT]*x + b[monthT]*y + c[monthT]*x*y + (d[monthT]*x*x) + (e[monthT]*y*y)
    

    # Now get SWE, although this is not returned by the function

    H_0swe = [8.37, 9.43,10.74,11.67,11.8,12.48,4.01,1.08,3.84,6.24,7.54,8.0]
    aswe = [-0.027,0.0058,0.1618,0.0841,-0.0043,0.2084,0.097,0.0712,0.0393,0.1158,0.0567,-0.054]
    bswe = [-0.34,-0.1309,0.0276,-0.1328,-0.4284,-0.5739,-0.493,-0.145,-0.2107,-0.2803,-0.3201,-0.365]
    cswe = [-0.0319,0.0017,0.0213,0.0081,-0.038,-0.0468,-0.0333,-0.0155,-0.0182,-0.0215,-0.0284,-0.0362]
    dswe = [-0.0056,-0.0021,0.0076,-0.0003,-0.0071,-0.0023,-0.0026,0.0014,-0.0053,0.0015,-0.0032,-0.0112]
    eswe = [-0.0005,-0.0072,-0.0125,-0.0301,-0.0063,-0.0253,-0.0343,0,-0.019,-0.0176,-0.0129,-0.0035]


    swe = H_0swe[monthT] + aswe[monthT]*x + bswe[monthT]*y + cswe[monthT]*x*y + dswe[monthT]*x*x + eswe[monthT]*y*y

    # Density in kg/m^3
    rho_s = 1000.*(swe/Hs)  
    #print(ma.mean(rho_s))

    # Could mask out bad regions (i.e. land) here if desired.
    # Hsw[where(region_maskG<9.6)]=np.nan
    # Hsw[where(region_maskG==14)]=np.nan
    # Hsw[where(region_maskG>15.5)]=np.nan

    # Could mask out bad regions (i.e. land) here if desired.
    #rho_s[where(region_maskG<9.6)]=np.nan
    #rho_s[where(region_maskG==14)]=np.nan
    #rho_s[where(region_maskG>15.5)]=np.nan

    # Convert snow depth to meters
    Hs=Hs/100.

    return Hs, rho_s

def get_psnlatslons(data_path, res=25):
    """ Get NSIDC polar stereographic grid data"""
    
    if (res==25):
        # 25 km grid
        mask_latf = open(data_path+'/psn25lats_v3.dat', 'rb')
        mask_lonf = open(data_path+'/psn25lons_v3.dat', 'rb')
        lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
        lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])
    elif (res==12):
        # 12.5 km grid
        mask_latf = open(data_path+'/psn12lats_v3.dat', 'rb')
        mask_lonf = open(data_path+'/psn12lons_v3.dat', 'rb')
        lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [896, 608])
        lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [896, 608])
    elif (res==6):
        # 12.5 km grid
        mask_latf = open(data_path+'/psn06lats_v3.dat', 'rb')
        mask_lonf = open(data_path+'/psn06lons_v3.dat', 'rb')
        lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [1792, 1216])
        lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [1792, 1216])

    return lats_mask, lons_mask




def assignRegionMask(dF, mapProj, ancDataPath='../Data/'):
    """
    Grab the NSIDC region mask and add to dataframe as a new column

    # 1   non-region oceans
    # 2   Sea of Okhotsk and Japan
    # 3   Bering Sea
    # 4   Hudson Bay
    # 5   Gulf of St. Lawrence
    # 6   Baffin Bay/Davis Strait/Labrador Sea
    # 7   Greenland Sea
    # 8   Barents Seas
    # 9   Kara
    # 10   Laptev
    # 11   E. Siberian
    # 12   Chukchi
    # 13   Beaufort
    # 14   Canadian Archipelago
    # 15   Arctic Ocean
    # 20   Land
    # 21   Coast

    Args:
        dF (data frame): original data frame
        mapProj (basemap instance): basemap map projection
          
    Returns:
        dF (data frame): data frame including ice type column (1 = multiyear ice, 0 = everything else)

    """
    region_mask, xptsI, yptsI = get_region_mask_sect(ancDataPath, mapProj, xypts_return=1)

    xptsI=xptsI.flatten()
    yptsI=yptsI.flatten()
    region_mask=region_mask.flatten()

    #iceTypeGs=[]
    regionFlags=ma.masked_all((size(dF['freeboard'].values)))
    for x in range(size(dF['freeboard'].values)):
        # Find nearest ice type
        dist=sqrt((xptsI-dF['xpts'].iloc[x])**2+(yptsI-dF['ypts'].iloc[x])**2)
        index_min = np.argmin(dist)
        regionFlags[x]=int(region_mask[index_min])
        
        # This is what I sometimes do but it appears slower in this case..
        # I checked and they gave the same answers
        # iceTypeG2 = griddata((xpts_type, ypts_type), ice_typeT2, (dF['xpts'].iloc[x], dF['ypts'].iloc[x]), method='nearest') 
        # print(iceTypeG)
        # iceTypeGs.append(iceTypeG)

    dF['region_flag'] = pd.Series(regionFlags, index=dF.index)

    return dF

def get_region_mask_sect(datapath, mplot, xypts_return=0):
    """ Get NSIDC section mask data """
    
    datatype='uint8'
    file_mask = datapath+'/sect_fixed_n.msk'
    # 1   non-region oceans
    # 2   Sea of Okhotsk and Japan
    # 3   Bering Sea
    # 4   Hudson Bay
    # 5   Gulf of St. Lawrence
    # 6   Baffin Bay/Davis Strait/Labrador Sea
    # 7   Greenland Sea
    # 8   Barents Seas
    # 9   Kara
    # 10   Laptev
    # 11   E. Siberian
    # 12   Chukchi
    # 13   Beaufort
    # 14   Canadian Archipelago
    # 15   Arctic Ocean
    # 20   Land
    # 21   Coast
    fd = open(file_mask, 'rb')
    region_mask = fromfile(file=fd, dtype=datatype)
    region_mask = reshape(region_mask, [448, 304])

    #xpts, ypts = mplot(lons_mask, lats_mask)
    if (xypts_return==1):
        mask_latf = open(datapath+'/psn25lats_v3.dat', 'rb')
        mask_lonf = open(datapath+'/psn25lons_v3.dat', 'rb')
        lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
        lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

        xpts, ypts = mplot(lons_mask, lats_mask)

        return region_mask, xpts, ypts
    else:
        return region_mask
    

def getDate(year, month, day):
    """ Get date string from year month and day"""

    return str(year)+'%02d' %month+'%02d' %day

def getNESOSIM(dF, fileSnow, outSnowVar='snow_depth_N', outDensityVar='snow_density_N'):
    """
    Load relevant NESOSIM snow data file and assign to freeboard values

    Args:
        dF (data frame): Pandas dataframe
        fileSnow (string): NESOSIM file path
        outSnowVar (string): Name of snow depth column
        outDensityVar (string): Name of snow density column

    Returns:
        dF (data frame): dataframe updated to include colocated NESOSIM (and dsitributed) snow data
    
    Versions:
        v2: Dropped basemap support and simplified
    """
    
    dateStrStart= getDate(dF.datetime.iloc[0].year, dF.datetime.iloc[0].month, dF.datetime.iloc[0].day)
    dateStrEnd= getDate(dF.datetime.iloc[-1].year, dF.datetime.iloc[-1].month, dF.datetime.iloc[-1].day)
    print('Check dates (should be within a day):', dateStrStart, dateStrEnd)
    
    dN = xr.open_dataset(fileSnow)

    # Get NESOSIM snow depth and density data for the date in the granule
    dNday = dN.sel(day=int(dateStrStart))
    
    # Get NESOSIM data for that data
    lonsN = np.array(dNday.longitude).flatten()
    latsN = np.array(dNday.latitude).flatten()

    # Get dates at start and end of freeboard file
    snowDepthNDay = np.array(dNday.snowDepth).flatten()
    snowDensityNDay = np.array(dNday.density).flatten()
    iceConcNDay = np.array(dNday.iceConc).flatten()
    
    # Remove data where snow depths less than 0 (masked).
    # Might need to chek if I need to apply any other masks here.
    mask=np.where((snowDepthNDay<0.01)|(snowDepthNDay>1)|(iceConcNDay<0.01)|np.isnan(snowDensityNDay))

    snowDepthNDay[mask]=np.nan
    snowDensityNDay[mask]=np.nan
    
    snowDepthNDay=snowDepthNDay
    snowDensityNDay=snowDensityNDay
    
    # I think it's better to declare array now so memory is allocated before the loop?
    snowDepthGISs=np.zeros((dF.shape[0]))
    snowDensityGISs=np.zeros((dF.shape[0]))
    
    # Should change this to an apply or lamda function 
    for x in range(dF.shape[0]):
        
        # Use nearest neighbor to find snow depth at IS2 point
        #snowDepthGISs[x] = griddata((xptsDay, yptsDay), snowDepthDay, (dF['xpts'].iloc[x], dF['ypts'].iloc[x]), method='nearest') 
        #snowDensityGISs[x] = griddata((xptsDay, yptsDay), densityDay, (dF['xpts'].iloc[x], dF['ypts'].iloc[x]), method='nearest')

        # Think this is the much faster way to find nearest neighbor!
        dist = np.sqrt((latsN-dF['lat'].iloc[x])**2+(lonsN-dF['lon'].iloc[x])**2)
        index_min = np.argmin(dist)
        snowDepthGISs[x]=snowDepthNDay[index_min]
        snowDensityGISs[x]=snowDensityNDay[index_min]
   
        
    dF[outSnowVar] = pd.Series(snowDepthGISs, index=dF.index)
    dF[outDensityVar] = pd.Series(snowDensityGISs, index=dF.index)

    return dF
    