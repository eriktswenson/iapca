
# Author: Erik Swenson (latest revision Sep 2023)

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# - read in monthly ERSSTv5 SST and NOAA interpolated OLR downloaded from:
# https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html
# https://downloads.psl.noaa.gov/Datasets/interp_OLR/olr.mon.mean.nc
# - compute Dec-Jan-Feb seasonal mean over tropical Pacific (100E-80W,20S-20N)
#   beginning in 1979 and concatenate together

def read(filename='sst.mnmean.nc',varname='sst',region=[100,280,-20,20],imon=np.array([0,1,2]),date=['1979-12','2022-11']):
# read in data
    ds = xr.open_dataset(filename).sel(time=slice(date[0],date[1]))
    lon = ds['lon'].sel(lon=slice(region[0],region[1])).data; nlon = len(lon)
    dslat = ds['lat']
    revlat = False
# check if latitudes reversed
    if dslat.data[1]<dslat.data[0]:
        revlat = True
        region = [region[0],region[1],region[3],region[2]]

    lat = dslat.sel(lat=slice(region[2],region[3])).data; nlat = len(lat)
    time = np.reshape(ds['time'].data,(-1,12))[:,imon[0]]
    var = ds[varname].sel(lon=slice(region[0],region[1]),lat=slice(region[2],region[3])).data
# seasonal mean
    var = np.mean(np.reshape(var,(-1,12,nlat,nlon))[:,imon,:,:],1)
# variable attributes
    varattrs = {'units': ds[varname].attrs['units'],'var_desc': ds[varname].attrs['var_desc'],'long_name': 'Seasonal mean '+ds[varname].attrs['var_desc'],'dataset': ds[varname].attrs['dataset']}
    ds.close()
# reverse latitudes if necessary
    if revlat:
        lat = np.flip(lat)
        var = np.flip(var,1)

    var = np.reshape(var,(-1,nlat*nlon))
# index for defined grid points
    idef = np.where(np.isnan(var[0,:])==False)[0]
    vardef = var[:,idef]
    n,m = vardef.shape

# subtract climatological mean
    climdef = np.mean(vardef,0)
    vardef = vardef - climdef

# latitude weighting
    g = np.sqrt(np.cos(lat*np.pi/180.0))
    ipole = np.where(abs(lat)==90.0)[0]; npole = len(ipole)
    if npole>0:
        g[ipole] = np.sqrt(1/nlon)

    g = np.reshape(np.repeat(g,nlon),(1,-1))
    gdef = g[:,idef]   

    return(vardef,lon,lat,idef,gdef,climdef,time,varattrs)


def readmerge():

# longitude/latitude for domains
    region = [100,280,-20,20]
    imon = np.array([0,1,2])

# read in SST and OLR
    X1, lon1, lat1, idef1, g1, clim1, time, varattrs1 = read(filename='sst.mnmean.nc',varname='sst',region=region,imon=imon)
    X2, lon2, lat2, idef2, g2, clim2, time, varattrs2 = read(filename='olr.mon.mean.nc',varname='olr',region=region,imon=imon)

# combine SST and OLR
    n = X1.shape[0]; undef = -999
    nlat1 = len(lat1); nlon1 = len(lon1); ndef1 = len(idef1)
    nlat2 = len(lat2); nlon2 = len(lon2); ndef2 = len(idef2)
    Y1 = g1*X1/np.sqrt(n); varX1 = np.sum(Y1**2)
    Y2 = g2*X2/np.sqrt(n); varX2 = np.sum(Y2**2)
    g2[:,:] = g2[:,:]*np.sqrt(varX1/varX2)
    Y2 = g2*X2/np.sqrt(n)
    Y = np.zeros((n,ndef1+ndef2)); g = np.zeros((1,ndef1+ndef2))
    Y[:,0:ndef1] = Y1[:,:]; Y[:,ndef1:(ndef1+ndef2)] = Y2[:,:]
    g[0,0:ndef1] = g1[:,:]; g[0,ndef1:(ndef1+ndef2)] = g2[:,:]
    C = np.zeros(ndef1+ndef2)
    C[0:ndef1] = g1*clim1/np.sqrt(n)
    C[ndef1:(ndef1+ndef2)] = g2*clim2/np.sqrt(n)

    idef = np.arange(ndef1+ndef2)
    idef[0:ndef1] = idef1[:]
    idef[ndef1:(ndef1+ndef2)] = nlat1*nlon1+idef2[:]

    return(Y+C,C,g,idef,ndef1,ndef2,lat1,lon1,lat2,lon2,varattrs1,varattrs2,time)


# - write out single variable to netcdf4 format

def write1(var,fname='null',varname='null',lon=[],lat=[],lev=[],time='2000-01-01',varattrs={},attrs={},levattrs={}):

    import xarray as xr
    import pandas as pd

# var has 3 dimensions (time,lat,lon) OR or 4 dimensions (time,lev,lat,lon)
    if len(var.shape)==3: n,nlat,nlon = var.shape; nlev = 0
    if len(var.shape)==4: n,nlev,nlat,nlon = var.shape

# default values if not provided
    if varname=='null': varname = 'var'
    if fname=='null': fname = varname+'.nc'
    if len(lon)==0: lon = np.arange(nlon)
    if len(lat)==0: lat = np.arange(nlat)
    if len(str(time[0]))==1: time = pd.date_range(start=time,periods=n,freq=pd.DateOffset(years=1))

    if nlev==0:
        dsstr = 'ds = xr.Dataset(data_vars=dict('+varname+'=(["time","lat","lon"],var)),'
        dsstr = dsstr+'coords=dict(lon=(["lon"],lon),lat=(["lat"],lat),time=time))'

    if nlev!=0:
        dsstr = 'ds = xr.Dataset(data_vars=dict('+varname+'=(["time","level","lat","lon"],var)),'
        if len(lev)==0:
            lev = 1000-np.arange(nlev)
        dsstr = dsstr+'coords=dict(lon=(["lon"],lon),lat=(["lat"],lat),level=(["level"],lev),time=time))'

    ldict = {}
    exec(dsstr,locals(),ldict)
    ds = ldict['ds']
    ds[varname].attrs = varattrs
    ds.attrs = attrs
# default coordinate attributes
    ds['lon'].attrs = dict(units='degrees_east',long_name='Longitude',standard_name='longitude',axis='X',coordinate_defines='center')
    ds['lat'].attrs = dict(units='degrees_north',long_name='Latitude',standard_name='latitude',axis='Y',coordinate_defines='center')
    ds['time'].attrs = dict(long_name='Time',standard_name='time',axis='T')
    if nlev!=0:
        if len(levattrs)==0:
            ds['level'].attrs = dict(units='hPa',long_name='Level',standard_name='level',axis='Z',coordinate_defines='point')
        else:
            ds['level'].attrs = levattrs

# write to netcdf4
    ds.to_netcdf(path=fname)


# - plot three patterns together

def plot(var,lon,lat,cint=1,adj=0,colors=['gray'],landcolor='peru',leg=['','',''],savenoplot=False):
   
    n,nlat,nlon = np.shape(var)
    undef = np.zeros((nlat,nlon))
    isnan = np.isnan(var[0,:,:])
    nnan = np.sum(isnan)
    if (nnan>0):
        undef[isnan] = -9999
        landlevels = [-10000,-1000]
        landcolors = [landcolor,landcolor,'white']
        plt.contourf(lon,lat,undef,landlevels,colors=landcolors)

    levels = np.array([-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10])*cint+adj
    lonr = lon[nlon-1]-lon[0]
    latr = lat[nlat-1]-lat[0]
    for k in range(0,n):
        plt.contour(lon,lat,var[k,:,:],levels,colors=colors[k])
        plt.text(lon[0]+lonr*k/n,lat[nlat-1]+latr*0.05,leg[k],color=colors[k],fontsize=13)

    if savenoplot:
        plt.savefig("plot3.pdf")
    else:
        plt.show()


