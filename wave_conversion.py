import xarray as xr
import numpy as np

qvr_dir = '/gws/nopw/j04/seasia/sferrett/era5_qvr_k2-40_p2-30/'

#set up parameters
y0=6.
g=9.8
beta=2.3e-11
radea=6.371e6
y0real= 2*np.pi*radea*y0/360.0   # convert trapping scale to metres
ce=2*y0real**2*beta
c_on_g=ce/g

def project_wave(year,plevs=[850,200]):
    '''
    Project the ERA5 reanalysis onto parabolic cylinder functions, to extract the Kelvin, WMRG, R1 and R2 wave components.
    
    Parameters:
    year: int, the year to process
    plevs: list ints, the pressure levels to use for the projection in hPa (default is [850,200])
            28 levels available 1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700, 650, 600, 
                550, 500, 450, 400, 350, 300, 250, 225, 200, 175, 150, 125, 100, 70
    
    Returns: 
    u,z,v: three xarray DataArrays, the eastward wind, geopotential height and northward wind components of the four wave types, 
    dimension mode is added where mode=0,1,2,3 corresponds to Kelvin, WMRG, R1 and R2 respectively.
    '''
    
    qvr_nc = xr.open_dataset(f'{qvr_dir}qvr_coefficient_{year}_k2-40_p2-30_28plev_era5.nc').sel(p=plevs)    
    qf_mode = qvr_nc.qn
    vf_mode = qvr_nc.vn
    rf_mode = qvr_nc.rn
    d = xr.open_dataset('Dn.nc')

    #Projection onto parabolic functions, loop over the four wave types.

    uf_wave=[]
    zf_wave=[]
    vf_wave=[]
    waves = ['Kelv','WMRG','R1','R2']

    for w in range(len(waves)):
        if w==0:
            #Kelvin
            uf_wave.append(0.5*qf_mode.isel(mode=0)*d.isel(mode=0).expand_dims(dict(mode=[w])))
            zf_wave.append(0.5*qf_mode.isel(mode=0)*d.isel(mode=0).expand_dims(dict(mode=[w]))*c_on_g)
        elif w==1:
            #WMRG
            uf_wave.append(0.5*qf_mode.isel(mode=1)*d.isel(mode=1).expand_dims(dict(mode=[w])))
            zf_wave.append(0.5*qf_mode.isel(mode=1)*d.isel(mode=1).expand_dims(dict(mode=[w]))*c_on_g)
            vf_wave.append( (vf_mode.isel(mode=0)*d.isel(mode=0)).expand_dims(dict(mode=[w])) )
        elif w==2:
            #R1
            uf_wave.append( 0.5*(qf_mode.isel(mode=2)*d.isel(mode=2)-rf_mode.isel(mode=0)*d.isel(mode=0)).expand_dims(dict(mode=[w])) )
            zf_wave.append( 0.5*(qf_mode.isel(mode=2)*d.isel(mode=2)+rf_mode.isel(mode=0)*d.isel(mode=0)).expand_dims(dict(mode=[w]))*c_on_g )
            vf_wave.append( (vf_mode.isel(mode=1)*d.isel(mode=1)).expand_dims(dict(mode=[w])) )
        elif w==3:
            #R2
            uf_wave.append( 0.5*(qf_mode.isel(mode=3)*d.isel(mode=3)-rf_mode.isel(mode=1)*d.isel(mode=1)).expand_dims(dict(mode=[w])) )
            zf_wave.append( 0.5*(qf_mode.isel(mode=3)*d.isel(mode=3)+rf_mode.isel(mode=1)*d.isel(mode=1)).expand_dims(dict(mode=[w]))*c_on_g )
            vf_wave.append( vf_mode.isel(mode=2)*d.isel(mode=2).expand_dims(dict(mode=[w])) )

    #concatenate the four wave types together along 'mode' dimension, for each variable
    #mode = 0,1,2,3 corresponds to Kelvin, WMRG, R1 and R2 respectively
    uzv = [xr.concat(ncs,'mode') for ncs in [uf_wave,zf_wave,vf_wave]]

    #rename the variables to the correct names for the output files, and return as a list of xarrays
    u,z,v = [nc.rename({'Dn':name}) for nc,name in zip(uzv,['eastward_wind','geopotential_height','northward_wind'])]
    return u,z,v

if __name__ == "__main__":
    out_dir = '/path/to/output/directory/'

    #example of how to run the function for a single year, change as needed
    u,z,v = project_wave(2010)
    u.to_netcdf(f'{out_dir}u_wave_2010.nc')
    z.to_netcdf(f'{out_dir}z_wave_2010.nc')
    v.to_netcdf(f'{out_dir}v_wave_2010.nc')
