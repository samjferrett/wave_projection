import xarray as xr

#set up parameters
y0=6.
g=9.8
beta=2.3e-11
radea=6.371e6
y0real= 2*np.pi*radea*y0/360.0   # convert trapping scale to metres
ce=2*y0real**2*beta
c_on_g=ce/g


#only calculating for 850hPa
qvr_nc = xr.open_dataset('qvr_coefficient_1980_k2-40_p2-30_28plev_era5.nc').sel(p=850)
qf_mode = qvr_nc.qn
vf_mode = qvr_nc.vn
rf_mode = qvr_nc.rn
d = xr.open_dataset('Dn.nc')

#Projection onto parabolic functions

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

uzv = [xr.concat(ncs,'mode') for ncs in [uf_wave,zf_wave,vf_wave]]
u,z,v = [nc.rename({'Dn':name}) for nc,name in zip(uzv,['eastward_wind','geopotential_height','northward_wind'])]
