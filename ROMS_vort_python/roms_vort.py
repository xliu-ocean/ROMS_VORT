import numpy as np
import xarray as xr

################### User Interface #########################################################
# grid information
N=15;            # vertical layers
# please note that this package supports only the vertical stretching Vstretching=1 &vtransform = 1   # vertical transform option (1=old(Vstretching=1, Vtransform=1)
                 #                            2=new(Vstretching=2, Vtransform=4) 
                 # please see https://www.myroms.org/wiki/Vertical_S-coordinate for detail)

# path of the data files
data_dir = '../../canyon/'  # direction of DIAG and AVG output
ideal_case = 1   # idealized or real-condition case (1=idealized 0=real)
if ideal_case==1:
    grdfile = data_dir+'roms_avg.nc'
    avgfile = data_dir+'roms_avg.nc'
    diafile = data_dir+'roms_dia.nc'
else:
    grdfile = data_dir+'ROMS_FILES/ocean_grd.nc'
    avgfile = data_dir+'ocean_avg.nc'
    diafile = data_dir+'ocean_dia.nc'
###########################################################################################
# load grid data (xarray)
ds_grid = xr.open_dataset(grdfile)

pm=ds_grid.pm.data
pn=ds_grid.pn.data
h =ds_grid.h.data
pm_p = 0.25*(pm[:-1,:-1]+pm[1:,:-1]+pm[:-1,1:]+pm[1:,1:])
pn_p = 0.25*(pn[:-1,:-1]+pn[1:,:-1]+pn[:-1,1:]+pn[1:,1:])
pm_u = 0.5*(pm[:,:-1]+pm[:,1:])
pn_v = 0.5*(pn[:-1,:]+pn[1:,:])
h_p = 0.25*(h[:-1,:-1]+h[1:,:-1]+h[:-1,1:]+h[1:,1:])
h_u = 0.5*(h[:,:-1]+h[:,1:])
h_v = 0.5*(h[:-1,:]+h[1:,:])
# fcor=ds_grid.pm.data
theta_s = ds_grid.theta_s.data
theta_b = ds_grid.theta_b.data
hc = ds_grid.hc.data
vtransform = ds_grid.Vtransform.data
vstretching = ds_grid.Vstretching.data

# if ideal_case:
#     xgrid = ds_grid.x_rho.data
#     ygrid = ds_grid.y_rho.data
#     xu = ds_grid.x_u.data
#     yu = ds_grid.y_u.data
#     xv = ds_grid.x_v.data
#     yv = ds_grid.y_v.data
#     xp = ds_grid.x_psi.data
#     yp = ds_grid.y_psi.data
# else:    
#     longrid = ds_grid.lon_rho.data
#     latgrid = ds_grid.lat_rho.data
#     lonu = ds_grid.lon_u.data
#     latu = ds_grid.lat_u.data
#     lonv = ds_grid.lon_v.data
#     latv = ds_grid.lat_v.data
#     lonp = ds_grid.lon_psi.data
#     latp = ds_grid.lat_psi.data

# load avg data
ds_avg = xr.open_dataset(avgfile)
ds_dia = xr.open_dataset(diafile)

def zlevs_ROMS_VORT(type,it0):
    global ds_grid, ds_avg
    global h
    global N
    theta_s = ds_grid.theta_s.data
    theta_b = ds_grid.theta_b.data
    hc = ds_grid.Tcline.data
    vtransform = ds_grid.Vtransform.data
    vstretching = ds_grid.Vstretching.data
    h = ds_grid.h.data
    zeta = ds_avg.zeta.isel(ocean_time=slice(it0,it0+1)).data
    (M,L) = h.shape
    sc_r = np.zeros(N);
    Cs_r = np.zeros(N);
    sc_w = np.zeros(N+1);
    Cs_w = np.zeros(N+1);
    def csf(sc,ts,tb):
        if ts > 0:
            csrf=(1-np.cosh(sc*ts))/(np.cosh(ts)-1)
        else:
            csrf=-sc**2
        if tb > 0:
            h = (np.exp(tb*csrf)-1)/(1-np.exp(-tb))
        else:
            h  = csrf
        return h
        
    if vtransform == 2 & vstretching == 4:
        ds=1./N
        if type=='w':
            sc_w[0] = -1.0
            sc_w[N] =  0
            Cs_w[0] = -1.0
            Cs_w[N] =  0
      
            sc_w[1:N] = ds*(np.arange(1,N)-N)
            Cs_w=csf(sc_w, theta_s,theta_b)
            Nw=N+1
        else:
            sc= ds*(np.arange(1,N+1)-N-0.5)
            Cs_r=csf(sc, theta_s,theta_b)
            sc_r=sc
            Nw = N
    else:
        cff1=1./np.sinh(theta_s)
        cff2=0.5/np.tanh(0.5*theta_s)
        if type=='w':
            sc=(np.arange(0,N+1)-N)/N
            Nw=N+1
        else:
            sc=(np.arange(1,N+1)-N-0.5)/N
            Nw = N
        Cs=(1.-theta_b)*cff1*np.sinh(theta_s*sc)+theta_b*(cff2*np.tanh(theta_s*(sc+0.5))-0.5)
    hinv=1./h
    z=np.zeros((Nw,M,L));
    if vtransform == 2 & vstretching == 4:
        if type=='w':
            cff1=Cs_w
            cff2=sc_w+1
            sc=sc_w
        else:
            cff1=Cs_r
            cff2=sc_r+1
            sc=sc_r
        h2=(h+hc)
        cff=hc*sc
        h2inv=1./h2
        for k in np.arange(0,Nw):
            z0=cff[k]+cff1[k]*h
            z[k,:,:]=z0*h/(h2) + zeta*(1.+z0*h2inv)
    else:
        cff1=Cs;
        cff2=sc+1;
        cff=hc*(sc-Cs);
        cff2=sc+1;
        for k in np.arange(0,Nw):
            z0=cff[k]+cff1[k]*h;
            z[k,:,:]=z0+zeta*(1.+z0*hinv);
    return z

def vort_2d_ROMS_VORT(dia_var,it0):
    global ds_dia, ds_grid
    global dvdx_2d, dudy_2d
    global pm_p, pn_p
        
    exec('u2d = ds_dia.ubar_'+dia_var+'.isel(ocean_time=slice('+"%d"%it0+','+"%d"%(it0+1)+')).data',globals())
    exec('v2d = ds_dia.vbar_'+dia_var+'.isel(ocean_time=slice('+"%d"%it0+','+"%d"%(it0+1)+')).data',globals())
    
    dvdx_2d = (v2d[0,:,1:]-v2d[0,:,:-1])*pn_p
    dudy_2d = (u2d[0,1:,:]-u2d[0,:-1,:])*pm_p    
        
    vort_2d=dvdx_2d-dudy_2d
    return vort_2d

def vort_2dtrans_ROMS_VORT(dia_var,it0):
    global ds_dia, ds_grid
    global dvdx_2d, dudy_2d
    global pm_p, pn_p, h_p, h_u, h_v
        
    exec('u2d = ds_dia.ubar_'+dia_var+'.isel(ocean_time=slice('+"%d"%it0+','+"%d"%(it0+1)+')).data',globals())
    exec('v2d = ds_dia.vbar_'+dia_var+'.isel(ocean_time=slice('+"%d"%it0+','+"%d"%(it0+1)+')).data',globals())

    dvdx_t_2d = (v2d[0,:,1:]*h_v[:,1:]-v2d[0,:,:-1]*h_v[:,:-1])*pm_p/h_p
    dudy_t_2d = (u2d[0,1:,:]*h_u[1:,:]-u2d[0,:-1,:]*h_u[:-1,:])*pn_p/h_p
        
    vort_t_2d=dvdx_t_2d-dudy_t_2d
    return vort_t_2d

def vort_3d_ROMS_VORT(dia_var,it0):
    global ds_dia, ds_grid
    global dvdx_2d, dudy_2d
    global vort_2d
    global pm_p,pn_p,pm_u,pn_v
    global N
            
    exec('u3d=(ds_dia.u_'+dia_var+'.isel(ocean_time=slice('+"%d"%it0+','+"%d"%(it0+1)+')).data)',globals())
    exec('v3d = ds_dia.v_'+dia_var+'.isel(ocean_time=slice('+"%d"%it0+','+"%d"%(it0+1)+')).data',globals())

    zr3d = zlevs_ROMS_VORT('r',it0)
    zw3d = zlevs_ROMS_VORT('w',it0)
    (ky1,kx1) = pm_p.shape
    kx = kx1+1 ; ky = ky1+1
    dzdx = (zr3d[:,:,1:]-zr3d[:,:,:-1])*pm_u
    dzdy = (zr3d[:,1:,:]-zr3d[:,:-1,:])*pn_v
    dzdx_p = 0.5*(dzdx[:,1:,:]+dzdx[:,:-1,:])
    dzdy_p = 0.5*(dzdy[:,:,1:]+dzdy[:,:,:-1])
    hz_1 = zw3d[1:,:,:]-zw3d[:-1,:,:]
    hz_1u = 0.5*(hz_1[:,:,1:]+hz_1[:,:,:-1])
    hz_1v = 0.5*(hz_1[:,1:,:]+hz_1[:,:-1,:])
    hz_2 = zr3d[1:,:,:]-zr3d[:-1,:,:]
    hz_2u = 0.5*(hz_2[:,:,1:]+hz_2[:,:,:-1])
    hz_2v = 0.5*(hz_2[:,1:,:]+hz_2[:,:-1,:])
    
    hz_2pu = np.zeros((N+1,ky,kx1))
    hz_2pu[1:-1,:,:]=hz_2u; hz_2pu[0,:,:]=hz_2pu[1,:,:]; hz_2pu[N,:,:] = hz_2pu[N-1,:,:]
    hz_2pv = np.zeros((N+1,ky1,kx))
    hz_2pv[1:-1,:,:]=hz_2v; hz_2pv[0,:,:]=hz_2pv[1,:,:]; hz_2pv[N,:,:] = hz_2pv[N-1,:,:]
    dp1dy01 = (u3d[0,:,1:,:]-u3d[0,:,:-1,:])*pn_p
    dp1dx01 = (v3d[0,:,:,1:]-v3d[0,:,:,:-1])*pm_p

    dp1dsy0 = u3d[0,1:,:,:]-u3d[0,:-1,:,:]
    dp1dsx0 = v3d[0,1:,:,:]-v3d[0,:-1,:,:]
    dp1dzy = np.zeros((N+1,ky,kx1))
    dp1dzx = np.zeros((N+1,ky1,kx))
    dp1dzy[1:-1,:,:] = dp1dsy0/hz_2u
    dp1dzx[1:-1,:,:] = dp1dsx0/hz_2v

# extrapolation of 1st and N+1 layers
    dp1dzy[-1,:,:] = dp1dzy[-2,:,:] + (hz_1u[-1,:,:]+hz_1u[-2,:,:])/   \
                    hz_1u[-2,:,:]*(dp1dzy[-1,:,:]-dp1dzy[-2,:,:])
    dp1dzy[0,:,:] = dp1dzy[2,:,:]+ (hz_1u[0,:,:]+hz_1u[1,:,:])/        \
                    hz_1u[1,:,:]*(dp1dzy[1,:,:]-dp1dzy[2,:,:])
    dp1dzx[-1,:,:] = dp1dzx[-2,:,:] + (hz_1v[-1,:,:]+hz_1v[-2,:,:])/   \
                    hz_1v[-2,:,:]*(dp1dzx[-1,:,:]-dp1dzx[-2,:,:])
    dp1dzx[0,:,:] = dp1dzx[2,:,:] + (hz_1v[0,:,:]+hz_1v[1,:,:])/       \
                    hz_1v[1,:,:]*(dp1dzx[1,:,:]-dp1dzx[2,:,:])

    dp1dzy1 = (hz_2pu[:-1,:,:]*dp1dzy[1:,:,:]+hz_2pu[1:,:,:]*dp1dzy[:-1,:,:])/\
              (hz_2pu[:-1,:,:]+hz_2pu[1:,:,:])
    dp1dzx1 = (hz_2pv[:-1,:,:]*dp1dzx[1:,:,:]+hz_2pv[1:,:,:]*dp1dzx[:-1,:,:])/\
              (hz_2pv[:-1,:,:]+hz_2pv[1:,:,:])
    dp1dzx1_u1 = 0.5*(dp1dzx1[:,:,1:]+dp1dzx1[:,:,:-1])
    dp1dzy1_v1 = 0.5*(dp1dzy1[:,1:,:]+dp1dzy1[:,:-1,:])

    dp1dy02 = dp1dzy1_v1*dzdy_p
    dp1dx02 = dp1dzx1_u1*dzdx_p

    dp1dy0=dp1dy01-dp1dy02
    dp1dx0=dp1dx01-dp1dx02
    
    vort_3d=dp1dx0-dp1dy0
    return vort_3d

# terms of 2D depth-averaged vorticity 
# vort_accel = vort_2d_ROMS_VORT("accel",it0)
# print(vort_2d)

# terms of 2D depth-averaged vorticity 
# vort_trans_accel = vort_2dtrans_ROMS_VORT("accel",it0)
# print(vort_trans_accel)
# vort_3d_accel = vort_3d_ROMS_VORT("accel",it0)


# print(vort_3d_accel)