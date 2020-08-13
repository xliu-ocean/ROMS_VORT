%% user interface of the ROMS Vorticity Diagnostic Package (ROMSVOR) %%%%%%

%%% add path of ROMSVOR
addpath('./ROMSVOR');
%%% specify the folder of AVG and DIA files
data_dir = 'canyon3d/';

%%% idea_case=1 if using the idealise tests of ROMS
idea_case = 1; 

%%% Total number of the s-coordinate layers used in your ROMS case. 
N=15;

%%% indices of the region of interest for residual check
iregion = 1:66;
jregion = 1:49;
% index of the time of interest for residual check
it0 = 12;

% files 
if idea_case
    grdfile = [data_dir,'roms_avg.nc'];
    avgfile = [data_dir,'roms_avg.nc'];
    diafile = [data_dir,'roms_dia.nc'];
else
    grdfile = [data_dir,'ROMS_FILES/ocean_grd_etopo1_v3_v2.1.nc.1'];
    avgfile = [data_dir,'netwv11t1a1b1_avg_2001dd_3_0014.nc'];
    diafile = [data_dir,'netwv11t1a1b1_dia_2001dd_3_0014.nc'];
end
%%%%%%%%%%%%% END of USER INTERFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loading grid 
disp('Loading Geometry information')
if idea_case
    xgrid = ncread(grdfile,'x_rho');
    ygrid = ncread(grdfile,'y_rho');
    xu = ncread(grdfile,'x_u');
    yu = ncread(grdfile,'y_u');
    xv = ncread(grdfile,'x_v');
    yv = ncread(grdfile,'y_v');
    xp = ncread(grdfile,'x_psi');
    yp = ncread(grdfile,'y_psi');
else
    xgrid = ncread(grdfile,'lon_rho');
    ygrid = ncread(grdfile,'lat_rho');
    xu = ncread(grdfile,'lon_u');
    yu = ncread(grdfile,'lat_u');
    xv = ncread(grdfile,'lon_v');
    yv = ncread(grdfile,'lat_v');
    xp = ncread(grdfile,'lon_psi');
    yp = ncread(grdfile,'lat_psi');
end
pm = ncread(grdfile,'pm');
pn = ncread(grdfile,'pn');
hdep = ncread(grdfile,'h');
fcor = ncread(grdfile,'f');
theta_s = ncread(grdfile,'theta_s');
theta_b = ncread(grdfile,'theta_b');
hc = ncread(grdfile,'Tcline');
vtransform = ncread(grdfile,'Vtransform');
vstretching = ncread(grdfile,'Vstretching');
%%% depth, dx and dy information 
[Nx,Ny]=size(hdep);

% m_proj('mercator','lon',[121 125],'lat',[24.5 26.5])
hdep_u = 0.5*(hdep(1:end-1,:)+hdep(2:end,:));
hdep_v = 0.5*(hdep(:,1:end-1)+hdep(:,2:end));
hdep_p = 0.25*(hdep(1:end-1,1:end-1)+hdep(2:end,1:end-1)+...
               hdep(1:end-1,2:end)+hdep(2:end,2:end));
pm_p=0.25*(pm(1:end-1,1:end-1)+pm(2:end,1:end-1)+pm(1:end-1,2:end)+pm(2:end,2:end));
pn_p=0.25*(pn(1:end-1,1:end-1)+pn(2:end,1:end-1)+pn(1:end-1,2:end)+pn(2:end,2:end));
pm_u = 0.5*(pm(1:end-1,:)+pm(2:end,:));
pn_v = 0.5*(pn(:,1:end-1)+pn(:,2:end));
%% loading data 
disp('Loading avg and dia data from NC files')
%%%% avg data
u3d    = ncread(avgfile,'u');
v3d    = ncread(avgfile,'v');
w3d    = ncread(avgfile,'w');
zeta2d = ncread(avgfile,'zeta');
ubar2d = ncread(avgfile,'ubar');
vbar2d = ncread(avgfile,'vbar');
% z information of rho pts
zr3d = zlevs_ROMSVOR(hdep,zeta2d(:,:,it0),theta_s,theta_b,hc,N,'r',vtransform,vstretching);
% z information of w pts
zw3d = zlevs_ROMSVOR(hdep,zeta2d(:,:,it0),theta_s,theta_b,hc,N,'w',vtransform,vstretching);
%%%% dia data
% terms of 2d momentum equations
%
% 2d_acce = 2d_hadv + 2d_fcor + 2d_pgrd + 2d_sstr + 2d_bstr + 2d_hvis
%
u2d_acce = ncread(diafile,'ubar_accel');
v2d_acce = ncread(diafile,'vbar_accel');
u2d_hadv = ncread(diafile,'ubar_hadv');
v2d_hadv = ncread(diafile,'vbar_hadv');
u2d_fcor = ncread(diafile,'ubar_cor');
v2d_fcor = ncread(diafile,'vbar_cor');
u2d_pgrd = ncread(diafile,'ubar_prsgrd');
v2d_pgrd = ncread(diafile,'vbar_prsgrd');
u2d_sstr = ncread(diafile,'ubar_sstr');
v2d_sstr = ncread(diafile,'vbar_sstr');
u2d_bstr = ncread(diafile,'ubar_bstr');
v2d_bstr = ncread(diafile,'vbar_bstr');
u2d_hvis = ncread(diafile,'ubar_hvisc');
v2d_hvis = ncread(diafile,'vbar_hvisc');
% terms of 3d momentum equations
%
% 3d_acce = 3d_hadv + 3d_vadv + 3d_fcor + 3d_pgrd + 3d_hvis + 3d_vvis
%
u3d_acce = ncread(diafile,'u_accel');
v3d_acce = ncread(diafile,'v_accel');
u3d_hadv = ncread(diafile,'u_hadv');
v3d_hadv = ncread(diafile,'v_hadv');
u3d_vadv = ncread(diafile,'u_vadv');
v3d_vadv = ncread(diafile,'v_vadv');
u3d_fcor = ncread(diafile,'u_cor');
v3d_fcor = ncread(diafile,'v_cor');
u3d_pgrd = ncread(diafile,'u_prsgrd');
v3d_pgrd = ncread(diafile,'v_prsgrd');
u3d_hvis = ncread(diafile,'u_hvisc');
v3d_hvis = ncread(diafile,'v_hvisc');
u3d_vvis = ncread(diafile,'u_vvisc');
v3d_vvis = ncread(diafile,'v_vvisc');
%%   
%%%%%%%%%%%%%%%%%%%%% Vorticity calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% there are 3 different forms of the vorticity equations
%%% 1. depth-averaged 2D vorticity equation 
%%% 2. 2D transport vorticity equation
%%% 3. depth-dependent 3D vorticity equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Calculating the 2d vorticity budget
%%% 1.1 get dvdx and dudy
[acce_dvdx_2d,acce_dudy_2d] = vort_2d_ROMSVOR(u2d_acce,v2d_acce,pm_p,pn_p);
[hadv_dvdx_2d,hadv_dudy_2d] = vort_2d_ROMSVOR(u2d_hadv,v2d_hadv,pm_p,pn_p);
[pgrd_dvdx_2d,pgrd_dudy_2d] = vort_2d_ROMSVOR(u2d_pgrd,v2d_pgrd,pm_p,pn_p);
[fcor_dvdx_2d,fcor_dudy_2d] = vort_2d_ROMSVOR(u2d_fcor,v2d_fcor,pm_p,pn_p);
[sstr_dvdx_2d,sstr_dudy_2d] = vort_2d_ROMSVOR(u2d_sstr,v2d_sstr,pm_p,pn_p);
[bstr_dvdx_2d,bstr_dudy_2d] = vort_2d_ROMSVOR(u2d_bstr,v2d_bstr,pm_p,pn_p);
[hvis_dvdx_2d,hvis_dudy_2d] = vort_2d_ROMSVOR(u2d_hvis,v2d_hvis,pm_p,pn_p);

%%% 1.2 get vorticity
acce_curl_2d = acce_dvdx_2d-acce_dudy_2d;
hadv_curl_2d = hadv_dvdx_2d-hadv_dudy_2d;
pgrd_curl_2d = pgrd_dvdx_2d-pgrd_dudy_2d;
fcor_curl_2d = fcor_dvdx_2d-fcor_dudy_2d;
bstr_curl_2d = bstr_dvdx_2d-bstr_dudy_2d;
sstr_curl_2d = sstr_dvdx_2d-sstr_dudy_2d;
hvis_curl_2d = hvis_dvdx_2d-hvis_dudy_2d;

%%% 1.3 terms in the 2D vorticity equation (with signs as follow)
%
% acce_curl_2d + adv_bPV + adv_rPV = JEBAR + curl_Ts + curl_Tb + curl_hvis + Res.
%
adv_bPV = -fcor_curl_2d;
adv_rPV = -hadv_curl_2d;
JEBAR   = pgrd_curl_2d;
curl_Ts = sstr_curl_2d;
curl_Tb = bstr_curl_2d;
curl_hvis = hvis_curl_2d;

%%% 1.4 residual of the 2d equation
Res_2d = acce_curl_2d+adv_bPV+adv_rPV-JEBAR-curl_Ts-curl_Tb-curl_hvis;
disp('Residual of 2d vorticity equations')
disp('amplitude of Residual %% accel  %%% adv_bPV %%% adv_rPV %%% JEBAR   %%% curl_Ts %%% curl_Tb %%% curl_hvis')
disp(['max(|var|)   ',num2str(nanmean(nanmean(abs(Res_2d(iregion,jregion,it0))))),'  ',...
      num2str(nanmean(nanmean(abs(acce_curl_2d(iregion,jregion,it0))))),'  ',...
      num2str(nanmean(nanmean(abs(adv_bPV(iregion,jregion,it0))))),'  ',...
      num2str(nanmean(nanmean(abs(adv_rPV(iregion,jregion,it0))))),'  ',...
      num2str(nanmean(nanmean(abs(JEBAR(iregion,jregion,it0))))),'  ',...
      num2str(nanmean(nanmean(abs(curl_Ts(iregion,jregion,it0))))),'  ',...
      num2str(nanmean(nanmean(abs(curl_Tb(iregion,jregion,it0))))),'  ',...
      num2str(nanmean(nanmean(abs(curl_hvis(iregion,jregion,it0)))))])
%% 2. Calculating the 2d vorticity transport equation 
%%% 2.1 get dvbar*D/dx and dubar*D/dy
[acce_dvdx_2dtrans,acce_dudy_2dtrans] = vort_2dtrans_ROMSVOR(u2d_acce,v2d_acce,...
                                          pm_p,pn_p,hdep_u,hdep_v,hdep_p);
[hadv_dvdx_2dtrans,hadv_dudy_2dtrans] = vort_2dtrans_ROMSVOR(u2d_hadv,v2d_hadv,...
                                          pm_p,pn_p,hdep_u,hdep_v,hdep_p);
[sstr_dvdx_2dtrans,sstr_dudy_2dtrans] = vort_2dtrans_ROMSVOR(u2d_sstr,v2d_sstr,...
                                          pm_p,pn_p,hdep_u,hdep_v,hdep_p);
[bstr_dvdx_2dtrans,bstr_dudy_2dtrans] = vort_2dtrans_ROMSVOR(u2d_bstr,v2d_bstr,...
                                          pm_p,pn_p,hdep_u,hdep_v,hdep_p);
[fcor_dvdx_2dtrans,fcor_dudy_2dtrans] = vort_2dtrans_ROMSVOR(u2d_fcor,v2d_fcor,...
                                          pm_p,pn_p,hdep_u,hdep_v,hdep_p);
[pgrd_dvdx_2dtrans,pgrd_dudy_2dtrans] = vort_2dtrans_ROMSVOR(u2d_pgrd,v2d_pgrd,...
                                          pm_p,pn_p,hdep_u,hdep_v,hdep_p);
[hvis_dvdx_2dtrans,hvis_dudy_2dtrans] = vort_2dtrans_ROMSVOR(u2d_hvis,v2d_hvis,...
                                          pm_p,pn_p,hdep_u,hdep_v,hdep_p);
                                      
%%% 2.2 get terms of the 2D vorticity transport equation
%%% acce_curl_2dtrans + hadv_2dtrans = curlTs_2dtrans + curlTb_2dtrans +
%%% BPT + hvis_curl_2dtrans + Res_2dtrans

acce_curl_2dtrans = acce_dvdx_2dtrans-acce_dudy_2dtrans;
hadv_2dtrans = -(hadv_dvdx_2dtrans-hadv_dudy_2dtrans);
curlTs_2dtrans = sstr_dvdx_2dtrans-sstr_dudy_2dtrans;
curlTb_2dtrans = bstr_dvdx_2dtrans-bstr_dudy_2dtrans;
BPT            = pgrd_dvdx_2dtrans-pgrd_dudy_2dtrans;
hvis_curl_2dtrans = hvis_dvdx_2dtrans-hvis_dudy_2dtrans;
diver_2dtrans  = fcor_dvdx_2dtrans-fcor_dudy_2dtrans; %% should be small due to continuity

%%% 2.3 residual of the 2dtrans equation
Res_2dtrans = acce_curl_2dtrans + hadv_2dtrans - curlTs_2dtrans - curlTb_2dtrans - BPT - hvis_curl_2dtrans - diver_2dtrans;
disp('Residual of 2d transport vorticity equations')
disp('amplitude of Residual %% accel   %% hadv_trans %% curl_Ts %%% curl_Tb %%% BPT     %%% curl_hvis')
disp(['max(|var|)   ',num2str(nanmean(nanmean(abs(Res_2dtrans(iregion,jregion,it0))))),'  ',...
      num2str(nanmean(nanmean(abs(acce_curl_2dtrans(iregion,jregion,it0))))),'  ',...
      num2str(nanmean(nanmean(abs(hadv_2dtrans(iregion,jregion,it0))))),'  ',...
      num2str(nanmean(nanmean(abs(curlTs_2dtrans(iregion,jregion,it0))))),'  ',...
      num2str(nanmean(nanmean(abs(curlTb_2dtrans(iregion,jregion,it0))))),'  ',...
      num2str(nanmean(nanmean(abs(BPT(iregion,jregion,it0))))),'  ',...
      num2str(nanmean(nanmean(abs(hvis_curl_2dtrans(iregion,jregion,it0)))))])
%% 3. Calculating the 3D vorticity equation 
%%% 3.1 get each term of (dv/dx) and (du/dy) 
[acce_dvdx_3d,acce_dudy_3d] = vort_3d_ROMSVOR(u3d_acce(:,:,:,it0),v3d_acce(:,:,:,it0),pm_p,pn_p,pm_u,pn_v,zr3d,zw3d,N);
[pgrd_dvdx_3d,pgrd_dudy_3d] = vort_3d_ROMSVOR(u3d_pgrd(:,:,:,it0),v3d_pgrd(:,:,:,it0),pm_p,pn_p,pm_u,pn_v,zr3d,zw3d,N);
[fcor_dvdx_3d,fcor_dudy_3d] = vort_3d_ROMSVOR(u3d_fcor(:,:,:,it0),v3d_fcor(:,:,:,it0),pm_p,pn_p,pm_u,pn_v,zr3d,zw3d,N);
[hadv_dvdx_3d,hadv_dudy_3d] = vort_3d_ROMSVOR(u3d_hadv(:,:,:,it0),v3d_hadv(:,:,:,it0),pm_p,pn_p,pm_u,pn_v,zr3d,zw3d,N);
[vadv_dvdx_3d,vadv_dudy_3d] = vort_3d_ROMSVOR(u3d_vadv(:,:,:,it0),v3d_vadv(:,:,:,it0),pm_p,pn_p,pm_u,pn_v,zr3d,zw3d,N);
[vvis_dvdx_3d,vvis_dudy_3d] = vort_3d_ROMSVOR(u3d_vvis(:,:,:,it0),v3d_vvis(:,:,:,it0),pm_p,pn_p,pm_u,pn_v,zr3d,zw3d,N);
[hvis_dvdx_3d,hvis_dudy_3d] = vort_3d_ROMSVOR(u3d_hvis(:,:,:,it0),v3d_hvis(:,:,:,it0),pm_p,pn_p,pm_u,pn_v,zr3d,zw3d,N);

%%% 3.2 get each term of 3D vorticity equation
%%%% accel_vort_3d + nonli_vort_3d + diver_vort_3d = vvisc_vort_3d +
%%%% hvisc_vort_3d
accel_vort_3d = acce_dvdx_3d-acce_dudy_3d;
nonli_vort_3d = -((hadv_dvdx_3d+vadv_dvdx_3d)-(hadv_dudy_3d+vadv_dudy_3d));
diver_vort_3d = -(fcor_dvdx_3d-fcor_dudy_3d);
vvisc_vort_3d = vvis_dvdx_3d-vvis_dudy_3d;
hvisc_vort_3d = hvis_dvdx_3d-hvis_dudy_3d;

%%% 3.3 residual of the 3d equation
Res_3d = accel_vort_3d + nonli_vort_3d + diver_vort_3d - vvisc_vort_3d - hvisc_vort_3d;
disp('Residual of 3d vorticity equations at surface layer')
disp('amplitude of Residual %% accel   %% nonlinear %% divergence %%% vviscosity %%% hviscosity')
disp(['max(|var|)   ',num2str(nanmean(nanmean(abs(Res_3d(iregion,jregion,N))))),'  ',num2str(nanmean(nanmean(abs(accel_vort_3d(iregion,jregion,N))))),...
      '  ',num2str(nanmean(nanmean(abs(nonli_vort_3d(iregion,jregion,N))))),'  ',num2str(nanmean(nanmean(abs(diver_vort_3d(iregion,jregion,N))))),'  ',...
      num2str(nanmean(nanmean(abs(vvisc_vort_3d(iregion,jregion,N))))),'  ',num2str(nanmean(nanmean(abs(hvisc_vort_3d(iregion,jregion,N)))))])
%% plotting each term in the vorticity eq
if idea_case == 1
clim0 = [-5 5]*3e-11;
cclines = (-5:0.1:5)*2e-10;
end
fontsize0 = 8;
figure;
colormap('jet')
subplot(231)
contourf(xp,yp,adv_bPV(:,:,it0),cclines,'linestyle','none')
hold on
set(gca,'clim',clim0,'fontsize',fontsize0)
% m_contour(longrid,latgrid,hdep,[200 200],'k','linewidth',1.5)
contour(xgrid,ygrid,hdep,[100,200,500,1000],'k');
% m_usercoast('ecs_high.mat','patch',[0.5 0.5 0.5])
% m_coast('patch',[.5 .5 .5])
% m_grid('box','on','linestyle','none','fontsize',fontsize0*1.5)
% m_text(121.6,25.05,'(a)','fontsize',fontsize0*2)
% set(gca,'xlim',[121.5 123.25],'ylim',[24.85 26.2])
title('adv(bPV)','fontsize',fontsize0*1.5)
subplot(232)
contourf(xp,yp,adv_rPV(:,:,it0),cclines,'linestyle','none')
hold on
set(gca,'clim',clim0,'fontsize',fontsize0)
% m_contour(longrid,latgrid,hdep,[200 200],'k','linewidth',1.5)
contour(xgrid,ygrid,hdep,[100,200,500,1000],'k');
% m_usercoast('ecs_high.mat','patch',[.5 .5 .5])
% m_coast('patch',[.5 .5 .5])
% m_text(121.6,25.05,'(b)','fontsize',fontsize0*2)
% set(gca,'xlim',[121.5 123.25],'ylim',[24.85 26.2])
title('adv(rPV)','fontsize',fontsize0*1.5)
subplot(233)
contourf(xp,yp,JEBAR(:,:,it0),cclines,'linestyle','none')
hold on
set(gca,'clim',clim0,'fontsize',fontsize0)
% m_contour(longrid,latgrid,hdep,[200 200],'k','linewidth',1.5)
contour(xgrid,ygrid,hdep,[100,200,500,1000],'k');
% m_usercoast('ecs_high.mat','patch',[.5 .5 .5])
% m_coast('patch',[.5 .5 .5])
% m_grid('box','on','linestyle','none','fontsize',fontsize0*1.5)
% m_text(121.6,25.05,'(c)','fontsize',fontsize0*2)
% set(gca,'xlim',[121.5 123.25],'ylim',[24.85 26.2])
title('JEBAR','fontsize',fontsize0*1.5)
subplot(234)
contourf(xp,yp,curl_Ts(:,:,it0),cclines,'linestyle','none')
hold on
set(gca,'clim',clim0,'fontsize',fontsize0)
% m_contour(longrid,latgrid,hdep,[200 200],'k','linewidth',1.5)
contour(xgrid,ygrid,hdep,[100,200,500,1000],'k');
% m_usercoast('ecs_high.mat','patch',[.5 .5 .5])
% m_coast('patch',[.5 .5 .5])
% m_grid('box','on','linestyle','none','fontsize',fontsize0*1.5)
% m_text(121.6,25.05,'(d)','fontsize',fontsize0*2)
% set(gca,'xlim',[121.5 123.25],'ylim',[24.85 26.2])
title('curl(\tau_s/H)','fontsize',fontsize0*1.5)
subplot(235)
contourf(xp,yp,curl_Tb(:,:,it0),cclines,'linestyle','none')
hold on
set(gca,'clim',clim0,'fontsize',fontsize0)
contour(xgrid,ygrid,hdep,[100,200,500,1000],'k');
% m_usercoast('ecs_high.mat','patch',[.5 .5 .5])
% m_coast('patch',[.5 .5 .5])
% m_grid('box','on','linestyle','none','fontsize',fontsize0*1.5)
% m_text(121.6,25.05,'(e)','fontsize',fontsize0*2)
% set(gca,'xlim',[121.5 123.25],'ylim',[24.85 26.2])
title('curl(\tau_b/H)','fontsize',fontsize0*1.5)
subplot(236)
contourf(xp,yp,JEBAR(:,:,it0)-adv_bPV(:,:,it0),cclines,'linestyle','none')
hold on
set(gca,'clim',clim0,'fontsize',fontsize0)
contour(xgrid,ygrid,hdep,[100,200,500,1000],'k');
% m_usercoast('ecs_high.mat','patch',[.5 .5 .5])
% m_coast('patch',[.5 .5 .5])
% m_grid('box','on','linestyle','none','fontsize',fontsize0*1.5)
% m_text(121.6,25.05,'(f)','fontsize',fontsize0*2)
% set(gca,'xlim',[121.5 123.25],'ylim',[24.85 26.2])
title('JEBAR-adv(bPV) ','fontsize',fontsize0*1.5)

clbar1 = colorbar('position',[0.925,0.3,0.015,0.4]);
set(gca,'fontsize',fontsize0)
set(clbar1,'fontsize',fontsize0*1.25)
%% integrating within upper 100m 
zw_p1a = 0.25*(zw3d(1:end-1,1:end-1,:)+zw3d(2:end,1:end-1,:)+...
    zw3d(1:end-1,2:end,:)+zw3d(2:end,2:end,:));
hz_1a = zw3d(:,:,2:end)-zw3d(:,:,1:end-1);
hz_1pa = 0.25*(hz_1a(1:end-1,1:end-1,:)+hz_1a(2:end,1:end-1,:)+...
    hz_1a(1:end-1,2:end,:)+hz_1a(2:end,2:end,:));
for ix = 1:Nx-1
    for iy = 1:Ny-1
        k00 = find(zw_p1a(ix,iy,1:end-1)<=-100&zw_p1a(ix,iy,2:end)>-100);
        if ~isempty(k00)
            k100(ix,iy) = k00;
        else
            k100(ix,iy) = 33;
        end
    end
end
for ix = 1:Nx-1
    for iy = 1:Ny-1
        if k100(ix,iy)~=33
           vort_p1_sum100a(ix,iy) = nansum(accel_vort_3d(ix,iy,k100(ix,iy):end).*...
               hz_1pa(ix,iy,k100(ix,iy):end),3)./hdep_p(ix,iy);
           vort_p2_sum100a(ix,iy) = nansum(diver_vort_3d(ix,iy,k100(ix,iy):end).*...
               hz_1pa(ix,iy,k100(ix,iy):end),3)./hdep_p(ix,iy);
           vort_p3_sum100a(ix,iy) = nansum(nonli_vort_3d(ix,iy,k100(ix,iy):end).*...
               hz_1pa(ix,iy,k100(ix,iy):end),3)./hdep_p(ix,iy);
           vort_p4_sum100a(ix,iy) = nansum(vvisc_vort_3d(ix,iy,k100(ix,iy):end).*...
               hz_1pa(ix,iy,k100(ix,iy):end),3)./hdep_p(ix,iy);
           %%%%%%%%%%%%%%%%%%%%%%%%
           if k100(ix,iy)~=1
              vort_p1_sum_bot100a(ix,iy) = nansum(accel_vort_3d(ix,iy,1:k100(ix,iy)-1).*...
                  hz_1pa(ix,iy,1:k100(ix,iy)-1),3)./hdep_p(ix,iy);
              vort_p2_sum_bot100a(ix,iy) = nansum(diver_vort_3d(ix,iy,1:k100(ix,iy)-1).*...
                  hz_1pa(ix,iy,1:k100(ix,iy)-1),3)./hdep_p(ix,iy);
              vort_p3_sum_bot100a(ix,iy) = nansum(nonli_vort_3d(ix,iy,1:k100(ix,iy)-1).*...
                  hz_1pa(ix,iy,1:k100(ix,iy)-1),3)./hdep_p(ix,iy);
              vort_p4_sum_bot100a(ix,iy) = nansum(vvisc_vort_3d(ix,iy,1:k100(ix,iy)-1).*...
                  hz_1pa(ix,iy,1:k100(ix,iy)-1),3)./hdep_p(ix,iy);
           else
              vort_p1_sum_bot100a(ix,iy) = 0;
              vort_p2_sum_bot100a(ix,iy) = 0;
              vort_p3_sum_bot100a(ix,iy) = 0;
              vort_p4_sum_bot100a(ix,iy) = 0;
           end
        else
           vort_p1_sum100a(ix,iy) = nansum(accel_vort_3d(ix,iy,1:end).*...
               hz_1pa(ix,iy,1:end),3)./sum(hz_1pa(ix,iy,1:end),3);
           vort_p2_sum100a(ix,iy) = nansum(diver_vort_3d(ix,iy,1:end).*...
               hz_1pa(ix,iy,1:end),3)./sum(hz_1pa(ix,iy,1:end),3);
           vort_p3_sum100a(ix,iy) = nansum(nonli_vort_3d(ix,iy,1:end).*...
               hz_1pa(ix,iy,1:end),3)./sum(hz_1pa(ix,iy,1:end),3);
           vort_p4_sum100a(ix,iy) = nansum(vvisc_vort_3d(ix,iy,1:end).*...
               hz_1pa(ix,iy,1:end),3)./sum(hz_1pa(ix,iy,1:end),3);
            vort_p1_sum_bot100a(ix,iy) = nan;
            vort_p2_sum_bot100a(ix,iy) = nan;
            vort_p3_sum_bot100a(ix,iy) = nan;
            vort_p4_sum_bot100a(ix,iy) = nan;
        end
    end
end
%%%%
% vort_total_sum100 = -vort_p2_sum100-vort_p3_sum100-vort_p4_sum100+vort_p5_sum100-vort_p6_sum100;
%% Plotting upper 100 m
cclines = (-4:0.04:4)*3e-10;
clim_s=[-1 1]*3e-10;
fontsize0 = 16/2;
ispace = 7;
iscale = 0.4;
% m_proj('mercator','lon',[121.5 124],'lat',[24.5 26.5])
figure
colormap('jet')
subplot(231)
contourf(xp,yp,squeeze(vort_p2_sum100a),cclines,'linestyle','none')
% m_contourf(lonp,latp,squeeze(f_p./hdep_p.*w_sum100)*1e9,cclines,'linestyle','none')
hold on
% m_contour(longrid,latgrid,hdep,[100,200,500,1000],'k');
set(gca,'clim',clim_s,'fontsize',fontsize0)
% m_usercoast('ecs_high.mat','patch',[.5 .5 .5])
title('divergence (0-100 m)','fontsize',fontsize0*1.5)
subplot(232)
contourf(xp,yp,squeeze(vort_p3_sum100a),cclines,'linestyle','none')
hold on
% m_contour(longrid,latgrid,hdep,[100,200,500,1000],'k');
set(gca,'clim',clim_s,'fontsize',fontsize0)
% m_usercoast('ecs_high.mat','patch',[.5 .5 .5])
title('nonlinear (0-100 m)','fontsize',fontsize0*1.5)
subplot(233)
contourf(xp,yp,squeeze(vort_p4_sum100a),cclines,'linestyle','none')
hold on
% m_contour(longrid,latgrid,hdep,[100,200,500,1000],'k');
set(gca,'clim',clim_s,'fontsize',fontsize0)
% m_usercoast('ecs_high.mat','patch',[.5 .5 .5])
title('dissipation (0-100 m)','fontsize',fontsize0*1.5)

subplot(234)
contourf(xp,yp,squeeze(vort_p2_sum_bot100a),cclines,'linestyle','none')
hold on
% m_contour(longrid,latgrid,hdep,[100,200,500,1000],'k');
set(gca,'clim',clim_s,'fontsize',fontsize0)
% m_usercoast('ecs_high.mat','patch',[.5 .5 .5])
title('divergence (100 m-bottom)','fontsize',fontsize0*1.5)
subplot(235)
contourf(xp,yp,squeeze(vort_p3_sum_bot100a),cclines,'linestyle','none')
hold on
% m_contour(longrid,latgrid,hdep,[100,200,500,1000],'k');
set(gca,'clim',clim_s,'fontsize',fontsize0)
% m_usercoast('ecs_high.mat','patch',[.5 .5 .5])
title('nonlinear (100 m-bottom)','fontsize',fontsize0*1.5)
subplot(236)
contourf(xp,yp,squeeze(vort_p4_sum_bot100a),cclines,'linestyle','none')
hold on
% m_contour(longrid,latgrid,hdep,[100,200,500,1000],'k');
set(gca,'clim',clim_s,'fontsize',fontsize0)
% m_usercoast('ecs_high.mat','patch',[.5 .5 .5])
title('dissipation (100 m-bottom)','fontsize',fontsize0*1.5)

% set(gcf,'Position',[466 801 1920/2 984/2])
clbar1 = colorbar('fontsize',fontsize0,'position',[0.9088    0.2988    0.0104    0.4339]);
% set(clbar1,'Label','x10^{-9}','fontsize',16)
% clbar1.Label.String = 'x10^{-9}';
set(clbar1,'fontsize',fontsize0*1.25,'ytick',(-1.2:0.2:1.2)*1e-9)

%% plot transport 2D vorticity equation
cclines = (-5:0.1:5)*3e-10;
clim_s = [-1.5 1.5]*3e-10;
fontsize0 = 8;
ispace = 7;
iscale = 0.4;
fg2 = figure ;
colormap('jet')
subplot(131)
contourf(xp,yp,hadv_2dtrans(:,:,it0),cclines,'linestyle','none')
hold on
contour(xgrid,ygrid,hdep,[100,200,500,1000],'k');
set(gca,'clim',clim_s,'fontsize',fontsize0)
% text(121.6,25.05,'(a)','fontsize',fontsize0*2)
title('advection','fontsize',fontsize0*1.5)
axis equal

subplot(132)
contourf(xp,yp,curlTs_2dtrans(:,:,it0)+curlTb_2dtrans(:,:,it0),cclines,'linestyle','none')
hold on
contour(xgrid,ygrid,hdep,[100,200,500,1000],'k');
set(gca,'clim',clim_s,'fontsize',fontsize0)
% m_usercoast('ecs_high.mat','patch',[.5 .5 .5])
% text(121.6,25.05,'(b)','fontsize',fontsize0*2)
title('stress','fontsize',fontsize0*1.5)
clbar1 = colorbar('SouthOutside','fontsize',fontsize0*1.5,'position',[0.3716    0.1228    0.2904    0.0237]);
% clbar1.Label.String = 'x10^{-9}';
axis equal

subplot(133)
contourf(xp,yp,BPT(:,:,it0),cclines,'linestyle','none')
hold on
contour(xgrid,ygrid,hdep,[100,200,500,1000],'k');
set(gca,'clim',clim_s,'fontsize',fontsize0)
% m_usercoast('ecs_high.mat','patch',[.5 .5 .5])
% text(121.6,25.05,'(c)','fontsize',fontsize0*2)
title('bottom torque','fontsize',fontsize0*1.5)
axis equal
% set(fg2,'position',[ 1          45        960         344]) 