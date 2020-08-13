%% user interface of the ROMS Vorticity Diagnostic Package (ROMSVOR) %%%%%%

%%% add path of ROMSVOR
addpath('./ROMSVOR');
%%% specify the folder of AVG and DIA files
data_dir = 'canyon3d/';

%%% idea_case=1 if using the idealise tests of ROMS
idea_case = 1; 

%%% parameters of the s-coordinate used in your ROMS case. 
theta_s=5;
theta_b=0.3; 
hc=90;
% theta_s=6.5;
% theta_b=2.5; 
% hc=400;
N=15;
vtransform=1;
vstretching=1;

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
pm=ncread(grdfile,'pm');
pn=ncread(grdfile,'pn');
hdep=ncread(grdfile,'h');
fcor=ncread(grdfile,'f');
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
t3d    = ncread(avgfile,'temp');
zeta2d = ncread(avgfile,'zeta');
ubar2d = ncread(avgfile,'ubar');
vbar2d = ncread(avgfile,'vbar');
% z information of rho pts
zr3d = zlevs_ROMSVOR(hdep,zeta2d(:,:,it0),theta_s,theta_b,hc,N,'r',vtransform,vstretching);
% z information of w pts
zw3d = zlevs_ROMSVOR(hdep,zeta2d(:,:,it0),theta_s,theta_b,hc,N,'w',vtransform,vstretching);
%
u3d_r=zeros(Nx,Ny,N);
v3d_r=zeros(Nx,Ny,N);
u3d_r(2:end-1,:,:) = 0.5*(u3d(1:end-1,:,:,it0)+u3d(2:end,:,:,it0));
v3d_r(:,2:end-1,:) = 0.5*(v3d(:,1:end-1,:,it0)+v3d(:,2:end,:,it0));
t3d_r = t3d(:,:,:,12);
w3d_r = w3d(:,:,:,12);

ubar2d_r=zeros(Nx,Ny);
vbar2d_r=zeros(Nx,Ny);
ubar2d_r(2:end-1,:) = 0.5*(ubar2d(1:end-1,:,it0)+ubar2d(2:end,:,it0));
vbar2d_r(:,2:end-1) = 0.5*(vbar2d(:,1:end-1,it0)+vbar2d(:,2:end,it0));
%%
h_inp = -50;
u_50r = vinterp_ROMSVOR(u3d_r,zr3d,h_inp);
v_50r = vinterp_ROMSVOR(v3d_r,zr3d,h_inp);
t_50r = vinterp_ROMSVOR(t3d  ,zr3d,h_inp);
w_50r = vinterp_ROMSVOR(w3d  ,zw3d,h_inp);

h_inp = -100;
u_100r = vinterp_ROMSVOR(u3d_r,zr3d,h_inp);
v_100r = vinterp_ROMSVOR(v3d_r,zr3d,h_inp);
t_100r = vinterp_ROMSVOR(t3d_r  ,zr3d,h_inp);
w_100r = vinterp_ROMSVOR(w3d  ,zw3d,h_inp);

h_inp = -200;
u_200r = vinterp_ROMSVOR(u3d_r,zr3d,h_inp);
v_200r = vinterp_ROMSVOR(v3d_r,zr3d,h_inp);
t_200r = vinterp_ROMSVOR(t3d  ,zr3d,h_inp);
w_200r = vinterp_ROMSVOR(w3d  ,zw3d,h_inp);

h_inp = -500;
u_500r = vinterp_ROMSVOR(u3d_r,zr3d,h_inp);
v_500r = vinterp_ROMSVOR(v3d_r,zr3d,h_inp);
t_500r = vinterp_ROMSVOR(t3d  ,zr3d,h_inp);
w_500r = vinterp_ROMSVOR(w3d  ,zw3d,h_inp);
%% plotting figure NEW Fig 1
xlim0 = [122 124];
ylim0 = [25 26.5];
% clim0 = [-2 2]*1e-9;
clim0 = [2.9 3.5];
fontsize0 = 16;
ispace = 1;
iscale = 50;
% figure('colormap',color_tbr)
figure
colormap('jet')
subplot(131)
contourf(xgrid/1000,ygrid/1000,t3d(:,:,end,it0),'linestyle','none');
hold on
quiver(xgrid(1:ispace:end,1:ispace:end)/1000,ygrid(1:ispace:end,1:ispace:end)/1000,...
       u3d_r(1:ispace:end,1:ispace:end,end)*iscale,v3d_r(1:ispace:end,1:ispace:end,end)*iscale,0,'k','linewidth',0.5)
set(gca,'clim',clim0,'fontsize',fontsize0)
axis equal
contour(xgrid/1000,ygrid/1000,hdep,[50,100,200,500,1000],'color',[0.4 0.4 0.4]);
title('0 m','fontsize',fontsize0*2)
subplot(132)
contourf(xgrid/1000,ygrid/1000,t_100r,'linestyle','none');
hold on
quiver(xgrid(1:ispace:end,1:ispace:end)/1000,ygrid(1:ispace:end,1:ispace:end)/1000,...
        u_100r(1:ispace:end,1:ispace:end)*iscale,v_100r(1:ispace:end,1:ispace:end)*iscale,0,'k','linewidth',1)
set(gca,'clim',clim0,'fontsize',fontsize0)
axis equal
contour(xgrid/1000,ygrid/1000,hdep,[50,100,200,500,1000],'color',[0.4 0.4 0.4]);
text(124.52,30.8,'(b)','fontsize',fontsize0*1.5)
clbar1 = colorbar('SouthOutside','fontsize',16);
title('100 m','fontsize',fontsize0*2)
set(clbar1,'position',[0.3625    0.2581    0.3078    0.0222])
set(clbar1,'fontsize',fontsize0+8)
subplot(133)
contourf(xgrid/1000,ygrid/1000,t_200r,'linestyle','none');
hold on
quiver(xgrid(1:ispace:end,1:ispace:end)/1000,ygrid(1:ispace:end,1:ispace:end)/1000,...
        u_200r(1:ispace:end,1:ispace:end)*iscale,v_200r(1:ispace:end,1:ispace:end)*iscale,0,'k','linewidth',1)
set(gca,'clim',clim0,'fontsize',fontsize0)
axis equal
contour(xgrid/1000,ygrid/1000,hdep,[50,100,200,500,1000],'color',[0.4 0.4 0.4]);
title('200 m','fontsize',fontsize0*2)
%% plot w (Fig.3 in the draft)
clim0 = [-1 1]*5e-4;
fontsize0 = 8;
figure;
% colormap('jet')
subplot(131)
hold on
contourf(xgrid,ygrid,w_50r,(-1:0.02:1)*1e-3,'linestyle','none')
set(gca,'clim',clim0,'fontsize',fontsize0)
contour(xgrid,ygrid,hdep,[100,200,500,1000],'k');
% text(121.6,25.05,'(a)','fontsize',fontsize0*1.5)
title('50 m','fontsize',fontsize0*2)
axis equal

subplot(132)
hold on
contourf(xgrid,ygrid,w_100r,(-1:0.02:1)*1e-3,'linestyle','none')
set(gca,'clim',clim0,'fontsize',fontsize0)
contour(xgrid,ygrid,hdep,[100,200,500,1000],'k');
% text(121.6,25.05,'(b)','fontsize',fontsize0*1.5)
title('100 m','fontsize',fontsize0*2)
axis equal
clbar1 = colorbar('SouthOutside','fontsize',16,'position',[0.3799    0.2449    0.2904    0.0237]);
set(clbar1,'fontsize',fontsize0+8)

subplot(133)
hold on
contourf(xgrid,ygrid,w_200r,(-1:0.02:1)*1e-3,'linestyle','none')
set(gca,'clim',clim0,'fontsize',fontsize0)
contour(xgrid,ygrid,hdep,[100,200,500,1000],'k');
text(121.6,25.05,'(c)','fontsize',fontsize0*1.5)
axis equal
title('200 m','fontsize',fontsize0*2)
