function [dvdx_t_2d,dudy_t_2d] = vort_2dtrans_ROMSVOR(u2d,v2d,pm_p,pn_p,hdep_u,hdep_v,hdep_p)
%%% get 1/D*d(vbar*D)/dx and 1/D*d(ubar*D)/dy 
%%% for the 2d transport vorticity equation

u2d_trans = u2d.*hdep_u; v2d_trans = v2d.*hdep_v;

dvdx_t_2d = (v2d_trans(2:end,:,:)-v2d_trans(1:end-1,:,:)).*pm_p./hdep_p;
dudy_t_2d = (u2d_trans(:,2:end,:)-u2d_trans(:,1:end-1,:)).*pn_p./hdep_p;