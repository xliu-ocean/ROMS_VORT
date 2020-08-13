function [dvdx_2d,dudy_2d] = vort_2d_ROMSVOR(u2d,v2d,pm_p,pn_p)
%%% get dvbar/dx and dubar/dy for the 2d vorticity equation

dvdx_2d = (v2d(2:end,:,:)-v2d(1:end-1,:,:)).*pm_p;
dudy_2d = (u2d(:,2:end,:)-u2d(:,1:end-1,:)).*pn_p;