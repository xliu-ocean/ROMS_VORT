function [dp1dx0,dp1dy0] = vort_3d_ROMSVOR(u3d,v3d,pm_p,pn_p,pm_u,pn_v,zr1,zw1,N)

% [kx,ky] = size(pm_1);
% kx1=kx-1;ky1=ky-1;
[kx1,ky1] = size(pm_p);
kx = kx1+1; ky = ky1+1;

dzdx = (zr1(2:end,:,:)-zr1(1:end-1,:,:)).*pm_u;
dzdy = (zr1(:,2:end,:)-zr1(:,1:end-1,:)).*pn_v;
dzdx_p = 0.5*(dzdx(:,1:end-1,:)+dzdx(:,2:end,:));
dzdy_p = 0.5*(dzdy(1:end-1,:,:)+dzdy(2:end,:,:));
hz_1 = zw1(:,:,2:end)-zw1(:,:,1:end-1);
hz_1u = 0.5*(hz_1(1:end-1,:,:)+hz_1(2:end,:,:));
hz_1v = 0.5*(hz_1(:,1:end-1,:)+hz_1(:,2:end,:));
% hz_1p = 0.25*(hz_1(1:end-1,1:end-1,:)+hz_1(1:end-1,2:end,:)+...
%               hz_1(2:end,1:end-1,:)+hz_1(2:end,2:end,:));
hz_2 = zr1(:,:,2:end)-zr1(:,:,1:end-1);
hz_2u = 0.5*(hz_2(1:end-1,:,:)+hz_2(2:end,:,:));
hz_2v = 0.5*(hz_2(:,1:end-1,:)+hz_2(:,2:end,:));

hz_2pu = zeros(kx1,ky,N+1);
hz_2pu(:,:,2:end-1)=hz_2u; hz_2pu(:,:,1) = hz_2pu(:,:,2); hz_2pu(:,:,N+1) = hz_2pu(:,:,N);
hz_2pv = zeros(kx,ky1,N+1);
hz_2pv(:,:,2:end-1)=hz_2v; hz_2pv(:,:,1) = hz_2pv(:,:,2); hz_2pv(:,:,N+1) = hz_2pv(:,:,N);

dp1dy01 = (u3d(:,2:end,:)-u3d(:,1:end-1,:)).*pn_p;
dp1dx01 = (v3d(2:end,:,:)-v3d(1:end-1,:,:)).*pm_p;

dp1dsy0 = u3d(:,:,2:end)-u3d(:,:,1:end-1);
dp1dsx0 = v3d(:,:,2:end)-v3d(:,:,1:end-1);
dp1dzy = zeros(kx1,ky,N+1);
dp1dzx = zeros(kx,ky1,N+1);
dp1dzy(:,:,2:end-1) = dp1dsy0./hz_2u;
dp1dzx(:,:,2:end-1) = dp1dsx0./hz_2v;

%%% extrapolation of 1st and N+1 layers
dp1dzy(:,:,N+1) = dp1dzy(:,:,N-1) + (hz_1u(:,:,N)+hz_1u(:,:,N-1))./...
                hz_1u(:,:,N-1).*(dp1dzy(:,:,N)-dp1dzy(:,:,N-1));
dp1dzy(:,:,1) = dp1dzy(:,:,3)+ (hz_1u(:,:,1)+hz_1u(:,:,2))./...
                hz_1u(:,:,2).*(dp1dzy(:,:,2)-dp1dzy(:,:,3));     
dp1dzx(:,:,N+1) = dp1dzx(:,:,N-1) + (hz_1v(:,:,N)+hz_1v(:,:,N-1))./...
                hz_1v(:,:,N-1).*(dp1dzx(:,:,N)-dp1dzx(:,:,N-1));
dp1dzx(:,:,1) = dp1dzx(:,:,3)+ (hz_1v(:,:,1)+hz_1v(:,:,2))./...
                hz_1v(:,:,2).*(dp1dzx(:,:,2)-dp1dzx(:,:,3));

dp1dzy1 = (hz_2pu(:,:,1:end-1).*dp1dzy(:,:,2:end)+hz_2pu(:,:,2:end).*dp1dzy(:,:,1:end-1))./ ...
        (hz_2pu(:,:,1:end-1)+hz_2pu(:,:,2:end));
dp1dzx1 = (hz_2pv(:,:,1:end-1).*dp1dzx(:,:,2:end)+hz_2pv(:,:,2:end).*dp1dzx(:,:,1:end-1))./ ...
        (hz_2pv(:,:,1:end-1)+hz_2pv(:,:,2:end));
dp1dzx1_u1 = 0.5*(dp1dzx1(2:end,:,:)+dp1dzx1(1:end-1,:,:));
dp1dzy1_v1 = 0.5*(dp1dzy1(:,2:end,:)+dp1dzy1(:,1:end-1,:));
% prsu_1 = uprsgrd(:,:,1)-dp1dzy(:,:,1).*0.5.*hz_1u(:,:,1);
% prsv_1 = vprsgrd(:,:,1)-dp1dzx(:,:,1).*0.5.*hz_1v(:,:,1);
% prsu_1 = dp1dzy1_v1(:,:,1).*hz_1p(:,:,1);
% prsv_1 = dp1dzx1_u1(:,:,1).*hz_1p(:,:,1);

dp1dy02 = dp1dzy1_v1.*dzdy_p;
dp1dx02 = dp1dzx1_u1.*dzdx_p;

dp1dy0=dp1dy01-dp1dy02;
dp1dx0=dp1dx01-dp1dx02;
end