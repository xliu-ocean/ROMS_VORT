function vnew = vinterp(var,z,depth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function  vnew = vinterp(var,z,depth)
%
% This function interpolate a 3D variable on a horizontal level of constant
% depth
%
% On Input:  
% 
%    var     Variable to process (3D matrix).
%    z       Depths (m) of RHO- or W-points (3D matrix).
%    depth   Slice depth (scalar; meters, negative).
% 
% On Output: 
%
%    vnew    Horizontal slice (2D matrix). 
%
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2002-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Modified by X. Liu for rearranging matrices
%  in ROMSVOR package
%  email: xh_liu@sio.org.cn
%  May of 2020
%  The ROMSVOR package is distributed on www.github.com/rayliuxh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Lp,Mp,N]=size(z);
%
% Find the grid position of the nearest vertical levels
%
a=z<depth;
levs=squeeze(sum(a,3));
levs(levs==N)=N-1;
warning off
mask=levs./levs;
warning on
[imat,jmat]=meshgrid((1:Mp),(1:Lp));
pos=(levs-1)*Lp*Mp+Lp*(imat-1)+jmat;
pos2=(levs)*Lp*Mp+Lp*(imat-1)+jmat;
pos(isnan(mask))=1;
%
% Do the interpolation
%
z1=z(pos2);
z2=z(pos);
v1=var(pos2);
v2=var(pos);
vnew=mask.*(((v1-v2)*depth+v2.*z1-v1.*z2)./(z1-z2));
return
