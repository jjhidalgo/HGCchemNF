function [F]= f_trans(par,por,c,ux,uz,Az,Azz,Bz,Bzz,Ax,Axx,Bx,Bxx)
%
% 
  Ra = par.Ra;

%
  dcx = zeros(size(c));
  dcxx = zeros(size(c));
  dcz = zeros(size(c));
  dczz = zeros(size(c));
%
% TRANSPORT EQUATION
  dcx = (Ax\(Bx*c'))';
  dcxx = (Axx\(Bxx*c'))';
  dcz = Az\(Bz*c);
  dczz = Azz\(Bzz*c);
%
  F = -(ux.*dcx + uz.*dcz)./por.por + (1/Ra)*( dcxx + dczz );

%
  clear dcx dcxx dcz dczz dispx dispxx dispz dispzz
%
end
