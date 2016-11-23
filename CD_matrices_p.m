function [Az,Bz,Azz,Bzz]= CD_matrices_p(Nz,dz)
%
% Computed the matrices for first and second derivatives.
% Compact schemes.
%
%
% First derivative (Z direction)
% Differentiation matrices of the compact scheme 6th order
%
  alfa = 1/3;
  a = 14/9;
  b = 1/9;
% c = 0
%
%
%
  R = [1 alfa zeros(1,Nz-2)];
  Az = toeplitz(R,R);
% Boundaries (periodic)
  Az(1,end) = alfa;
  Az(Nz,1) = alfa;
%
  C = [0 a/2 b/4  zeros(1,Nz-5) -b/4 -a/2]/dz;
  Bz = toeplitz(-C,C);
  %
% Boundaries
  Az = sparse(Az);
  Bz = sparse(Bz);
%
% Second derivative
%
  alfa = 2/11;
% beta =0;
  a = 12/11;
  b = 3/11;
% c = 0;
%
  ddz = dz*dz;
%
%
  R = [1 alfa zeros(1,Nz-2)];
  Azz = toeplitz(R,R);
% Boundaries
  Azz(1,end) = alfa;
  Azz(Nz,1) = alfa;
%
  C = [-2*(a+b/4) a b/4  zeros(1,Nz-5) b/4 a]/ddz;
  Bzz = toeplitz(C,C);
  Azz = sparse(Azz);
  Bzz = sparse(Bzz);
%
end