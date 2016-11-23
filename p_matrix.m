function [Am,T]= p_matrix(grid,par,K,conc)
%
% Computes flow system matrix.
%
% Assumes dz=dx=h.
%
  Nz = grid.Nz;
  Nx = grid.Nx;
  %h = grid.h;
  dx = grid.dx;
  dy = grid.dz; %a small inconsistency to spice it up...
  dx2 = dx*dx;
  dy2 = dy*dy;
%  
  Np = Nz*Nx;    
  mu = ones(size(K.kperm));  
  
  if (abs(par.R)>0)

    mu = (exp(par.R*(par.cmax-conc)))./K.kperm;
    
  elseif(K.isHet)
    
    mu = 1./K.kperm;

  end
%
% Transmisibility matrices.
%
  T.Tx = zeros(Nz,Nx+1);
  T.Tx(:,2:Nx) = (2*dy)./(mu(:,1:Nx-1) + mu(:,2:Nx));
%
  T.Ty = zeros(Nz+1,Nx);
  T.Ty(2:Nz,:) = (2*dx)./(mu(1:Nz-1,:)+mu(2:Nz,:));
%
  T.Tx1 = reshape(T.Tx(:,1:Nx),Np,1);
  T.Tx2 = reshape(T.Tx(:,2:Nx+1),Np,1);
%
  T.Ty1 = reshape(T.Ty(1:Nz,:),Np,1);
  T.Ty2 = reshape(T.Ty(2:Nz+1,:),Np,1);
% 
  Ty11 = T.Ty1;
  Ty22 = T.Ty2;
%
% Assemble system of equations
%
  Dp = zeros(Np,1);
  Dp = T.Tx1/dx2 + T.Ty1/dy2 + T.Tx2/dx2 + T.Ty2/dy2;
%
  Am = spdiags(...
   [-T.Tx2/dx2, ...
    -Ty22/dy2,  Dp, -Ty11/dy2,  ...
    -T.Tx1/dx2], ...
   [-Nz, -1, 0, 1, Nz], Np, Np);
% 
end