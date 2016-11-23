function [Am,T]= p_matrix_p(grid,par,K,conc)
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
  T.Ty(2:Nz,:) = (2*dx)./(mu(1:Nz-1,:) + mu(2:Nz,:));
%
  T.Typ = (2*dx)./(mu(Nz,:)+mu(1,:));
  
  T.TypUp = zeros(Np,1);%T.TypUp(end-Nx+1:end,1) = T.Typ(:);
  T.TypDw = zeros(Np,1);%T.TypDw(1:Nx,1) = T.Typ(:);
  %T.Ty(1,:) = (2*dx)./(mu(Nz,:)+mu(1,:));
  %T.Ty(Nz,:) = T.Ty(1,:);
 
  ixup = 1:Nz:Np-Nz+1;
  ixup = ixup + Nz-1;
  ixdw = 1:Nz:Np-Nz+1;
    
  T.TypUp(ixup,1) = T.Typ(:);
  T.TypDw(ixdw,1) = T.Typ(:);
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
%
% Assemble system of equations
%
  Dp = zeros(Np,1);
  Dp = T.Tx1/dx2 + T.Ty1/dy2 + T.Tx2/dx2 + T.Ty2/dy2;
%
% Periodic
%
  Dp(1:Nz:Np,1) = Dp(1:Nz:Np,1) + T.Typ(:)/dy2;
  Dp(Nz:Nz:Np,1) = Dp(Nz:Nz:Np,1) + T.Typ(:)/dy2;
%
  %Am = spdiags(...
  % [-T.Tx2/dx2, ...
  %  -T.TypDw/dy2,...
  %  -Ty22/dy2,  Dp, -Ty11/dy2,  ...
  %  -T.TypUp/dy2, ...
  %  -T.Tx1/dx2], ...
  % [-Nz,-Nz+1, -1, 0, 1,Nz-1, Nz], Np, Np);
 Am = spdiags(...
   [-T.Tx2/dx2, ...
    -T.TypDw/dy2,...
    -Ty22/dy2,  [1; Dp(2:end)], [-Ty11(1)/dy2; 0;-Ty11(3:end)/dy2],  ...
    [-T.TypUp(1:Nz-1)/dy2; 0; -T.TypUp(Nz+1:end)/dy2], ...
    [-T.Tx1(1:Nz)/dx2; 0; -T.Tx1(Nz+2:end)/dx2]], ...
   [-Nz,-Nz+1, -1, 0, 1,Nz-1, Nz], Np, Np);
 
end