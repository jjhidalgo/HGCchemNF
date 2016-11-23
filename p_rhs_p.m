function [ux_avg,uz_avg]= p_rhs_p(grid,par,conc,Am,T,t)

% Computes right hand side of flow system and solves.
%
%
  Nx = grid.Nx;
  Nz = grid.Nz;
  dx = grid.dx;
  dz = grid.dz;
  Np = grid.Nz*grid.Nx;
%
% Density
%
% If horizontal, only average in z direction is computed.
%
  if (abs(grid.angle)>1e-5)
    [rhoavg_z, rhoavg_x] = rho_average(grid,par,conc);
    rhoavg_x1 = reshape(rhoavg_x(:,1:Nx),Np,1); %rho j-1/2  (x direction)
    rhoavg_x2 = reshape(rhoavg_x(:,2:Nx+1),Np,1); %rho j+1/2
  else
    [rhoavg_z] = rho_average(grid,par,conc);
  end
%
  rhoavg_z1 = reshape(rhoavg_z(1:Nz,:),Np,1); %rho i-1/2  (z direction)
  rhoavg_z2 = reshape(rhoavg_z(2:Nz+1,:),Np,1); %rho i+1/2  
%
% BC
%
  pO = 0;
%
  uO = 0;
%
% Source term
%
  S = zeros(Np,1);
  p = zeros(Np,1);
%
% Bouyancy
%
  S = (rhoavg_z2.*T.Ty2 - rhoavg_z1.*T.Ty1)*grid.cos/dx;
%
  if (abs(grid.angle)>1e-5)   
    S = S + (rhoavg_x2.*T.Tx2 - rhoavg_x1.*T.Tx1)*grid.sin/dz; 
  end
%
% Boundary conditions
%
  S(1:Nz,1) = S(1:Nz,1) + uO/dx; %Neuman BC x=0;
  S(1,1) =1; %to fix one point.
%
% Solve system of equations
%
% spparms('spumoni',2);
  p = Am\S;
  p = reshape(p,Nz,Nx);
  p = full(p);
%
% Compute velocities
%
  ux = zeros(Nz,Nx+1);
  uy = zeros(Nz+1,Nx);
% 
  ux(:,1) = uO;
  ux(:,Nx+1) = 0;
  ux(:,2:Nx) = -T.Tx(:,2:Nx).*(p(:,2:Nx) - p(:,1:Nx-1))/dx;
  if (abs(grid.angle)>1e-5)
    ux(:,1:Nx) = ux(:,1:Nx) -T.Tx(:,1:Nx).*rhoavg_x(:,1:Nx)*grid.sin;
  end
%
  uy(2:Nz,:) =  -T.Ty(2:Nz,:).*([p(2:Nz,:)]-[p(1:Nz-1,:)])/dz ...
                -T.Ty(2:Nz,:).*rhoavg_z(2:Nz,:)*grid.cos;
%
% periodic
   uy(1,:) =  -T.Typ.*([p(1,:)]-[p(Nz,:)])/dz ...
                -T.Typ.*rhoavg_z(1,:)*grid.cos;
   uy(Nz+1,:) =  uy(1,:);              
%
%
% Average velocity in the cell.
%
  ux_avg = zeros(Nz,Nx);
  uz_avg = zeros(Nz,Nx);
  ux_avg = 0.5*( ux(:,1:Nx)+ux(:,2:Nx+1) )/dz;
  uz_avg = 0.5*( uy(1:Nz,:)+uy(2:Nz+1,:) )/dx;
%
end
