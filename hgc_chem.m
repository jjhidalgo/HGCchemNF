close all
clear all
%
% Initializes data structures/
%
[grid, par, por,K, t,saveopt, restart] = hgc_initialize();
% Geometry 
%
grid.A = 1;
grid.angle = -90; %( in degrees);
%
% Discretization
% (cells) 
% if Nx=0 then it is computed so that dx=dz.
%
grid.Nz = 200;
% (cells) 
% if Nx=0 then it is computed so that dx=dz.
%
if (grid.Nx==0)
  grid.Nx = round(grid.A*grid.Nz);
end
%
% Time
%
t.Tmax = 3.0;
t.Tpar = 0.01;
t.inj = 99;
%
% Physical parameters
%
par.IsPeriodic = true();
par.Ra = 500;
par.R = 0.0;
par.G = 1.0;
par.Np = 0.11;
par.denslaw = 1; %(1:PG+water; 2:Salt+water; 0: rho constant);
par.reactlaw = 0; %(0:None; 1:Calcite-cubic);
par.u1 = 10; %end members
par.u2 = 0.49;
par.displ = 0.00;
par.dispt = 0.0;
par.xini = 0.5; %initial position of the front.
par.noise = 0.2; %amplitud of noise in the to intial solution (0.01).
%
% Permeability
%
K.load = false();%
K.file = 'K-0001.mat';
K.isHet = true();
K.isReact = true();
%
% Output options
%
saveopt.conc = 2; %(0:No, 1: Save; 2:Show on screen, 3:show&save). 
saveopt.vel = 0;
saveopt.por = 3;
%
restart.do = false();
restart.Tmax = 10.0;
restart.file = '';
%
[grid, par, por, K, restart] = hgc_default(grid, par, por, K, restart);
%
% Porosity
por.ini = 0.2*ones(grid.Nz,grid.Nx);
%por.char = 0.2;
%por.ini = nthroot(kk,3);
por.por = ones(grid.Nz,grid.Nx);
%load por.mat
%por.ini = bb;
%clear bb
%
% Permeability
K.kperm = ones(grid.Nz,grid.Nx);
%K.kperm  = kk;
%clear kk
%clear kperm

%
% Initial condicitons
%
c = c_ini(grid,par); c=flipud(c);
%save('ini.mat','c')
%load ini.mat
%
% Default values derived from given data.
% 
[telapsed] = gravitycurrent(grid,par,por,K,c,t,saveopt,restart);
%
return