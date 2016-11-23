function [telapsed]=gravitycurrent(grid,par,por,K,c,t,saveopt,restart)
%
% grid
%  Nx    --> Cells in x direction
%  Nz    --> Idem z direction
%  Lx
%  Lz
%  A       --> Lx/Lz.
%
% par
%   Pe    --> Peclet number
%   R     --> Mobility ratio R = log(mu2/mu1). Fluid 1 up, fluid 2 down.
%   denslaw --> Density law. 1:PG+water; 2:Salt+water;
%
% K (permeability)
%  isHet
%  var_lnk 
%  corr_lenx 
%  corr_lenz 
%
% t
%   Tmax  --> Simulation time 
%   Tpar  --> Data saved and pictured every Tpar secs
%
% c    --> Initial concentration (then concentration at any time).
%
% saveopts
%   conc -->  Saves concentration data (0:No, 1: Yes; 2:Show on screen, 3:show&save).
%   vel  -->  Saves velocity data (0:No, 1: Yes; 2:Show on screen, 3:show&save)
%
tic;
tstart = tic;
%
% Differentiation matrices for the CDs
%
[Ax,Bx,Axx,Bxx] = CD_matrices(grid.Nx,grid.dx);
%
if par.IsPeriodic
  [Az,Bz,Azz,Bzz] = CD_matrices_p(grid.Nz,grid.dz);
else
  [Az,Bz,Azz,Bzz] = CD_matrices(grid.Nz,grid.dz);
end
%
% Auxiliar strings
strR = num2str(par.R);
strRa = num2str(par.Ra);
%
% Matrices and scalars initialization
%
p = zeros(size(c));
ux = zeros(size(c));
uz = zeros(size(c));
umax = 0;
%
%%%t = time_update(t,(t,grid,par,0,0);
%
post_res = [];
%
if (restart.do)
  
  [time,t,c,ux,uz,umax,por,K] = restart_simulation(restart);

else

  Nframes = floor(t.Tmax/t.Tpar) + 1; %If mod(Tm/Tp)~=0, last frame not pictured  
    
  t.time = zeros(1,Nframes);
  t.time(t.iframe+1) = 0;
  
  saveopt.now = true();
  saveopt.post = true();
  save_data(saveopt,grid,par,por,K,t,c,ux,uz,restart);
  %[post_res] = postprocess(c,t,grid,par,post_res);
  pos_res = 1;

  [t saveopt] = time_update(t,grid,par,por,saveopt,umax);
  
end
%
is_first_time = true();
%
% Iteration
%
while t.timesc <= t.Tmax
     
     c = transport(par,por,t,c,ux,uz,Az,Azz,Bz,Bzz,Ax,Axx,Bx,Bxx);
     
     if (K.isReact)
       por.por = porosity(par,por,c,t,Ax,Bx,Az,Bz);
       K.kperm = permeability(por);
     end
%
% Computes flow matrix when viscosity is not constant
% and solves for pressure and velocity.
% If viscosity is constant, flow matrix is computed just once.
%

    if (abs(par.R)>0 || K.isReact || is_first_time)
        if par.IsPeriodic
          [Am,Trans] = p_matrix_p(grid,par,K,c);
          [ux,uz] = p_rhs_p(grid,par,c,Am,Trans,t);
        else
          [Am,Trans] = p_matrix(grid,par,K,c);
          [ux,uz] = p_rhs(grid,par,c,Am,Trans,t);
        end
        is_first_time = false();
    else
      if par.IsPeriodic
        [ux,uz] = p_rhs_p(grid,par,c,Am,Trans,t);
      else
        [ux,uz] = p_rhs(grid,par,c,Am,Trans,t);
      end
    end

    umax = max(max(sqrt(ux.^2 + uz.^2)));
%
% Saves output
%
    if (saveopt.now)
        
        save_data(saveopt,grid,par,por,K,t,c,ux,uz,restart);
        %[post_res] = postprocess(c,t,grid,par,post_res);
        
        t.time(t.iframe+1) = t.timesc;
        disp(strcat(['time = ' num2str(t.timesc)]));

    end
    
    if (saveopt.post && ~saveopt.now)
     
      % edit this function to customize postprocess.
      % [post_res] = postprocess(c,t,grid,par,post_res);     

    end
%
% Next time step value
%
    [t saveopt] = time_update(t,grid,par,por,saveopt,umax);
    
end % while

  telapsed = toc(tstart); %Computation time 

  cfile = strcat('Ra',strRa,'-R',strR,'.mat');

  time = getfield(t,'time');
  save(cfile,'grid','time','telapsed','t','par','K');
    
  save(restart.file,'ux','uz','umax','por','K','t','c'); 

  %save_postprocess(c,t,grid,par,post_res);
  
  disp('Elapsed time');
  disp(strcat([num2str(telapsed/3600) 'hours']));
  disp('End of simulation.');

end %function
