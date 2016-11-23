function [rhoavg_z,rhoavg_x] = rho_average(grid,par,conc)  

 
  switch (par.denslaw)
    
    case (0) %constant density
      
      r = zeros(size(conc));
      %ratio = 0;
      
    case (1) %PG+water
      
      c2 = conc.^2;
      c3 = conc.*c2;
      r = (6.18725*c3 -17.8607*c2 + 8.071574*conc);
      %ratio = 1.0247/0.0086;
      
    case 2 %linear

      r = -3.6*conc;
      %ratio = 1;

    case -2 %linear -1

      r = -1.0*conc;
    
    case 3 %canonical rho = c

      r = conc;
    case 7 %tobias
      
       r =-23.333*conc.*conc + 10*conc;
      
  end %switch
  
  
  xNz = zeros(grid.Nz,1);
  
  rhoavg_x = ([xNz r] + [r xNz])/2;
  
  if par.IsPeriodic
    rhoavg_z = ([r(1,:); r] + [r; r(grid.Nz,:)])/2;
  else
    zNx = zeros(1,grid.Nx);
    rhoavg_z = ([zNx; r] + [r; zNx])/2;
  end
  

end
