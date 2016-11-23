function por_out = porosity(par,por,c,t,Ax,Bx,Az,Bz)

  dcx = zeros(size(c));
  dcz = zeros(size(c));
%
% Gradient of c
%
  gradc2 = ((Ax\(Bx*c'))').^2 + (Az\(Bz*c)).^2;
%
% reaction rate
%
  
  rr = react_rate(par,c);
  
  por_out = por.por.* exp(-t.dt.*rr.*(par.Np/par.Ra).*gradc2);
%
% The maximum porosity allowed is  (1/por.ini)
%
  ix = por_out>1./por.ini;
  por_out(ix) = 1.0./por.ini(ix);
  %ix = por_out>1./por.char;
  %por_out(ix) = 1.0./por.char;

  por_out(por_out<1e-10) = 1e-10;
  
end