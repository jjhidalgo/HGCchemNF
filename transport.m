function c = transport(par,por,t,c,ux,uz,Az,Azz,Bz,Bzz,Ax,Axx,Bx,Bxx);

  for istage = 1:3;
        
    [F] = f_trans(par,por,c,ux,uz,Az,Azz,Bz,Bzz,Ax,Axx,Bx,Bxx);
        
    if istage==1;
        c = c + t.dt*( (8/15)*F );
        F1 = F;

    elseif istage==2;
        c = c + t.dt*( (5/12)*F + (-17/60)*F1 );
        F1 = F;
    else
        c = c + t.dt*( (3/4)*F + (-5/12)*F1 );
    end
%
% Boundary conditions
%
    if (t.timesc<=t.inj)
      c(:,1) = 1;
    else
     c(:,1) = c(:,2);
    end
    
    c(:,end) = c(:,end-1);
    
    if ~par.IsPeriodic
      c(1,:) = c(2,:);
      c(end,:) = c(end-1,:);
    end

  end %for      
  
end