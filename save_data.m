function saved = save_data(saveopt,grid,par,por_in,K,t,c,ux,uz,restart)
%
% plot and/or saves data according to defined options.
%

  saved = false();
%
% Auxiliar strings
  c_path = './conc/';
  vel_path = './vel/';
  por_path = './por/';
  strR = strcat(['R' num2str(par.R)]);
  strRa = strcat(['Ra' num2str(par.Ra)]);
  num = addzeros(t.iframe+1);
  timesc = t.timesc;
  strtitle = strcat(['t = ' num2str(t.timesc)]);
%
% Concentration
%
  if (saveopt.conc==1 || saveopt.conc==3)
  
    if (t.iframe+1==1)
      mkdir('./conc')
    end

    cfile = strcat(c_path,strRa,'-',strR,'-',num,'.mat');
    save(cfile,'c','grid','timesc','par');
    

  end %if saveconc

  if (saveopt.conc>=2)
    fig=figure(1);
    set(gcf,'Renderer', 'zbuffer')
    s1=subplot(2,1,1);
    surf(c);
    view(2);
    shading('interp');
    axis('equal','off');
    caxis([0 1]);
    title(strtitle);
    drawnow
  end
% Velocity
%
  if (saveopt.vel==1 || saveopt.vel==3)
 
    if (t.iframe+1==1)
      mkdir('./vel')
    end
  
    velfile = strcat(vel_path,strRa,'-',strR,'-',num,'.mat');
    save(velfile,'ux','uz','grid','timesc','par');

  end
%
  if (saveopt.por==1 || saveopt.por==3)
 
    if (t.iframe+1==1)
      mkdir('./por')
    end
    
    por = por_in.por;
    porfile = strcat(por_path,strRa,'-',strR,'-',num,'.mat');
    save(porfile,'por','grid','timesc','par');

  end
  
  if saveopt.por>=2
    s2=subplot(2,1,2);
    surf(por_in.por-por_in.ini);
    maxdiff = max(max(por_in.por-por_in.ini));
    mindiff = min(min(por_in.por-por_in.ini));
    view(2);shading('interp');axis('equal','off');
    title('porosity');
    caxis([mindiff maxdiff+1e-9]);
    drawnow
  end
%    

  if (saveopt.vel>=2)
    figure(2)
    %mysurf(xx,zz,sqrt(ux.*ux+uz.*uz)/sqrt(max(max(ux.*ux+uz.*uz))));
    quiver(ux(1:7:end,1:7:end),uz(1:7:end,1:7:end),0);
    axis('equal')

    title(strtitle);
    drawnow
  end
%
% restart file
  if saveopt.restart
    % we need to update umax, dt and timesc so that
    % there is no jump in the marching time loop after restart.
    umax = max(max(sqrt(ux.^2 + uz.^2)));
    t.dt = time_update_dt(grid,par,por_in,umax);
    t.timesc = t.timesc + t.dt;
    por = por_in;
    save(restart.file,'por','K','ux','uz','umax','t','c'); 
  end

  saved  = true();

end
