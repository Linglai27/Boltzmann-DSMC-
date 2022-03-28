function [Vout, u_dir, sigma] = Collision(V)   
        N = size(V,1);
        V1 = V(1:N/2,:);
        V2 = V(N/2+1:N,:);
        u = V1 - V2;
        u_perp=u(:,1).^2+u(:,2).^2;
        u_abs=sqrt(u_perp+u(:,3).^2);
        u_perp=sqrt(u_perp);
        u_dir=u./u_abs;
        
    % Maxwell molecules collisions (sigma-direction is uniformly distributed over the unit sphere)
        COS=1-2*rand(N/2,1);     % so that sin(theta)d(theta) is uniformly distributed over [0,pi] and the direction (theta,phi) is uniformly distributed over unit sphere
        ONEminusCOS=(1-COS)/2;   % include 1/2 here rather than  in the formulas for D later,    1/2 = m_apha_beta/m_alpha
        SIN=sqrt(1-COS.^2)/2;    % include 1/2 here rather than  in the formulas for D later,    1/2 = m_apha_beta/m_alpha
        
        phi=2*pi*rand(N/2,1);   % azimutal angle
    %abs_SINphi=u_abs.*sin(phi); 
        COSphi=cos(phi); %z_COSphi=u(:,3).*COSphi;
        
    %h=[ -(u(:,1).*z_COSphi - u(:,2).*abs_SINphi)./u_perp    -(u(:,2).*z_COSphi + u(:,1).*abs_SINphi)./u_perp    (u_perp.*COSphi) ];       
    %D=h.*SIN + u.*ONEminusCOS;
    %D=[-(u(:,1).*z_COSphi - u(:,2).*abs_SINphi)./u_perp    -(u(:,2).*z_COSphi + u(:,1).*abs_SINphi)./u_perp    u_perp.*COSphi].*SIN + u.*ONEminusCOS;
          
    %D=[([-u(:,1) -u(:,2)].*(u(:,3).*COSphi)   + [u(:,2) -u(:,1)].*(u_abs.*sin(phi)))./u_perp    u_perp.*COSphi].*SIN + u.*ONEminusCOS;           
        D=[([u(:,1) u(:,2)].*(u(:,3).*COSphi)   + [-u(:,2) u(:,1)].*(u_abs.*sin(phi)))./u_perp    -u_perp.*COSphi].*SIN - u.*ONEminusCOS;           
    
    %D=[u_perp.*COSphi    -(u(:,2).*u(:,1).*COSphi + u(:,3).*abs_SINphi)./u_perp    -(u(:,1).*z_COSphi - u(:,2).*abs_SINphi)./u_perp ].*SIN + u.*ONEminusCOS;

    %%sigma = [SIN.*COSphi SIN.*sin(phi) COS];
    
    % Alternatively 2*D=-du=u-|u|sigma, where sigma is a unit vector given by (theta,phi) angles so:
    % sigma=[COSphi.*SIN*2  sin(phi).*SIN*2  COS];
    % D2=(u-u_abs.*sigma)/2;
    
    %Vout=V+[-D;D]; 
        Vout=V+[D;-D]; 
        sigma = (u+2*D)./u_abs;
    end
       