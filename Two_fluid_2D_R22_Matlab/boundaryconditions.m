function [FBC,GBC]=boundaryconditions(W,Wfx,Wfy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                            %
%       DEFINE FLUXES                                        %
%                                                            %
%       1. Supersonic outflow                                %
%       2. Supersonic inflow                                 %
%       3. Subsonic outflow                                  %
%       4. Subsonic inflow                                   %
%       5. Solid wall                                        %
%                                                            %
%       l=left                                               %
%       r=right                                              %
%       a=above                                              %
%       b=below                                              %
%                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global NX NY
global gamma1 gamma2

%      l r a b
%BC = [1,1,1,1]
BC_L = 0;     % Left boundary
BC_R = 0;     % right boundary
BC_A = 5;     % above boundary
BC_B = 5;     % below boundary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      NO BOUNDARY                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Left Boundary
if BC_L==0
    FBC(:,1,1) = W(:,1,3).*W(:,1,4);
    FBC(:,1,2) = W(:,1,3).*W(:,1,4).^2+W(:,1,6);
    FBC(:,1,3) = W(:,1,3).*W(:,1,4).*W(:,1,5);
    FBC(:,1,4) = (W(:,1,8)/(gamma1-1)+(1-W(:,1,8))/(gamma2-1)+1).*W(:,1,4).*W(:,1,6)+1/2*W(:,1,3).*W(:,1,4).*(W(:,1,4).^2+W(:,1,5).^2);
    FBC(:,1,5) = W(:,1,7).*W(:,1,3).*W(:,1,4);
    FBC(:,1,6) = gamma1/(gamma1-1)*W(:,1,4).*W(:,1,6).*W(:,1,8)+1/2*W(:,1,3).*W(:,1,4).*(W(:,1,4).^2+W(:,1,5).^2).*W(:,1,8);
end

%Right Boundary
if BC_R==0
    FBC(:,NX+1,1) = W(:,NX,3).*W(:,NX,4);
    FBC(:,NX+1,2) = W(:,NX,3).*W(:,NX,4).^2+W(:,NX,6);
    FBC(:,NX+1,3) = W(:,NX,3).*W(:,NX,4).*W(:,NX,5);
    FBC(:,NX+1,4) = (W(:,NX,8)/(gamma1-1)+(1-W(:,NX,8))/(gamma2-1)+1).*W(:,NX,4).*W(:,NX,6)+1/2*W(:,NX,3).*W(:,NX,4).*(W(:,NX,4).^2+W(:,NX,5).^2);
    FBC(:,NX+1,5) = W(:,NX,7).*W(:,NX,3).*W(:,NX,4);
    FBC(:,NX+1,6) = gamma1/(gamma1-1).*W(:,NX,4).*W(:,NX,6).*W(:,NX,8)+1/2*W(:,NX,3).*W(:,NX,4).*(W(:,NX,4).^2+W(:,NX,5).^2).*W(:,NX,8);
end

%Upper Boundary
if BC_A==0 
    GBC(1,:,1) = W(1,:,3).*W(1,:,5);
    GBC(1,:,2) = W(1,:,3).*W(1,:,5).*W(1,:,4);
    GBC(1,:,3) = W(1,:,3).*W(1,:,5).^2+W(1,:,6);
    GBC(1,:,4) = (W(1,:,8)/(gamma1-1)+(1-W(1,:,8))/(gamma2-1)+1).*W(1,:,5).*W(1,:,6)+1/2*W(1,:,3).*W(1,:,5).*(W(1,:,4).^2+W(1,:,5).^2);
    GBC(1,:,5) = W(1,:,7).*W(1,:,3).*W(1,:,5);
    GBC(1,:,6) = gamma1/(gamma1-1)*W(1,:,5).*W(1,:,6).*W(1,:,8)+1/2*W(1,:,3).*W(1,:,5).*(W(1,:,4).^2+W(1,:,5).^2).*W(1,:,8);
end

%Lower Boundary
if BC_B==0
    GBC(NY+1,:,1) = W(NY,:,3).*W(NY,:,5);
    GBC(NY+1,:,2) = W(NY,:,3).*W(NY,:,5).*W(NY,:,4);
    GBC(NY+1,:,3) = W(NY,:,3).*W(NY,:,5).^2+W(NY,:,6);
    GBC(NY+1,:,4) = (W(NY,:,8)/(gamma1-1)+(1-W(NY,:,8))/(gamma2-1)+1).*W(NY,:,5).*W(NY,:,6)+1/2*W(NY,:,3).*W(NY,:,5).*(W(NY,:,4).^2+W(NY,:,5).^2);
    GBC(NY+1,:,5) = W(NY,:,7).*W(NY,:,3).*W(NY,:,5);
    GBC(NY+1,:,6) = gamma1/(gamma1-1)*W(NY,:,5).*W(NY,:,6).*W(NY,:,8)+1/2*W(NY,:,3).*W(NY,:,5).*(W(NY,:,4).^2+W(NY,:,5).^2).*W(NY,:,8);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Supersonic outflow                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Left Boundary
if BC_L==1
    FBC(:,1,1) = Wfx(:,1,3).*Wfx(:,1,4);
    FBC(:,1,2) = Wfx(:,1,3).*Wfx(:,1,4).^2+Wfx(:,1,6);
    FBC(:,1,3) = Wfx(:,1,3).*Wfx(:,1,4).*Wfx(:,1,5);
    FBC(:,1,4) = (Wfx(:,1,8)/(gamma1-1)+(1-Wfx(:,1,8))/(gamma2-1)+1).*Wfx(:,1,4).*Wfx(:,1,6)+1/2*Wfx(:,1,3).*Wfx(:,1,4).*(Wfx(:,1,4).^2+Wfx(:,1,5).^2);
    FBC(:,1,5) = Wfx(:,1,7).*Wfx(:,1,3).*Wfx(:,1,4);
    FBC(:,1,6) = gamma1/(gamma1-1)*Wfx(:,1,4).*Wfx(:,1,6).*Wfx(:,1,8)+1/2*Wfx(:,1,3).*Wfx(:,1,4).*(Wfx(:,1,4).^2+Wfx(:,1,5).^2).*Wfx(:,1,8);
end

%Right Boundary
if BC_R==1
    FBC(:,NX+1,1) = Wfx(:,2*NX,3).*Wfx(:,2*NX,4);
    FBC(:,NX+1,2) = Wfx(:,2*NX,3).*Wfx(:,2*NX,4).^2+Wfx(:,2*NX,6);
    FBC(:,NX+1,3) = Wfx(:,2*NX,3).*Wfx(:,2*NX,4).*Wfx(:,2*NX,5);
    FBC(:,NX+1,4) = (Wfx(:,2*NX,8)/(gamma1-1)+(1-Wfx(:,2*NX,8))/(gamma2-1)+1)*Wfx(:,2*NX,4).*Wfx(:,2*NX,6)+1/2*Wfx(:,2*NX,3).*Wfx(:,2*NX,4).*(Wfx(:,2*NX,4).^2+Wfx(:,2*NX,5).^2);
    FBC(:,NX+1,5) = Wfx(:,2*NX,7).*Wfx(:,2*NX,3).*Wfx(:,2*NX,4);
    FBC(:,NX+1,6) = gamma1/(gamma1-1)*Wfx(:,2*NX,4).*Wfx(:,2*NX,6).*Wfx(:,2*NX,8)+1/2*Wfx(:,2*NX,3).*Wfx(:,2*NX,4).*(Wfx(:,2*NX,4).^2+Wfx(:,2*NX,5).^2).*Wfx(:,2*NX,8);
end

%Upper Boundary
if BC_A==1
    GBC(1,:,1) = Wfy(1,:,3).*Wfy(1,:,5);
    GBC(1,:,2) = Wfy(1,:,3).*Wfy(1,:,5).*Wfy(1,:,4);
    GBC(1,:,3) = Wfy(1,:,3).*Wfy(1,:,5).^2+Wfy(1,:,6);
    GBC(1,:,4) = (Wfy(1,:,8)/(gamma1-1)+(1-Wfy(1,:,8))/(gamma2-1)+1).*Wfy(1,:,5).*Wfy(1,:,6)+1/2*Wfy(1,:,3).*Wfy(1,:,5).*(Wfy(1,:,4).^2+Wfy(1,:,5).^2);
    GBC(1,:,5) = Wfy(1,:,7)*Wfy(1,:,3)*Wfy(1,:,5);
    GBC(1,:,6) = gamma1/(gamma1-1)*Wfy(1,:,5).*Wfy(1,:,6).*Wfy(1,:,8)+1/2*Wfy(1,:,3).*Wfy(1,:,5).*(Wfy(1,:,4).^2+Wfy(1,:,5).^2).*Wfy(1,:,8);
end

%Lower Boundary
if BC_B==1
    GBC(NY+1,:,1) = Wfy(2*NY,:,3).*Wfy(2*NY,:,5);
    GBC(NY+1,:,2) = Wfy(2*NY,:,3).*Wfy(2*NY,:,5).*Wfy(2*NY,:,4);
    GBC(NY+1,:,3) = Wfy(2*NY,:,3).*Wfy(2*NY,:,5).^2+Wfy(2*NY,:,6);
    GBC(NY+1,:,4) = (Wfy(2*NY,:,8)/(gamma1-1)+(1-Wfy(2*NY,:,8))/(gamma2-1)+1).*Wfy(2*NY,:,5).*Wfy(2*NY,:,6)+1/2*Wfy(2*NY,:,3).*Wfy(2*NY,:,5).*(Wfy(2*NY,:,4).^2+Wfy(2*NY,:,5).^2);
    GBC(NY+1,:,5) = Wfy(2*NY,:,7).*Wfy(2*NY,:,3).*Wfy(2*NY,:,5);
    GBC(NY+1,:,6) = gamma1/(gamma1-1)*Wfy(2*NY,:,5).*Wfy(2*NY,:,6).*Wfy(2*NY,:,8)+1/2*Wfy(2*NY,:,3).*Wfy(2*NY,:,5).*(Wfy(2*NY,:,4).^2+Wfy(2*NY,:,5).^2).*Wfy(2*NY,:,8);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Supersonic inflow                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho   = 1.92691;
u     = -0.33361;
v     = 0e0;
p     = 1.56980e0;
beta  = 1;
alpha = 1;

%Left Boundary
if BC_L==2
    FBC(:,1,1) = rho*u;
    FBC(:,1,2) = rho*u^2+p;
    FBC(:,1,3) = rho*u*v;
    FBC(:,1,4) = (alpha/(gamma1-1)+(1-alpha)/(gamma2-1)+1)*u*p+1/2*rho*u*(u^2+v^2);
    FBC(:,1,5) = beta*rho*u;
    FBC(:,1,6) = gamma1/(gamma1-1)*u*p*alpha+1/2*rho*u*(u^2+v^2)*beta;
end

%Right Boundary
if BC_R==2
    FBC(:,NX+1,1) = rho*u;
    FBC(:,NX+1,2) = rho*u^2+p;
    FBC(:,NX+1,3) = rho*u*v;
    FBC(:,NX+1,4) = (alpha/(gamma1-1)+(1-alpha)/(gamma2-1)+1)*u*p+1/2*rho*u*(u^2+v^2);
    FBC(:,NX+1,5) = beta*rho*u;
    FBC(:,NX+1,6) = gamma1/(gamma1-1)*u*p*alpha+1/2*rho*u*(u^2+v^2)*beta;
end

%Upper Boundary
if BC_A==2
    GBC(1,:,1)  = rho*v;
    GBC(1,:,2)  = rho*v*u;
    GBC(1,:,3)  = rho*v^2+p;
    GBC(1,:,4)  = (alpha/(gamma1-1)+(1-alpha)/(gamma2-1)+1)*v*p+1/2*rho*v*(u^2+v^2);
    GBC(1,:,5)  = beta*rho*v;
    GBC(1,:,6)  = gamma1/(gamma1-1)*v*p*alpha+1/2*rho*v*(u^2+v^2)*beta;
end

%Lower Boundary
if BC_B==2
    GBC(NY+1,:,1)  = rho*v;
    GBC(NY+1,:,2)  = rho*v*u;
    GBC(NY+1,:,3)  = rho*v^2+p;
    GBC(NY+1,:,4)  = (alpha/(gamma1-1)+(1-alpha)/(gamma2-1)+1)*v*p+1/2*rho*v*(u^2+v^2);
    GBC(NY+1,:,5)  = beta*rho*v;
    GBC(NY+1,:,6)  = gamma1/(gamma1-1)*v*p*alpha+1/2*rho*v*(u^2+v^2)*beta;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Subsonic outflow                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pflux = 1;

%Left Boundary
if BC_L==3
    for i=1:NY
        rho   = Wfx(i,1,3);
        u     = Wfx(i,1,4);
        v     = Wfx(i,1,5);
        p     = Wfx(i,1,6);
        beta  = Wfx(i,1,7);
        alpha = Wfx(i,1,8);

        NRK4 = 10;

        dp   = (pflux - p)/NRK4;

        for k=1:NRK4
            [du_k1,drho_k1,dalpha_k1] = RK4_p(rho          ,p     ,alpha            ,dp,+1);
            [du_k2,drho_k2,dalpha_k2] = RK4_p(rho+drho_k1/2,p+dp/2,alpha+dalpha_k1/2,dp,+1);
            [du_k3,drho_k3,dalpha_k3] = RK4_p(rho+drho_k2/2,p+dp/2,alpha+dalpha_k2/2,dp,+1);
            [du_k4,drho_k4,dalpha_k4] = RK4_p(rho+drho_k3  ,p+dp  ,alpha+dalpha_k3  ,dp,+1);
            u     = u     + (du_k1     + 2*du_k2     + 2*du_k3     + du_k4    )/6;
            rho   = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
            alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
        end

        p = pflux;

        FBC(i,1,1) = rho*u;
        FBC(i,1,2) = rho*u^2+p;
        FBC(i,1,3) = rho*u*v;
        FBC(i,1,4) = (alpha/(gamma1-1)+(1-alpha)/(gamma2-1)+1)*u*p+1/2*rho*u*(u^2+v^2);
        FBC(i,1,5) = beta*rho*u;
        FBC(i,1,6) = gamma1/(gamma1-1)*u*p*alpha+1/2*rho*u*(u^2+v^2)*beta;
    end
end
   
%Right Boundary
if BC_R==3
    for i=1:NY
        rho   = Wfx(i,2*NX,3);
        u     = Wfx(i,2*NX,4);
        v     = Wfx(i,2*NX,5);
        p     = Wfx(i,2*NX,6);
        beta  = Wfx(i,2*NX,7);
        alpha = Wfx(i,2*NX,8);

        NRK4 = 10;

        dp   = (pflux - p)/NRK4;

        for k=1:NRK4
            [du_k1,drho_k1,dalpha_k1] = RK4_p(rho          ,p     ,alpha            ,dp,-1);
            [du_k2,drho_k2,dalpha_k2] = RK4_p(rho+drho_k1/2,p+dp/2,alpha+dalpha_k1/2,dp,-1);
            [du_k3,drho_k3,dalpha_k3] = RK4_p(rho+drho_k2/2,p+dp/2,alpha+dalpha_k2/2,dp,-1);
            [du_k4,drho_k4,dalpha_k4] = RK4_p(rho+drho_k3  ,p+dp  ,alpha+dalpha_k3  ,dp,-1);
            u     = u     + (du_k1     + 2*du_k2     + 2*du_k3     + du_k4    )/6;
            rho   = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
            alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
        end

        p = pflux;

        FBC(i,NX+1,1) = rho*u;
        FBC(i,NX+1,2) = rho*u^2+p;
        FBC(i,NX+1,3) = rho*u*v;
        FBC(i,NX+1,4) = (alpha/(gamma1-1)+(1-alpha)/(gamma2-1)+1)*u*p+1/2*rho*u*(u^2+v^2);
        FBC(i,NX+1,5) = beta*rho*u;
        FBC(i,NX+1,6) = gamma1/(gamma1-1)*u*p*alpha+1/2*rho*u*(u^2+v^2)*beta;
    end
end
          
%Upper Boundary
if BC_A==3
    for j=1:NX
    rho   = Wfy(1,j,3);
    u     = Wfy(1,j,4);
    v     = Wfy(1,j,5);
    p     = Wfy(1,j,6);
    beta  = Wfy(1,j,7);
    alpha = Wfy(1,j,8);

    NRK4 = 10;

    dp   = (pflux - p)/NRK4;

    for k=1:NRK4
            [dv_k1,drho_k1,dalpha_k1] = RK4_p(rho          ,p     ,alpha            ,dp,+1);
            [dv_k2,drho_k2,dalpha_k2] = RK4_p(rho+drho_k1/2,p+dp/2,alpha+dalpha_k1/2,dp,+1);
            [dv_k3,drho_k3,dalpha_k3] = RK4_p(rho+drho_k2/2,p+dp/2,alpha+dalpha_k2/2,dp,+1);
            [dv_k4,drho_k4,dalpha_k4] = RK4_p(rho+drho_k3  ,p+dp  ,alpha+dalpha_k3  ,dp,+1);
            v     = v     + (dv_k1     + 2*dv_k2     + 2*dv_k3     + dv_k4    )/6;
            rho   = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
            alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
    end

    p = pflux;

    GBC(1,j,1) = rho*v;
    GBC(1,j,2) = rho*v*u;
    GBC(1,j,3) = rho*v^2+p;
    GBC(1,j,4) = (alpha/(gamma1-1)+(1-alpha)/(gamma2-1)+1)*v*p+1/2*rho*v*(u^2+v^2);
    GBC(1,j,5) = beta*rho*v;
    GBC(1,j,6) = gamma1/(gamma1-1)*v*p*alpha+1/2*rho*v*(u^2+v^2)*beta;
    end
end

%Lower Boundary
if BC_B==3
    for j=1:NX
    rho   = Wfy(2*NY,j,3);
    u     = Wfy(2*NY,j,4);
    v     = Wfy(2*NY,j,5);
    p     = Wfy(2*NY,j,6);
    beta  = Wfy(2*NY,j,7);
    alpha = Wfy(2*NY,j,8);

    NRK4 = 10;

    dp   = (pflux - p)/NRK4;

    for k=1:NRK4
            [dv_k1,drho_k1,dalpha_k1] = RK4_p(rho          ,p     ,alpha            ,dp,-1);
            [dv_k2,drho_k2,dalpha_k2] = RK4_p(rho+drho_k1/2,p+dp/2,alpha+dalpha_k1/2,dp,-1);
            [dv_k3,drho_k3,dalpha_k3] = RK4_p(rho+drho_k2/2,p+dp/2,alpha+dalpha_k2/2,dp,-1);
            [dv_k4,drho_k4,dalpha_k4] = RK4_p(rho+drho_k3  ,p+dp  ,alpha+dalpha_k3  ,dp,-1);
            v     = v     + (dv_k1     + 2*dv_k2     + 2*dv_k3     + dv_k4    )/6;
            rho   = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
            alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
    end

    p = pflux;

    GBC(NY+1,j,1) = rho*v;
    GBC(NY+1,j,2) = rho*v*u;
    GBC(NY+1,j,3) = rho*v^2+p;
    GBC(NY+1,j,4) = (alpha/(gamma1-1)+(1-alpha)/(gamma2-1)+1)*v*p+1/2*rho*v*(u^2+v^2);
    GBC(NY+1,j,5) = beta*rho*v;
    GBC(NY+1,j,6) = gamma1/(gamma1-1)*v*p*alpha+1/2*rho*v*(u^2+v^2)*beta;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Subsonic inflow                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhoflux   = 1.92691;
uflux     = -0.33361;
vflux     = 0e0;
betaflux  = 1;
alphaflux = 1;

%Left Boundary
if BC_L==4
    for i=1:NY
        rho   = Wfx(i,1,3);
        u     = Wfx(i,1,4);
        p     = Wfx(i,1,6);
        alpha = Wfx(i,1,8);

        NRK4 = 10;

        du   = (uflux - u)/NRK4;

        for k=1:NRK4
        [dp_k1,drho_k1,dalpha_k1] = RK4_u(rho          ,p        ,alpha            ,du,+1);
        [dp_k2,drho_k2,dalpha_k2] = RK4_u(rho+drho_k1/2,p+dp_k1/2,alpha+dalpha_k1/2,du,+1);
        [dp_k3,drho_k3,dalpha_k3] = RK4_u(rho+drho_k2/2,p+dp_k2/2,alpha+dalpha_k2/2,du,+1);
        [dp_k4,drho_k4,dalpha_k4] = RK4_u(rho+drho_k3  ,p+dp_k3  ,alpha+dalpha_k3  ,du,+1);
        p     = p     + (dp_k1     + 2*dp_k2     + 2*dp_k3     + dp_k4    )/6;
        rho   = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
        alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
        end

        rho   = rhoflux;
        u     = uflux;
        v     = vflux;
        beta  = betaflux;
        alpha = alphaflux;

        FBC(i,1,1) = rho*u;
        FBC(i,1,2) = rho*u^2+p;
        FBC(i,1,3) = rho*u*v;
        FBC(i,1,4) = (alpha/(gamma1-1)+(1-alpha)/(gamma2-1)+1)*u*p+1/2*rho*u*(u^2+v^2);
        FBC(i,1,5) = beta*rho*u;
        FBC(i,1,6) = gamma1/(gamma1-1)*u*p*alpha+1/2*rho*u*(u^2+v^2)*beta;
    end
end

%Right Boundary
if BC_R==4
    for i=1:NY
        rho   = Wfx(i,2*NX,3);
        u     = Wfx(i,2*NX,4);
        p     = Wfx(i,2*NX,6);
        alpha = Wfx(i,2*NX,8);

        NRK4 = 10;

        du   = (uflux - u)/NRK4;

        for k=1:NRK4
        [dp_k1,drho_k1,dalpha_k1] = RK4_u(rho          ,p        ,alpha            ,du,-1);
        [dp_k2,drho_k2,dalpha_k2] = RK4_u(rho+drho_k1/2,p+dp_k1/2,alpha+dalpha_k1/2,du,-1);
        [dp_k3,drho_k3,dalpha_k3] = RK4_u(rho+drho_k2/2,p+dp_k2/2,alpha+dalpha_k2/2,du,-1);
        [dp_k4,drho_k4,dalpha_k4] = RK4_u(rho+drho_k3  ,p+dp_k3  ,alpha+dalpha_k3  ,du,-1);
        p     = p     + (dp_k1     + 2*dp_k2     + 2*dp_k3     + dp_k4    )/6;
        rho   = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
        alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
        end

        rho   = rhoflux;
        u     = uflux;
        v     = vflux;
        beta  = betaflux;
        alpha = alphaflux;

        FBC(i,NX+1,1) = rho*u;
        FBC(i,NX+1,2) = rho*u^2+p;
        FBC(i,NX+1,3) = rho*u*v;
        FBC(i,NX+1,4) = (alpha/(gamma1-1)+(1-alpha)/(gamma2-1)+1)*u*p+1/2*rho*u*(u^2+v^2);
        FBC(i,NX+1,5) = beta*rho*u;
        FBC(i,NX+1,6) = gamma1/(gamma1-1)*u*p*alpha+1/2*rho*u*(u^2+v^2)*beta;
    end
end

%Upper Boundary
if BC_A==4
    for j=1:NX
        rho   = Wfy(1,j,3);
        v     = Wfy(1,j,5);
        p     = Wfy(1,j,6);
        alpha = Wfy(1,j,8);

        NRK4 = 10;

        dv   = (vflux - v)/NRK4;

        for k=1:NRK4
        [dp_k1,drho_k1,dalpha_k1] = RK4_u(rho          ,p        ,alpha            ,dv,+1);
        [dp_k2,drho_k2,dalpha_k2] = RK4_u(rho+drho_k1/2,p+dp_k1/2,alpha+dalpha_k1/2,dv,+1);
        [dp_k3,drho_k3,dalpha_k3] = RK4_u(rho+drho_k2/2,p+dp_k2/2,alpha+dalpha_k2/2,dv,+1);
        [dp_k4,drho_k4,dalpha_k4] = RK4_u(rho+drho_k3  ,p+dp_k3  ,alpha+dalpha_k3  ,dv,+1);
        p     = p     + (dp_k1     + 2*dp_k2     + 2*dp_k3     + dp_k4    )/6;
        rho   = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
        alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
        end

        rho   = rhoflux;
        u     = uflux;
        v     = vflux;
        beta  = betaflux;
        alpha = alphaflux;

        GBC(1,j,1) = rho*v;
        GBC(1,j,2) = rho*v*u;
        GBC(1,j,3) = rho*v^2+p;
        GBC(1,j,4) = (alpha/(gamma1-1)+(1-alpha)/(gamma2-1)+1)*v*p+1/2*rho*v*(u^2+v^2);
        GBC(1,j,5) = beta*rho*v;
        GBC(1,j,6) = gamma1/(gamma1-1)*v*p*alpha+1/2*rho*v*(u^2+v^2)*beta;
    end
end

%Lower Boundary
if BC_B==4
    for j=1:NX
        rho   = Wfy(2*NY,j,3);
        v     = Wfy(2*NY,j,5);
        p     = Wfy(2*NY,j,6);
        alpha = Wfy(2*NY,j,8);

        NRK4 = 10;

        dv   = (vflux - v)/NRK4;

        for k=1:NRK4
        [dp_k1,drho_k1,dalpha_k1] = RK4_u(rho          ,p        ,alpha            ,dv,-1);
        [dp_k2,drho_k2,dalpha_k2] = RK4_u(rho+drho_k1/2,p+dp_k1/2,alpha+dalpha_k1/2,dv,-1);
        [dp_k3,drho_k3,dalpha_k3] = RK4_u(rho+drho_k2/2,p+dp_k2/2,alpha+dalpha_k2/2,dv,-1);
        [dp_k4,drho_k4,dalpha_k4] = RK4_u(rho+drho_k3  ,p+dp_k3  ,alpha+dalpha_k3  ,dv,-1);
        p     = p     + (dp_k1     + 2*dp_k2     + 2*dp_k3     + dp_k4    )/6;
        rho   = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
        alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
        end

        rho   = rhoflux;
        u     = uflux;
        v     = vflux;
        beta  = betaflux;
        alpha = alphaflux;

        GBC(NY+1,j,:) = rho*v;
        GBC(NY+1,j,:) = rho*v*u;
        GBC(NY+1,j,:) = rho*v^2+p;
        GBC(NY+1,j,:) = (alpha/(gamma1-1)+(1-alpha)/(gamma2-1)+1)*v*p+1/2*rho*v*(u^2+v^2);
        GBC(NY+1,j,:) = beta*rho*v;
        GBC(NY+1,j,:) = gamma1/(gamma1-1)*v*p*alpha+1/2*rho*v*(u^2+v^2)*beta;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Solid wall                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uflux = 0;
vflux = 0;

%Left Boundary
if BC_L==5
    for i=1:NY
        rho   = Wfx(i,1,3);
        u     = Wfx(i,1,4);
        v     = Wfx(i,1,5);
        p     = Wfx(i,1,6);
        beta  = Wfx(i,1,7);
        alpha = Wfx(i,1,8);

        NRK4 = 10;

        du   = (uflux - u)/NRK4;

        for k=1:NRK4
        [dp_k1,drho_k1,dalpha_k1] = RK4_u(rho          ,p        ,alpha            ,du,+1);
        [dp_k2,drho_k2,dalpha_k2] = RK4_u(rho+drho_k1/2,p+dp_k1/2,alpha+dalpha_k1/2,du,+1);
        [dp_k3,drho_k3,dalpha_k3] = RK4_u(rho+drho_k2/2,p+dp_k2/2,alpha+dalpha_k2/2,du,+1);
        [dp_k4,drho_k4,dalpha_k4] = RK4_u(rho+drho_k3  ,p+dp_k3  ,alpha+dalpha_k3  ,du,+1);
        p     = p     + (dp_k1     + 2*dp_k2     + 2*dp_k3     + dp_k4    )/6;
        rho   = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
        alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
        end

        u = uflux;

        FBC(i,1,1) = rho*u;
        FBC(i,1,2) = rho*u^2+p;
        FBC(i,1,3) = rho*u*v;
        FBC(i,1,4) = (alpha/(gamma1-1)+(1-alpha)/(gamma2-1)+1)*u*p+1/2*rho*u*(u^2+v^2);
        FBC(i,1,5) = beta*rho*u;
        FBC(i,1,6) = gamma1/(gamma1-1)*u*p*alpha+1/2*rho*u*(u^2+v^2)*beta;
    end
end

%Right Boundary
if BC_R==5
    for i=1:NY
        rho   = Wfx(i,2*NX,3);
        u     = Wfx(i,2*NX,4);
        v     = Wfx(i,2*NX,5);
        p     = Wfx(i,2*NX,6);
        beta  = Wfx(i,2*NX,7);
        alpha = Wfx(i,2*NX,8);

        NRK4 = 10;

        du   = (uflux - u)/NRK4;

        for k=1:NRK4
        [dp_k1,drho_k1,dalpha_k1] = RK4_u(rho          ,p        ,alpha            ,du,-1);
        [dp_k2,drho_k2,dalpha_k2] = RK4_u(rho+drho_k1/2,p+dp_k1/2,alpha+dalpha_k1/2,du,-1);
        [dp_k3,drho_k3,dalpha_k3] = RK4_u(rho+drho_k2/2,p+dp_k2/2,alpha+dalpha_k2/2,du,-1);
        [dp_k4,drho_k4,dalpha_k4] = RK4_u(rho+drho_k3  ,p+dp_k3  ,alpha+dalpha_k3  ,du,-1);
        p     = p     + (dp_k1     + 2*dp_k2     + 2*dp_k3     + dp_k4    )/6;
        rho   = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
        alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
        end

        u = uflux;

        FBC(i,NX+1,1) = rho*u;
        FBC(i,NX+1,2) = rho*u^2+p;
        FBC(i,NX+1,3) = rho*u*v;
        FBC(i,NX+1,4) = (alpha/(gamma1-1)+(1-alpha)/(gamma2-1)+1)*u*p+1/2*rho*u*(u^2+v^2);
        FBC(i,NX+1,5) = beta*rho*u;
        FBC(i,NX+1,6) = gamma1/(gamma1-1)*u*p*alpha+1/2*rho*u*(u^2+v^2)*beta;
    end
end

%Upper Boundary
if BC_A==5
    for j=1:NX
        rho   = Wfy(1,j,3);
        u     = Wfy(1,j,4);
        v     = Wfy(1,j,5);
        p     = Wfy(1,j,6);
        beta  = Wfy(1,j,7);
        alpha = Wfy(1,j,8);

        NRK4 = 10;

        dv   = (vflux - v)/NRK4;

        for k=1:NRK4
        [dp_k1,drho_k1,dalpha_k1] = RK4_u(rho          ,p        ,alpha            ,dv,+1);
        [dp_k2,drho_k2,dalpha_k2] = RK4_u(rho+drho_k1/2,p+dp_k1/2,alpha+dalpha_k1/2,dv,+1);
        [dp_k3,drho_k3,dalpha_k3] = RK4_u(rho+drho_k2/2,p+dp_k2/2,alpha+dalpha_k2/2,dv,+1);
        [dp_k4,drho_k4,dalpha_k4] = RK4_u(rho+drho_k3  ,p+dp_k3  ,alpha+dalpha_k3  ,dv,+1);
        p     = p     + (dp_k1     + 2*dp_k2     + 2*dp_k3     + dp_k4    )/6;
        rho   = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
        alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
        end

        v = vflux;

        GBC(1,j,1) = rho*v;
        GBC(1,j,2) = rho*v*u;
        GBC(1,j,3) = rho*v^2+p;
        GBC(1,j,4) = (alpha/(gamma1-1)+(1-alpha)/(gamma2-1)+1)*v*p+1/2*rho*v*(u^2+v^2);
        GBC(1,j,5) = beta*rho*v;
        GBC(1,j,6) = gamma1/(gamma1-1)*v*p*alpha+1/2*rho*v*(u^2+v^2)*beta;
    end
end

%Lower Boundary
if BC_B==5
    for j=1:NX
        rho   = Wfy(2*NY,j,3);
        u     = Wfy(2*NY,j,4);
        v     = Wfy(2*NY,j,5);
        p     = Wfy(2*NY,j,6);
        beta  = Wfy(2*NY,j,7);
        alpha = Wfy(2*NY,j,8);

        NRK4 = 10;

        dv   = (vflux - v)/NRK4;

        for k=1:NRK4
        [dp_k1,drho_k1,dalpha_k1] = RK4_u(rho          ,p        ,alpha            ,dv,-1);
        [dp_k2,drho_k2,dalpha_k2] = RK4_u(rho+drho_k1/2,p+dp_k1/2,alpha+dalpha_k1/2,dv,-1);
        [dp_k3,drho_k3,dalpha_k3] = RK4_u(rho+drho_k2/2,p+dp_k2/2,alpha+dalpha_k2/2,dv,-1);
        [dp_k4,drho_k4,dalpha_k4] = RK4_u(rho+drho_k3  ,p+dp_k3  ,alpha+dalpha_k3  ,dv,-1);
        p     = p     + (dp_k1     + 2*dp_k2     + 2*dp_k3     + dp_k4    )/6;
        rho   = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
        alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
        end

        v = vflux;

        GBC(NY+1,j,1) = rho*v;
        GBC(NY+1,j,2) = rho*v*u;
        GBC(NY+1,j,3) = rho*v^2+p;
        GBC(NY+1,j,4) = (alpha/(gamma1-1)+(1-alpha)/(gamma2-1)+1)*v*p+1/2*rho*v*(u^2+v^2);
        GBC(NY+1,j,5) = beta*rho*v;
        GBC(NY+1,j,6) = gamma1/(gamma1-1)*v*p*alpha+1/2*rho*v*(u^2+v^2)*beta;
    end
end



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [du,drho,dalpha] = RK4_p(rho,p,alpha,dp,lr)
%------------------------------------------------%
% RK4 step for left wave                         %
%------------------------------------------------%

global gamma1 gamma2

a      = sqrt(1/(alpha/gamma1+(1-alpha)/gamma2)*p/rho);
b      = alpha*(1-alpha)*(gamma1-gamma2)/((1-alpha)*gamma1+alpha*gamma2);
du     = lr*dp/(rho*a);
dalpha =    dp*b/(rho*a^2);
drho   =    dp/(a^2);

end % function RK4_p


function [dp,drho,dalpha]=RK4_u(rho,p,alpha,du,lr)

global gamma1 gamma2

a      = sqrt(1/(alpha/gamma1+(1-alpha)/gamma2)*p/rho);
b      = alpha*(1-alpha)*(gamma1-gamma2)/((1-alpha)*gamma1+alpha*gamma2);
dp     = lr*du*rho*a;
drho   = lr*du*rho/a;
dalpha = lr*du*b/a;

end % function RK4_u