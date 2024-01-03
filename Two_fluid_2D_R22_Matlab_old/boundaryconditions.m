function [FBC,GBC]=boundaryconditions(NX,NY,W,Wfx,Wfy,gamma1,gamma2)

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
    FBC(:,1,2) = W(:,1,3).*W(:,1,4).^(2e0)+W(:,1,6);
    FBC(:,1,3) = W(:,1,3).*W(:,1,4).*W(:,1,5);
    FBC(:,1,4) = (W(:,1,8)/(gamma1-1e0)+(1e0-W(:,1,8))/(gamma2-1e0)+1e0).*W(:,1,4).*W(:,1,6)+0.5e0*W(:,1,3).*W(:,1,4).*(W(:,1,4).^(2e0)+W(:,1,5).^(2e0));
    FBC(:,1,5) = W(:,1,7).*W(:,1,3).*W(:,1,4);
    FBC(:,1,6) = gamma1/(gamma1-1e0)*W(:,1,4).*W(:,1,6).*W(:,1,8)+1e0/2e0*W(:,1,3).*W(:,1,4).*(W(:,1,4).^(2e0)+W(:,1,5).^(2e0)).*W(:,1,8);
end

%Right Boundary
if BC_R==0
    FBC(:,NX+1,1) = W(:,NX,3).*W(:,NX,4);
    FBC(:,NX+1,2) = W(:,NX,3).*W(:,NX,4).^(2e0)+W(:,NX,6);
    FBC(:,NX+1,3) = W(:,NX,3).*W(:,NX,4).*W(:,NX,5);
    FBC(:,NX+1,4) = (W(:,NX,8)/(gamma1-1e0)+(1e0-W(:,NX,8))/(gamma2-1e0)+1e0).*W(:,NX,4).*W(:,NX,6)+0.5e0*W(:,NX,3).*W(:,NX,4).*(W(:,NX,4).^(2e0)+W(:,NX,5).^(2e0));
    FBC(:,NX+1,5) = W(:,NX,7).*W(:,NX,3).*W(:,NX,4);
    FBC(:,NX+1,6) = gamma1/(gamma1-1e0).*W(:,NX,4).*W(:,NX,6).*W(:,NX,8)+1e0/2e0*W(:,NX,3).*W(:,NX,4).*(W(:,NX,4).^(2e0)+W(:,NX,5)^(2e0)).*W(:,NX,8);
end

%Upper Boundary
if BC_A==0 
    GBC(1,:,1) = W(1,:,3).*W(1,:,5);
    GBC(1,:,2) = W(1,:,3).*W(1,:,5).*W(1,:,4);
    GBC(1,:,3) = W(1,:,3).*W(1,:,5).^(2e0)+W(1,:,6);
    GBC(1,:,4) = (W(1,:,8)/(gamma1-1e0)+(1e0-W(1,:,8))/(gamma2-1e0)+1e0).*W(1,:,5).*W(1,:,6)+0.5e0*W(1,:,3).*W(1,:,5).*(W(1,:,4).^(2e0)+W(1,:,5).^(2e0));
    GBC(1,:,5) = W(1,:,7).*W(1,:,3).*W(1,:,5);
    GBC(1,:,6) = gamma1/(gamma1-1e0)*W(1,:,5).*W(1,:,6).*W(1,:,8)+1e0/2e0*W(1,:,3).*W(1,:,5).*(W(1,:,4).^(2e0)+W(1,:,5).^(2e0)).*W(1,:,8);
end

%Lower Boundary
if BC_B==0
    GBC(NY+1,:,1) = W(NY,:,3).*W(NY,:,5);
    GBC(NY+1,:,2) = W(NY,:,3).*W(NY,:,5).*W(NY,:,4);
    GBC(NY+1,:,3) = W(NY,:,3).*W(NY,:,5).^(2e0)+W(NY,:,6);
    GBC(NY+1,:,4) = (W(NY,:,8)/(gamma1-1e0)+(1e0-W(NY,:,8))/(gamma2-1e0)+1e0).*W(NY,:,5).*W(NY,:,6)+0.5e0*W(NY,:,3).*W(NY,:,5).*(W(NY,:,4).^(2e0)+W(NY,:,5).^(2e0));
    GBC(NY+1,:,5) = W(NY,:,7).*W(NY,:,3).*W(NY,:,5);
    GBC(NY+1,:,6) = gamma1/(gamma1-1e0)*W(NY,:,5).*W(NY,:,6).*W(NY,:,8)+1e0/2e0*W(NY,:,3).*W(NY,:,5).*(W(NY,:,4).^(2e0)+W(NY,:,5).^(2e0)).*W(NY,:,8);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Supersonic outflow                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Left Boundary
if (BC_L==1) 
    FBC(:,1,1) = Wfx(:,1,3).*Wfx(:,1,4);
    FBC(:,1,2) = Wfx(:,1,3).*Wfx(:,1,4).^(2e0)+Wfx(:,1,6);
    FBC(:,1,3) = Wfx(:,1,3).*Wfx(:,1,4).*Wfx(:,1,5);
    FBC(:,1,4) = (Wfx(:,1,8)/(gamma1-1e0)+(1e0-Wfx(:,1,8))/(gamma2-1e0)+1e0).*Wfx(:,1,4).*Wfx(:,1,6)+0.5e0*Wfx(:,1,3).*Wfx(:,1,4).*(Wfx(:,1,4).^(2e0)+Wfx(:,1,5).^(2e0));
    FBC(:,1,5) = Wfx(:,1,7).*Wfx(:,1,3).*Wfx(:,1,4);
    FBC(:,1,6) = gamma1/(gamma1-1e0)*Wfx(:,1,4).*Wfx(:,1,6).*Wfx(:,1,8)+1e0/2e0*Wfx(:,1,3).*Wfx(:,1,4).*(Wfx(:,1,4).^(2e0)+Wfx(:,1,5).^(2e0)).*Wfx(:,1,8);
end

%Right Boundary
if (BC_R==1) 
    FBC(:,NX+1,1) = Wfx(:,2*NX,3).*Wfx(:,2*NX,4);
    FBC(:,NX+1,2) = Wfx(:,2*NX,3).*Wfx(:,2*NX,4).^(2e0)+Wfx(:,2*NX,6);
    FBC(:,NX+1,3) = Wfx(:,2*NX,3).*Wfx(:,2*NX,4).*Wfx(:,2*NX,5);
    FBC(:,NX+1,4) = (Wfx(:,2*NX,8)/(gamma1-1e0)+(1e0-Wfx(:,2*NX,8))/(gamma2-1e0)+1e0)*Wfx(:,2*NX,4).*Wfx(:,2*NX,6)+0.5e0*Wfx(:,2*NX,3).*Wfx(:,2*NX,4).*(Wfx(:,2*NX,4).^(2e0)+Wfx(:,2*NX,5).^(2e0));
    FBC(:,NX+1,5) = Wfx(:,2*NX,7).*Wfx(:,2*NX,3).*Wfx(:,2*NX,4);
    FBC(:,NX+1,6) = gamma1/(gamma1-1e0)*Wfx(:,2*NX,4).*Wfx(:,2*NX,6).*Wfx(:,2*NX,8)+1e0/2e0*Wfx(:,2*NX,3).*Wfx(:,2*NX,4).*(Wfx(:,2*NX,4).^(2e0)+Wfx(:,2*NX,5).^(2e0)).*Wfx(:,2*NX,8);
end

%Upper Boundary
if (BC_A==1) 
    GBC(1,:,1) = Wfy(1,:,3).*Wfy(1,:,5);
    GBC(1,:,2) = Wfy(1,:,3).*Wfy(1,:,5).*Wfy(1,:,4);
    GBC(1,:,3) = Wfy(1,:,3).*Wfy(1,:,5).^(2e0)+Wfy(1,:,6);
    GBC(1,:,4) = (Wfy(1,:,8)/(gamma1-1e0)+(1e0-Wfy(1,:,8))/(gamma2-1e0)+1e0).*Wfy(1,:,5).*Wfy(1,:,6)+0.5e0*Wfy(1,:,3).*Wfy(1,:,5).*(Wfy(1,:,4).^(2e0)+Wfy(1,:,5).^(2e0));
    GBC(1,:,5) = Wfy(1,:,7)*Wfy(1,:,3)*Wfy(1,:,5);
    GBC(1,:,6) = gamma1/(gamma1-1e0)*Wfy(1,:,5).*Wfy(1,:,6).*Wfy(1,:,8)+1e0/2e0*Wfy(1,:,3).*Wfy(1,:,5).*(Wfy(1,:,4).^(2e0)+Wfy(1,:,5).^(2e0)).*Wfy(1,:,8);
end

%Lower Boundary
if (BC_B==1) 
    GBC(NY+1,:,1) = Wfy(2*NY,:,3).*Wfy(2*NY,:,5);
    GBC(NY+1,:,2) = Wfy(2*NY,:,3).*Wfy(2*NY,:,5).*Wfy(2*NY,:,4);
    GBC(NY+1,:,3) = Wfy(2*NY,:,3).*Wfy(2*NY,:,5).^(2e0)+Wfy(2*NY,:,6);
    GBC(NY+1,:,4) = (Wfy(2*NY,:,8)/(gamma1-1e0)+(1e0-Wfy(2*NY,:,8))/(gamma2-1e0)+1e0).*Wfy(2*NY,:,5).*Wfy(2*NY,:,6)+0.5e0*Wfy(2*NY,:,3).*Wfy(2*NY,:,5).*(Wfy(2*NY,:,4).^(2e0)+Wfy(2*NY,:,5).^(2e0));
    GBC(NY+1,:,5) = Wfy(2*NY,:,7).*Wfy(2*NY,:,3).*Wfy(2*NY,:,5);
    GBC(NY+1,:,6) = gamma1/(gamma1-1e0)*Wfy(2*NY,:,5).*Wfy(2*NY,:,6).*Wfy(2*NY,:,8)+1e0/2e0*Wfy(2*NY,:,3).*Wfy(2*NY,:,5).*(Wfy(2*NY,:,4).^(2e0)+Wfy(2*NY,:,5).^(2e0)).*Wfy(2*NY,:,8);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Supersonic inflow                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho   = 1.92691e0;
u     = -0.33361e0;
v     = 0e0;
p     = 1.56980e0;
beta  = 1e0;
alpha = 1e0;

%Left Boundary
if BC_L==2
    FBC(:,1,1) = rho*u;
    FBC(:,1,2) = rho*u^(2e0)+p;
    FBC(:,1,3) = rho*u*v;
    FBC(:,1,4) = (alpha/(gamma1-1e0)+(1e0-alpha)/(gamma2-1e0)+1e0)*u*p+0.5e0*rho*u*(u^(2e0)+v^(2e0));
    FBC(:,1,5) = beta*rho*u;
    FBC(:,1,6) = gamma1/(gamma1-1e0)*u*p*alpha+1e0/2e0*rho*u*(u^(2e0)+v^(2e0))*beta;
end

%Right Boundary
if BC_R==2
    FBC(:,NX+1,1) = rho*u;
    FBC(:,NX+1,2) = rho*u^(2e0)+p;
    FBC(:,NX+1,3) = rho*u*v;
    FBC(:,NX+1,4) = (alpha/(gamma1-1e0)+(1e0-alpha)/(gamma2-1e0)+1e0)*u*p+0.5e0*rho*u*(u^(2e0)+v^(2e0));
    FBC(:,NX+1,5) = beta*rho*u;
    FBC(:,NX+1,6) = gamma1/(gamma1-1e0)*u*p*alpha+1e0/2e0*rho*u*(u^(2e0)+v^(2e0))*beta;
end

%Upper Boundary
if (BC_A==2) 
    GBC(1,:,1)  = rho*v;
    GBC(1,:,2)  = rho*v*u;
    GBC(1,:,3)  = rho*v^(2e0)+p;
    GBC(1,:,4)  = (alpha/(gamma1-1e0)+(1e0-alpha)/(gamma2-1e0)+1e0)*v*p+0.5e0*rho*v*(u^(2e0)+v^(2e0));
    GBC(1,:,5)  = beta*rho*v;
    GBC(1,:,6)  = gamma1/(gamma1-1e0)*v*p*alpha+1e0/2e0*rho*v*(u^(2e0)+v^(2e0))*beta;
end

%Lower Boundary
if (BC_B==2) 
    GBC(NY+1,:,1)  = rho*v;
    GBC(NY+1,:,2)  = rho*v*u;
    GBC(NY+1,:,3)  = rho*v^(2e0)+p;
    GBC(NY+1,:,4)  = (alpha/(gamma1-1e0)+(1e0-alpha)/(gamma2-1e0)+1e0)*v*p+0.5e0*rho*v*(u^(2e0)+v^(2e0));
    GBC(NY+1,:,5)  = beta*rho*v;
    GBC(NY+1,:,6)  = gamma1/(gamma1-1e0)*v*p*alpha+1e0/2e0*rho*v*(u^(2e0)+v^(2e0))*beta;
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

        NRK4 = 500;

        dpr   = (pflux - p)/NRK4;

        for k=1:NRK4
%             CALL RK4_k1_R_ux(rho,p,alpha,gamma1,gamma2,dpr,u_k1,rho_k1,alpha_k1)
%             CALL RK4_k2_R_ux(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dpr,u_k2,rho_k2,alpha_k2)
%             CALL RK4_k3_R_ux(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dpr,u_k3,rho_k3,alpha_k3)
%             CALL RK4_k4_R_ux(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dpr,u_k4,rho_k4,alpha_k4)
            u = u + (u_k1 + 2e0*u_k2 + 2e0*u_k3 + u_k4)/6e0;
            rho = rho + (rho_k1 + 2e0*rho_k2 + 2e0*rho_k3 + rho_k4)/6e0;
            alpha = alpha + (alpha_k1 + 2e0*alpha_k2 + 2e0*alpha_k3 + alpha_k4)/6e0;
        end

        p = pflux;

        FBC(i,1,1) = rho*u;
        FBC(i,1,2) = rho*u^(2e0)+p;
        FBC(i,1,3) = rho*u*v;
        FBC(i,1,4) = (alpha/(gamma1-1e0)+(1e0-alpha)/(gamma2-1e0)+1e0)*u*p+0.5e0*rho*u*(u^(2e0)+v^(2e0));
        FBC(i,1,5) = beta*rho*u;
        FBC(i,1,6) = gamma1/(gamma1-1e0)*u*p*alpha+1e0/2e0*rho*u*(u^(2e0)+v^(2e0))*beta;
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

        NRK4 = 500;

        dpl   = (pflux - p)/NRK4;

        for k=1:NRK4
        %     CALL RK4_k1_L_ux(rho,p,alpha,gamma1,gamma2,dpl,u_k1,rho_k1,alpha_k1)
        %     CALL RK4_k2_L_ux(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dpl,u_k2,rho_k2,alpha_k2)
        %     CALL RK4_k3_L_ux(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dpl,u_k3,rho_k3,alpha_k3)
        %     CALL RK4_k4_L_ux(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dpl,u_k4,rho_k4,alpha_k4)
            u = u + (u_k1 + 2e0*u_k2 + 2e0*u_k3 + u_k4)/6e0;
            rho = rho + (rho_k1 + 2e0*rho_k2 + 2e0*rho_k3 + rho_k4)/6e0;
            alpha = alpha + (alpha_k1 + 2e0*alpha_k2 + 2e0*alpha_k3 + alpha_k4)/6e0;
        end

        p = pflux;

        FBC(i,NX+1,1) = rho*u;
        FBC(i,NX+1,2) = rho*u^(2e0)+p;
        FBC(i,NX+1,3) = rho*u*v;
        FBC(i,NX+1,4) = (alpha/(gamma1-1e0)+(1e0-alpha)/(gamma2-1e0)+1e0)*u*p+0.5e0*rho*u*(u^(2e0)+v^(2e0));
        FBC(i,NX+1,5) = beta*rho*u;
        FBC(i,NX+1,6) = gamma1/(gamma1-1e0)*u*p*alpha+1e0/2e0*rho*u*(u^(2e0)+v^(2e0))*beta;
    end
end
          
%Upper Boundary
if (BC_A==3) 
    for j=1:NX
    rho   = Wfy(1,j,3);
    u     = Wfy(1,j,4);
    v     = Wfy(1,j,5);
    p     = Wfy(1,j,6);
    beta  = Wfy(1,j,7);
    alpha = Wfy(1,j,8);

    NRK4 = 500;

    dpr   = (pflux - p)/NRK4;

    for k=1:NRK4
%         CALL RK4_k1_R_vy(rho,p,alpha,gamma1,gamma2,dpr,v_k1,rho_k1,alpha_k1)
%         CALL RK4_k2_R_vy(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dpr,v_k2,rho_k2,alpha_k2)
%         CALL RK4_k3_R_vy(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dpr,v_k3,rho_k3,alpha_k3)
%         CALL RK4_k4_R_vy(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dpr,v_k4,rho_k4,alpha_k4)
        v = v + (v_k1 + 2e0*v_k2 + 2e0*v_k3 + v_k4)/6e0;
        rho = rho + (rho_k1 + 2e0*rho_k2 + 2e0*rho_k3 + rho_k4)/6e0;
        alpha = alpha + (alpha_k1 + 2e0*alpha_k2 + 2e0*alpha_k3 + alpha_k4)/6e0;
    end

    p = pflux;

    GBC(1,j,1) = rho*v;
    GBC(1,j,2) = rho*v*u;
    GBC(1,j,3) = rho*v^(2e0)+p;
    GBC(1,j,4) = (alpha/(gamma1-1e0)+(1e0-alpha)/(gamma2-1e0)+1e0)*v*p+0.5e0*rho*v*(u^(2e0)+v^(2e0));
    GBC(1,j,5) = beta*rho*v;
    GBC(1,j,6) = gamma1/(gamma1-1e0)*v*p*alpha+1e0/2e0*rho*v*(u^(2e0)+v^(2e0))*beta;
    end
end

%Lower Boundary
if (BC_B==3) 
    for j=1:NX
    rho   = Wfy(2*NY,j,3);
    u     = Wfy(2*NY,j,4);
    v     = Wfy(2*NY,j,5);
    p     = Wfy(2*NY,j,6);
    beta  = Wfy(2*NY,j,7);
    alpha = Wfy(2*NY,j,8);

    NRK4 = 500;

    dpl   = (pflux - p)/NRK4;

    for k=1:NRK4
%         CALL RK4_k1_L_vy(rho,p,alpha,gamma1,gamma2,dpl,v_k1,rho_k1,alpha_k1)
%         CALL RK4_k2_L_vy(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dpl,v_k2,rho_k2,alpha_k2)
%         CALL RK4_k3_L_vy(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dpl,v_k3,rho_k3,alpha_k3)
%         CALL RK4_k4_L_vy(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dpl,v_k4,rho_k4,alpha_k4)
        v = v + (v_k1 + 2e0*v_k2 + 2e0*v_k3 + v_k4)/6e0;
        rho = rho + (rho_k1 + 2e0*rho_k2 + 2e0*rho_k3 + rho_k4)/6e0;
        alpha = alpha + (alpha_k1 + 2e0*alpha_k2 + 2e0*alpha_k3 + alpha_k4)/6e0;
    end

    p = pflux;

    GBC(NY+1,j,1) = rho*v;
    GBC(NY+1,j,2) = rho*v*u;
    GBC(NY+1,j,3) = rho*v^(2e0)+p;
    GBC(NY+1,j,4) = (alpha/(gamma1-1e0)+(1e0-alpha)/(gamma2-1e0)+1e0)*v*p+0.5e0*rho*v*(u^(2e0)+v^(2e0));
    GBC(NY+1,j,5) = beta*rho*v;
    GBC(NY+1,j,6) = gamma1/(gamma1-1e0)*v*p*alpha+1e0/2e0*rho*v*(u^(2e0)+v^(2e0))*beta;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Subsonic inflow                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhoflux   = 1.92691e0;
uflux     = -0.33361e0;
vflux     = 0e0;
betaflux  = 1e0;
alphaflux = 1e0;

%Left Boundary
if BC_L==4
    for i=1:NY
        rho   = Wfx(i,1,3);
        u     = Wfx(i,1,4);
        p     = Wfx(i,1,6);
        alpha = Wfx(i,1,8);

        NRK4 = 500;

        dur   = (uflux - u)/NRK4;

        for k=1:NRK4
        %     CALL RK4_k1_R_x(rho,p,alpha,gamma1,gamma2,dur,p_k1,rho_k1,alpha_k1)
        %     CALL RK4_k2_R_x(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dur,p_k2,rho_k2,alpha_k2)
        %     CALL RK4_k3_R_x(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dur,p_k3,rho_k3,alpha_k3)
        %     CALL RK4_k4_R_x(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dur,p_k4,rho_k4,alpha_k4)
            p = p + (p_k1 + 2e0*p_k2 + 2e0*p_k3 + p_k4)/6e0;
            rho = rho + (rho_k1 + 2e0*rho_k2 + 2e0*rho_k3 + rho_k4)/6e0;
            alpha = alpha + (alpha_k1 + 2e0*alpha_k2 + 2e0*alpha_k3 + alpha_k4)/6e0;
        end

        rho   = rhoflux;
        u     = uflux;
        v     = vflux;
        beta  = betaflux;
        alpha = alphaflux;

        FBC(i,1,1) = rho*u;
        FBC(i,1,2) = rho*u^(2e0)+p;
        FBC(i,1,3) = rho*u*v;
        FBC(i,1,4) = (alpha/(gamma1-1e0)+(1e0-alpha)/(gamma2-1e0)+1e0)*u*p+0.5e0*rho*u*(u^(2e0)+v^(2e0));
        FBC(i,1,5) = beta*rho*u;
        FBC(i,1,6) = gamma1/(gamma1-1e0)*u*p*alpha+1e0/2e0*rho*u*(u^(2e0)+v^(2e0))*beta;
    end
end

%Right Boundary
if (BC_R==4) 
    for i=1:NY
        rho   = Wfx(i,2*NX,3);
        u     = Wfx(i,2*NX,4);
        p     = Wfx(i,2*NX,6);
        alpha = Wfx(i,2*NX,8);

        NRK4 = 500;

        dul   = (uflux - u)/NRK4;

        for k=1:NRK4
        %     CALL RK4_k1_L_x(rho,p,alpha,gamma1,gamma2,dul,p_k1,rho_k1,alpha_k1)
        %     CALL RK4_k2_L_x(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dul,p_k2,rho_k2,alpha_k2)
        %     CALL RK4_k3_L_x(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dul,p_k3,rho_k3,alpha_k3)
        %     CALL RK4_k4_L_x(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dul,p_k4,rho_k4,alpha_k4)
            p = p + (p_k1 + 2e0*p_k2 + 2e0*p_k3 + p_k4)/6e0;
            rho = rho + (rho_k1 + 2e0*rho_k2 + 2e0*rho_k3 + rho_k4)/6e0;
            alpha = alpha + (alpha_k1 + 2e0*alpha_k2 + 2e0*alpha_k3 + alpha_k4)/6e0;
        end

        rho   = rhoflux;
        u     = uflux;
        v     = vflux;
        beta  = betaflux;
        alpha = alphaflux;

        FBC(i,NX+1,1) = rho*u;
        FBC(i,NX+1,2) = rho*u^(2e0)+p;
        FBC(i,NX+1,3) = rho*u*v;
        FBC(i,NX+1,4) = (alpha/(gamma1-1e0)+(1e0-alpha)/(gamma2-1e0)+1e0)*u*p+0.5e0*rho*u*(u^(2e0)+v^(2e0));
        FBC(i,NX+1,5) = beta*rho*u;
        FBC(i,NX+1,6) = gamma1/(gamma1-1e0)*u*p*alpha+1e0/2e0*rho*u*(u^(2e0)+v^(2e0))*beta;
    end
end

%Upper Boundary
if (BC_A==4) 
    for j=1:NX
        rho   = Wfy(1,j,3);
        v     = Wfy(1,j,5);
        p     = Wfy(1,j,6);
        alpha = Wfy(1,j,8);

        NRK4 = 500;

        dvr   = (vflux - v)/NRK4;

        for k=1:NRK4
        %     CALL RK4_k1_R_y(rho,p,alpha,gamma1,gamma2,dvr,p_k1,rho_k1,alpha_k1)
        %     CALL RK4_k2_R_y(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dvr,p_k2,rho_k2,alpha_k2)
        %     CALL RK4_k3_R_y(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dvr,p_k3,rho_k3,alpha_k3)
        %     CALL RK4_k4_R_y(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dvr,p_k4,rho_k4,alpha_k4)
            p = p + (p_k1 + 2e0*p_k2 + 2e0*p_k3 + p_k4)/6e0;
            rho = rho + (rho_k1 + 2e0*rho_k2 + 2e0*rho_k3 + rho_k4)/6e0;
            alpha = alpha + (alpha_k1 + 2e0*alpha_k2 + 2e0*alpha_k3 + alpha_k4)/6e0;
        end

        rho   = rhoflux;
        u     = uflux;
        v     = vflux;
        beta  = betaflux;
        alpha = alphaflux;

        GBC(1,j,1) = rho*v;
        GBC(1,j,2) = rho*v*u;
        GBC(1,j,3) = rho*v^(2e0)+p;
        GBC(1,j,4) = (alpha/(gamma1-1e0)+(1e0-alpha)/(gamma2-1e0)+1e0)*v*p+0.5e0*rho*v*(u^(2e0)+v^(2e0));
        GBC(1,j,5) = beta*rho*v;
        GBC(1,j,6) = gamma1/(gamma1-1e0)*v*p*alpha+1e0/2e0*rho*v*(u^(2e0)+v^(2e0))*beta;
    end
end

%Lower Boundary
if (BC_B==4) 
    for j=1:NX
        rho   = Wfy(2*NY,j,3);
        v     = Wfy(2*NY,j,5);
        p     = Wfy(2*NY,j,6);
        alpha = Wfy(2*NY,j,8);

        NRK4 = 500;

        dvl   = (vflux - v)/NRK4;

        for k=1:NRK4
        %     CALL RK4_k1_L_y(rho,p,alpha,gamma1,gamma2,dvl,p_k1,rho_k1,alpha_k1)
        %     CALL RK4_k2_L_y(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dvl,p_k2,rho_k2,alpha_k2)
        %     CALL RK4_k3_L_y(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dvl,p_k3,rho_k3,alpha_k3)
        %     CALL RK4_k4_L_y(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dvl,p_k4,rho_k4,alpha_k4)
            p = p + (p_k1 + 2e0*p_k2 + 2e0*p_k3 + p_k4)/6e0;
            rho = rho + (rho_k1 + 2e0*rho_k2 + 2e0*rho_k3 + rho_k4)/6e0;
            alpha = alpha + (alpha_k1 + 2e0*alpha_k2 + 2e0*alpha_k3 + alpha_k4)/6e0;
        end

        rho   = rhoflux;
        u     = uflux;
        v     = vflux;
        beta  = betaflux;
        alpha = alphaflux;

        GBC(NY+1,j,:) = rho*v;
        GBC(NY+1,j,:) = rho*v*u;
        GBC(NY+1,j,:) = rho*v^(2e0)+p;
        GBC(NY+1,j,:) = (alpha/(gamma1-1e0)+(1e0-alpha)/(gamma2-1e0)+1e0)*v*p+0.5e0*rho*v*(u^(2e0)+v^(2e0));
        GBC(NY+1,j,:) = beta*rho*v;
        GBC(NY+1,j,:) = gamma1/(gamma1-1e0)*v*p*alpha+1e0/2e0*rho*v*(u^(2e0)+v^(2e0))*beta;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Solid wall                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uflux = 0e0;
vflux = 0e0;

%Left Boundary
if (BC_L==5) 
    for i=1:NY
        rho   = Wfx(i,1,3);
        u     = Wfx(i,1,4);
        v     = Wfx(i,1,5);
        p     = Wfx(i,1,6);
        beta  = Wfx(i,1,7);
        alpha = Wfx(i,1,8);

        NRK4 = 500;

        dur   = (uflux - u)/NRK4;

        for k=1:NRK4
        %     CALL RK4_k1_R_x(rho,p,alpha,gamma1,gamma2,dur,p_k1,rho_k1,alpha_k1)
        %     CALL RK4_k2_R_x(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dur,p_k2,rho_k2,alpha_k2)
        %     CALL RK4_k3_R_x(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dur,p_k3,rho_k3,alpha_k3)
        %     CALL RK4_k4_R_x(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dur,p_k4,rho_k4,alpha_k4)
            p = p + (p_k1 + 2e0*p_k2 + 2e0*p_k3 + p_k4)/6e0;
            rho = rho + (rho_k1 + 2e0*rho_k2 + 2e0*rho_k3 + rho_k4)/6e0;
            alpha = alpha + (alpha_k1 + 2e0*alpha_k2 + 2e0*alpha_k3 + alpha_k4)/6e0;
        end

        u = uflux;

        FBC(i,1,1) = rho*u;
        FBC(i,1,2) = rho*u^(2e0)+p;
        FBC(i,1,3) = rho*u*v;
        FBC(i,1,4) = (alpha/(gamma1-1e0)+(1e0-alpha)/(gamma2-1e0)+1e0)*u*p+0.5e0*rho*u*(u^(2e0)+v^(2e0));
        FBC(i,1,5) = beta*rho*u;
        FBC(i,1,6) = gamma1/(gamma1-1e0)*u*p*alpha+1e0/2e0*rho*u*(u^(2e0)+v^(2e0))*beta;
    end
end

%Right Boundary
if (BC_R==5) 
    for i=1:NY
        rho   = Wfx(i,2*NX,3);
        u     = Wfx(i,2*NX,4);
        v     = Wfx(i,2*NX,5);
        p     = Wfx(i,2*NX,6);
        beta  = Wfx(i,2*NX,7);
        alpha = Wfx(i,2*NX,8);

        NRK4 = 500;

        dul   = (uflux - u)/NRK4;

        for k=1:NRK4
        %     CALL RK4_k1_L_x(rho,p,alpha,gamma1,gamma2,dul,p_k1,rho_k1,alpha_k1)
        %     CALL RK4_k2_L_x(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dul,p_k2,rho_k2,alpha_k2)
        %     CALL RK4_k3_L_x(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dul,p_k3,rho_k3,alpha_k3)
        %     CALL RK4_k4_L_x(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dul,p_k4,rho_k4,alpha_k4)
            p = p + (p_k1 + 2e0*p_k2 + 2e0*p_k3 + p_k4)/6e0;
            rho = rho + (rho_k1 + 2e0*rho_k2 + 2e0*rho_k3 + rho_k4)/6e0;
            alpha = alpha + (alpha_k1 + 2e0*alpha_k2 + 2e0*alpha_k3 + alpha_k4)/6e0;
        end

        u = uflux;

        FBC(i,NX+1,1) = rho*u;
        FBC(i,NX+1,2) = rho*u^(2e0)+p;
        FBC(i,NX+1,3) = rho*u*v;
        FBC(i,NX+1,4) = (alpha/(gamma1-1e0)+(1e0-alpha)/(gamma2-1e0)+1e0)*u*p+0.5e0*rho*u*(u^(2e0)+v^(2e0));
        FBC(i,NX+1,5) = beta*rho*u;
        FBC(i,NX+1,6) = gamma1/(gamma1-1e0)*u*p*alpha+1e0/2e0*rho*u*(u^(2e0)+v^(2e0))*beta;
    end
end


%Upper Boundary
if (BC_A==5) 
    for j=1:NX
        rho   = Wfy(1,j,3);
        u     = Wfy(1,j,4);
        v     = Wfy(1,j,5);
        p     = Wfy(1,j,6);
        beta  = Wfy(1,j,7);
        alpha = Wfy(1,j,8);

        NRK4 = 500;

        dvr   = (vflux - v)/NRK4;

        for k=1:NRK4
        %     CALL RK4_k1_R_y(rho,p,alpha,gamma1,gamma2,dvr,p_k1,rho_k1,alpha_k1)
        %     CALL RK4_k2_R_y(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dvr,p_k2,rho_k2,alpha_k2)
        %     CALL RK4_k3_R_y(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dvr,p_k3,rho_k3,alpha_k3)
        %     CALL RK4_k4_R_y(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dvr,p_k4,rho_k4,alpha_k4)
            p = p + (p_k1 + 2e0*p_k2 + 2e0*p_k3 + p_k4)/6e0;
            rho = rho + (rho_k1 + 2e0*rho_k2 + 2e0*rho_k3 + rho_k4)/6e0;
            alpha = alpha + (alpha_k1 + 2e0*alpha_k2 + 2e0*alpha_k3 + alpha_k4)/6e0;
        end

        v = vflux;

        GBC(1,j,1) = rho*v;
        GBC(1,j,2) = rho*v*u;
        GBC(1,j,3) = rho*v^(2e0)+p;
        GBC(1,j,4) = (alpha/(gamma1-1e0)+(1e0-alpha)/(gamma2-1e0)+1e0)*v*p+0.5e0*rho*v*(u^(2e0)+v^(2e0));
        GBC(1,j,5) = beta*rho*v;
        GBC(1,j,6) = gamma1/(gamma1-1e0)*v*p*alpha+1e0/2e0*rho*v*(u^(2e0)+v^(2e0))*beta;
    end
end

%Lower Boundary
if (BC_B==5) 
    for j=1:NX
        rho   = Wfy(2*NY,j,3);
        u     = Wfy(2*NY,j,4);
        v     = Wfy(2*NY,j,5);
        p     = Wfy(2*NY,j,6);
        beta  = Wfy(2*NY,j,7);
        alpha = Wfy(2*NY,j,8);

        NRK4 = 500;

        dvl   = (vflux - v)/NRK4;

        for k=1:NRK4
        %     CALL RK4_k1_L_y(rho,p,alpha,gamma1,gamma2,dvl,p_k1,rho_k1,alpha_k1)
        %     CALL RK4_k2_L_y(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dvl,p_k2,rho_k2,alpha_k2)
        %     CALL RK4_k3_L_y(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dvl,p_k3,rho_k3,alpha_k3)
        %     CALL RK4_k4_L_y(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dvl,p_k4,rho_k4,alpha_k4)
            p = p + (p_k1 + 2e0*p_k2 + 2e0*p_k3 + p_k4)/6e0;
            rho = rho + (rho_k1 + 2e0*rho_k2 + 2e0*rho_k3 + rho_k4)/6e0;
            alpha = alpha + (alpha_k1 + 2e0*alpha_k2 + 2e0*alpha_k3 + alpha_k4)/6e0;
        end

        v = vflux;

        GBC(NY+1,j,1) = rho*v;
        GBC(NY+1,j,2) = rho*v*u;
        GBC(NY+1,j,3) = rho*v^(2e0)+p;
        GBC(NY+1,j,4) = (alpha/(gamma1-1e0)+(1e0-alpha)/(gamma2-1e0)+1e0)*v*p+0.5e0*rho*v*(u^(2e0)+v^(2e0));
        GBC(NY+1,j,5) = beta*rho*v;
        GBC(NY+1,j,6) = gamma1/(gamma1-1e0)*v*p*alpha+1e0/2e0*rho*v*(u^(2e0)+v^(2e0))*beta;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%