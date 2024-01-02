function [Fout,Gout] = BoundaryConditions(Fin,Gin,Wfx,Wfy,gamma,NX,NY,BC)

Fout = Fin;
Gout = Gin;

BC_L = BC(1);
BC_R = BC(2);
BC_A = BC(3);
BC_B = BC(4);

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
%       L = left                                             %
%       R = right                                            %
%       A = above                                            %
%       B = below                                            %
%                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                      Supersonic outflow                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Left Boundary
if BC_L==1
    Fout(:,1,1) = Wfx(:,1,1).*Wfx(:,1,2);
    Fout(:,1,2) = Wfx(:,1,1).*Wfx(:,1,2).^2+Wfx(:,1,4);
    Fout(:,1,3) = Wfx(:,1,1).*Wfx(:,1,2).*Wfx(:,1,3);
    Fout(:,1,4) = Wfx(:,1,2).*Wfx(:,1,4)*gamma/(gamma-1)+1/2*Wfx(:,1,1).*Wfx(:,1,2).*(Wfx(:,1,2).^2+Wfx(:,1,3).^2);
end

% Right Boundary
if BC_R==1
    Fout(:,NX+1,1) = Wfx(:,end,1).*Wfx(:,end,2);
    Fout(:,NX+1,2) = Wfx(:,end,1).*Wfx(:,end,2).^2+Wfx(:,end,4);
    Fout(:,NX+1,3) = Wfx(:,end,1).*Wfx(:,end,2).*Wfx(:,end,3);
    Fout(:,NX+1,4) = Wfx(:,end,2).*Wfx(:,end,4)*gamma/(gamma-1)+1/2*Wfx(:,end,1).*Wfx(:,end,2).*(Wfx(:,end,2).^2+Wfx(:,end,3).^2);
end

% upper boundary
if BC_A==1
    Gout(1,:,1) = Wfy(1,:,1).*Wfy(1,:,3);
    Gout(1,:,2) = Wfy(1,:,1).*Wfy(1,:,3).*Wfy(1,:,2);
    Gout(1,:,3) = Wfy(1,:,1).*Wfy(1,:,3).^2+Wfy(1,:,4);
    Gout(1,:,4) = Wfy(1,:,3).*Wfy(1,:,4)*gamma/(gamma-1)+1/2*Wfy(1,:,1).*Wfy(1,:,3).*(Wfy(1,:,2).^2+Wfy(1,:,3).^2);
end

% lower Boundary
if BC_B==1
    Gout(NY+1,:,1) = Wfy(end,:,1).*Wfy(end,:,3);
    Gout(NY+1,:,2) = Wfy(end,:,1).*Wfy(end,:,3).*Wfy(end,:,2);
    Gout(NY+1,:,3) = Wfy(end,:,1).*Wfy(end,:,3).^2+Wfy(end,:,4);
    Gout(NY+1,:,4) = Wfy(end,:,3).*Wfy(end,:,4)*gamma/(gamma-1)+1/2*Wfy(end,:,1).*Wfy(end,:,3).*(Wfy(end,:,2).^2+Wfy(end,:,3).^2);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                      Supersonic inflow                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhoB = 1;
uB   = 2.9;
vB   = 0;
pB   = 1/gamma;

%Left Boundary
if BC_L==2
    Fout(:,1,1) = rhoB*uB;
    Fout(:,1,2) = rhoB*uB^2+pB;
    Fout(:,1,3) = rhoB*uB*vB;
    Fout(:,1,4) = uB*pB*gamma/(gamma-1)+1/2*rhoB*uB*(uB^2+vB^2);
end

%Right Boundary
if BC_R==2
    Fout(:,NX+1,1) = rhoB*uB;
    Fout(:,NX+1,2) = rhoB*uB^2+pB;
    Fout(:,NX+1,3) = rhoB*uB*vB;
    Fout(:,NX+1,4) = uB*pB*gamma/(gamma-1)+1/2*rhoB*uB*(uB^2+vB^2);
end

rhoB = 1.74;
uB   = 2.6193;
vB   = 0.5063;
pB   = 1.5282;


%upper Boundary
if BC_A==2
    Gout(1,:,1) = rhoB*vB;
    Gout(1,:,2) = rhoB*vB*uB;
    Gout(1,:,3) = rhoB*vB^2+pB;
    Gout(1,:,4) = vB*pB*gamma/(gamma-1)+1/2*rhoB*vB*(uB^2+vB^2);
end

%lower Boundary
if BC_B==2
    Gout(NY+1,:,1) = rhoB*vB;
    Gout(NY+1,:,2) = rhoB*vB*uB;
    Gout(NY+1,:,3) = rhoB*vB^2+pB;
    Gout(NY+1,:,4) = vB*pB*gamma/(gamma-1)+1/2*rhoB*vB*(uB^2+vB^2);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                      Subsonic outflow                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pB = 1;

%Left Boundary
if BC_L==3
    % Osher
    rho = Wfx(:,1,1).*(pB/Wfx(:,1,4)).^(1/gamma);
    a   = sqrt(gamma*pB./rho);
    a4  = sqrt(gamma*Wfx(:,1,4)./Wfx(:,1,1));
    u   = Wfx(:,1,2)-2/(gamma-1)*(a4-a);
    v   = Wfx(:,1,3);

    Fout(:,1,1) = rho.*u;
    Fout(:,1,2) = rho.*u.^2+pB;
    Fout(:,1,3) = rho.*u.*v;
    Fout(:,1,4) = u.*pB*gamma/(gamma-1)+1/2*rho.*u.*(u.^2+v.^2);
end

%Right Boundary
if BC_R==3
    % Osher
    rho = Wfx(:,end,1).*(pB/Wfx(:,end,4)).^(1/gamma);
    a   = sqrt(gamma*pB./rho);
    a1  = sqrt(gamma*Wfx(:,end,4)./Wfx(:,end,1));
    u   = Wfx(:,end,2)-2/(gamma-1)*(a1-a);
    v   = Wfx(:,end,3);

    Fout(:,NX+1,1) = rho.*u;
    Fout(:,NX+1,2) = rho.*u.^2+pB;
    Fout(:,NX+1,3) = rho.*u.*v;
    Fout(:,NX+1,4) = u.*pB*gamma/(gamma-1)+1/2*rho.*u.*(u.^2+v.^2);
end

%upper Boundary
if BC_A==3
    % Osher
    rho = Wfy(1,:,1).*(pB/Wfy(1,:,4)).^(1/gamma);
    a   = sqrt(gamma*pB./rho);
    a4  = sqrt(gamma*Wfy(1,:,4)./Wfy(1,:,1));
    v   = Wfy(1,:,2)-2/(gamma-1)*(a4-a);
    u   = Wfy(1,:,3);

    Gout(1,:,1) = rho.*v;
    Gout(1,:,2) = rho.*v.*u;
    Gout(1,:,3) = rho.*v.^2+pB;
    Gout(1,:,4) = u.*pB*gamma/(gamma-1)+1/2*rho.*v.*(u.^2+v.^2);
end

%lower Boundary
if BC_B==3
    % Osher
    rho = Wfy(end,:,1).*(pB/Wfy(end,:,4)).^(1/gamma);
    a   = sqrt(gamma*pB./rho);
    a1  = sqrt(gamma*Wfy(end,:,4)./Wfy(end,:,1));
    v   = Wfy(end,:,3)-2/(gamma-1)*(a1-a);
    u   = Wfy(end,:,2);

    Gout(NY+1,:,1) = rho.*v;
    Gout(NY+1,:,2) = rho.*v.*u;
    Gout(NY+1,:,3) = rho.*v.^2+pB;
    Gout(NY+1,:,4) = v.*pB*gamma/(gamma-1)+1/2*rho.*v.*(u.^2+v.^2);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                      Subsonic inflow                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhoB   = 1.0;
uB     = 1.0;

%Left Boundary
if BC_L==4
    % Osher
    a4 = sqrt(gamma*(Wfx(:,1,4)./Wfx(:,1,1)));
    a3 = a4-(gamma-1)/2*(Wfx(:,1,2)-uB);
    p  = (Wfx(:,1,4)./Wfx(:,1,1).^gamma*(gamma/a3.^2).^gamma).^(1/(1-gamma));
    v  = Wfx(:,1,3);

    Fout(:,1,1) = rhoB.*uB;
    Fout(:,1,2) = rhoB.*uB.^2+p;
    Fout(:,1,3) = rhoB.*uB.*v;
    Fout(:,1,4) = uB.*p*gamma/(gamma-1)+1/2*rhoB.*uB.*(uB.^2+v.^2);
end


%Right Boundary
if BC_R==4
    % Osher
    a1 = sqrt(gamma*(Wfx(:,end,4)./Wfx(:,end,1)));
    a2 = a1+(gamma-1)/2*(Wfx(:,end,2)-uB);
    p  = (Wfx(:,end,4)./Wfx(:,end,1).^gamma*(gamma/a2.^2).^gamma).^(1/(1-gamma));
    v  = Wfx(:,end,3);

    Fout(:,NX+1,1) = rhoB.*uB;
    Fout(:,NX+1,2) = rhoB.*uB.^2+p;
    Fout(:,NX+1,3) = rhoB.*uB.*v;
    Fout(:,NX+1,4) = uB.*p*gamma/(gamma-1)+1/2*rhoB.*uB.*(uB.^2+v.^2);
end

rhoB   = 1.0;
vB     = 1.0;

%upper Boundary
if BC_A==4
    % Osher
    a4 = sqrt(gamma*(Wfy(1,:,4)./Wfy(1,:,1)));
    a3 = a4-(gamma-1)/2*(Wfy(1,:,3)-vB);
    p  = (Wfy(1,:,4)./Wfy(1,:,1).^gamma*(gamma/a3.^2).^gamma).^(1/(1-gamma));
    u = Wfy(1,:,2);

    Gout(1,:,1) = rhoB.*vB;
    Gout(1,:,2) = rhoB.*vB.*u;
    Gout(1,:,3) = rhoB.*vB.^2+p;
    Gout(1,:,4) = vB.*p*gamma/(gamma-1)+1/2*rhoB.*vB.*(u.^2+vB.^2);
end


%lower Boundary
if BC_B==4
    % Osher
    a1 = sqrt(gamma*(Wfx(end,:,4)./Wfy(end,:,1)));
    a2 = a1+(gamma-1)/2*(Wfy(end,:,3)-vB);
    p  = (Wfy(end,:,4)./Wfy(end,:,1).^gamma*(gamma/a2.^2).^gamma).^(1/(1-gamma));
    u  = Wfy(end,:,2);
    
    Gout(NY+1,:,1) = rhoB.*vB;
    Gout(NY+1,:,2) = rhoB.*vB.*u;
    Gout(NY+1,:,3) = rhoB.*vB.^2+p;
    Gout(NY+1,:,4) = vB.*p*gamma/(gamma-1)+1/2*rhoB.*vB.*(u.^2+vB.^2);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                           Solid wall                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Left Boundary
if BC_L==5
    % Osher
    a4 = sqrt(gamma*Wfx(:,1,4)./Wfx(:,1,1));
    a  = a4-(gamma-1)/2*Wfx(:,1,2);
    p  = (Wfx(:,1,4)./(Wfx(:,1,1).^gamma).*(gamma./a.^2).^gamma).^(1/(1-gamma));

    Fout(:,1,1) = 0;
    Fout(:,1,2) = p;
    Fout(:,1,3) = 0;
    Fout(:,1,4) = 0;
end


%Right Boundary
if BC_R==5
    % Osher
    a1 = sqrt(gamma*Wfx(:,end,4)./Wfx(:,end,1));
    a  = a1+(gamma-1)/2*Wfx(:,end,2);
    p  = (Wfx(:,end,4)./(Wfx(:,end,1).^gamma).*(gamma./a.^2).^gamma).^(1/(1-gamma));

    Fout(:,NX+1,1) = 0;
    Fout(:,NX+1,2) = p;
    Fout(:,NX+1,3) = 0;
    Fout(:,NX+1,4) = 0;
end

%upper Boundary
if BC_A==5
    % Osher
    a4 = sqrt(gamma*Wfy(1,:,4)./Wfy(1,:,1));
    a  = a4-(gamma-1)/2*Wfy(1,:,3);
    p  = (Wfy(1,:,4)./(Wfy(1,:,1).^gamma).*(gamma./a.^2).^gamma).^(1/(1-gamma));

    Gout(1,:,1) = 0;
    Gout(1,:,2) = 0;
    Gout(1,:,3) = p;
    Gout(1,:,4) = 0;
end


%lower Boundary
if BC_B==5
    % Osher
    a1 = sqrt(gamma*Wfy(end,:,4)./Wfy(end,:,1));
    a  = a1+(gamma-1)/2*Wfy(end,:,3);
    p  = (Wfy(end,:,4)./(Wfy(end,:,1).^gamma).*(gamma./a.^2).^gamma).^(1/(1-gamma));

    Gout(NY+1,:,1) = 0;
    Gout(NY+1,:,2) = 0;
    Gout(NY+1,:,3) = p;
    Gout(NY+1,:,4) = 0;
end