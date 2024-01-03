function [Fout,Gout] = BoundaryConditions(Fin,Gin,Wfx,Wfy,gamma,NX,NY,BC)

BC = [2 1 2 5];

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

% Right Boundary
Fout(:,NX+1,1) = Wfx(:,end,1).*Wfx(:,end,2);
Fout(:,NX+1,2) = Wfx(:,end,1).*Wfx(:,end,2).^2+Wfx(:,end,4);
Fout(:,NX+1,3) = Wfx(:,end,1).*Wfx(:,end,2).*Wfx(:,end,3);
Fout(:,NX+1,4) = Wfx(:,end,2).*Wfx(:,end,4)*gamma/(gamma-1)+1/2*Wfx(:,end,1).*Wfx(:,end,2).*(Wfx(:,end,2).^2+Wfx(:,end,3).^2);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                      Supersonic inflow                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhoB = 1;
uB   = 2.9;
vB   = 0;
pB   = 1/gamma;

%Left Boundary
Fout(:,1,1) = rhoB*uB;
Fout(:,1,2) = rhoB*uB^2+pB;
Fout(:,1,3) = rhoB*uB*vB;
Fout(:,1,4) = uB*pB*gamma/(gamma-1)+1/2*rhoB*uB*(uB^2+vB^2);

%upper Boundary
Gout(1,:,1) = 0;
Gout(1,:,2) = 0;
Gout(1,:,3) = pB;
Gout(1,:,4) = 0;

% rhoB   = 1.0;
% vB     = -0.5;
% uB     = sqrt(2.9^2-vB^2);
% pB   = 1/gamma;
% 
% rr = ceil(1/2*NX):NX;
% 
% Gout(NY+1,rr,1) = rhoB*vB;
% Gout(NY+1,rr,2) = rhoB*vB*uB;
% Gout(NY+1,rr,3) = rhoB*vB^2+pB;
% Gout(NY+1,rr,4) = vB*pB*gamma/(gamma-1)+1/2*rhoB*vB*(uB^2+vB^2);
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                      Subsonic inflow                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhoB =  1.0;
vB   = -0.5;

rr = ceil(1/2*NX):NX;

% Osher
a1 = sqrt(gamma*(Wfy(end,rr,4)./Wfy(end,rr,1)));
a2 = a1+(gamma-1)/2*(Wfy(end,rr,3)-vB);
p  = (Wfy(end,rr,4)./Wfy(end,rr,1).^gamma.*(gamma./a2.^2).^gamma).^(1/(1-gamma));
u  = Wfy(end,rr,2);

Gout(NY+1,rr,1) = rhoB.*vB;
Gout(NY+1,rr,2) = rhoB.*vB.*u;
Gout(NY+1,rr,3) = rhoB.*vB.^2+p;
Gout(NY+1,rr,4) = vB.*p*gamma/(gamma-1)+1/2*rhoB.*vB.*(u.^2+vB.^2);
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                           Solid wall                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%lower Boundary

ss = 1:floor(1/2*NX);

% Osher
a1 = sqrt(gamma*Wfy(end,ss,4)./Wfy(end,ss,1));
a  = a1+(gamma-1)/2*Wfy(end,ss,3);
p  = (Wfy(end,ss,4)./(Wfy(end,ss,1).^gamma).*(gamma./a.^2).^gamma).^(1/(1-gamma));

Gout(NY+1,ss,1) = 0;
Gout(NY+1,ss,2) = 0;
Gout(NY+1,ss,3) = p;
Gout(NY+1,ss,4) = 0;