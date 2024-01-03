function [Fout] = BoundaryConditions(Fin,Wf,gamma1,gamma2,pi1,pi2,N,BC_L,BC_R)

Fout = Fin;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                          open tube                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Left Boundary
if BC_L==1
Fout(:,1) = [Wf(1,1)*Wf(2,1);
             Wf(1,1)*Wf(2,1)^2+Wf(3,1);
             Wf(5,1)*gamma1*(Wf(3,1)+pi1)*Wf(2,1)/(gamma1-1)+(1-Wf(5,1))*gamma2*(Wf(3,1)+pi2)*Wf(2,1)/(gamma2-1)+1/2*Wf(1,1)*Wf(2,1)^3;
             Wf(1,1)*Wf(2,1)*Wf(4,1);
             Wf(5,1)*gamma1*(Wf(3,1)+pi1)*Wf(2,1)/(gamma1-1)+1/2*Wf(1,1)*Wf(2,1)^3*Wf(4,1) ];
end

% Right Boundary
if BC_R==1
Fout(:,N+1) = [Wf(1,2*N)*Wf(2,2*N);
               Wf(1,2*N)*Wf(2,2*N)^2+Wf(3,2*N);
               Wf(5,2*N)*gamma1*(Wf(3,2*N)+pi1)*Wf(2,2*N)/(gamma1-1)+(1-Wf(5,2*N))*gamma2*(Wf(3,2*N)+pi2)*Wf(2,2*N)/(gamma2-1)+1/2*Wf(1,2*N)*Wf(2,2*N)^3;
               Wf(1,2*N)*Wf(2,2*N)*Wf(4,2*N);
               Wf(5,2*N)*gamma1*(Wf(3,2*N)+pi1)*Wf(2,2*N)/(gamma1-1)+1/2*Wf(1,2*N)*Wf(2,2*N)^3*Wf(4,2*N) ];
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                           Solid wall                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Left Boundary
if BC_L==2

    gamma = 1/(Wf(5,1)/gamma1+(1-Wf(5,1))/gamma2);  %% NIET ECHT TWO-FLUID
    
    a4 = sqrt(gamma*Wf(3,1)/Wf(1,1));
    a  = a4-(gamma-1)/2*Wf(2,1);
    p  = (Wf(3,1)/(Wf(1,1)^gamma)*(gamma/a^2)^gamma)^(1/(1-gamma));

    Fout(1,1) = 0;
    Fout(2,1) = p;
    Fout(3,1) = 0;
    Fout(4,1) = 0;
    Fout(5,1) = 0;
end


%Right Boundary       % NIET KLAAR VOOR STIFFENED GAS LAW %
if BC_R==2

        gamma = 1/(Wf(5,1)/gamma1+(1-Wf(5,1))/gamma2); %% NIET ECHT TWO-FLUID
        
    a1 = sqrt(gamma*Wf(3,2*N)/Wf(1,2*N));
    a  = a1+(gamma-1)/2*Wf(2,2*N);
    p  = (Wf(3,2*N)/(Wf(1,2*N)^gamma)*(gamma/a^2)^gamma)^(1/(1-gamma));

    Fout(1,N+1) = 0;
    Fout(2,N+1) = p;
    Fout(3,N+1) = 0;
    Fout(5,N+1) = 0;
    Fout(4,N+1) = 0;
end