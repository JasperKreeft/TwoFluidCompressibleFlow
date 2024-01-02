function [Fout] = BoundaryConditions(Fin,Wf,gamma,N,BC_L,BC_R)

Fout = Fin;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                          open tube                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Left Boundary
if BC_L==1
Fout(1,1) = Wf(1,1)*Wf(2,1);
Fout(2,1) = Wf(1,1)*Wf(2,1)^2+Wf(3,1);
Fout(3,1) = Wf(3,1)*Wf(2,1)*gamma/(gamma-1)+1/2*Wf(1,1)*Wf(2,1)^3;
end

% Right Boundary
if BC_R==1
Fout(1,N+1) = Wf(1,2*N)*Wf(2,2*N);
Fout(2,N+1) = Wf(1,2*N)*Wf(2,2*N)^2+Wf(3,2*N);
Fout(3,N+1) = Wf(3,2*N)*Wf(2,2*N)*gamma/(gamma-1)+1/2*Wf(1,2*N)*Wf(2,2*N)^3;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                           Solid wall                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Left Boundary
if BC_L==2

% linear Osher
% a = sqrt(gamma*Wf(3,1)/Wf(1,1));
% p = Wf(3,1)-Wf(1,1)*a*Wf(2,1);

a4 = sqrt(gamma*Wf(3,1)/Wf(1,1));
a  = a4-(gamma-1)/2*Wf(2,1);
p  = (Wf(3,1)/(Wf(1,1)^gamma)*(gamma/a^2)^gamma)^(1/(1-gamma));

Fout(1,1) = 0;
Fout(2,1) = p;
Fout(3,1) = 0;
end


%Right Boundary
if BC_R==2

% linear Osher
% a = sqrt(gamma*Wf(3,2*N)/Wf(1,2*N));
% p = Wf(3,2*N)+Wf(1,2*N)*a*Wf(2,2*N);

a1 = sqrt(gamma*Wf(3,2*N)/Wf(1,2*N));
a  = a1+(gamma-1)/2*Wf(2,2*N);
p  = (Wf(3,2*N)/(Wf(1,2*N)^gamma)*(gamma/a^2)^gamma)^(1/(1-gamma));

Fout(1,N+1) = 0;
Fout(2,N+1) = p;
Fout(3,N+1) = 0;
end