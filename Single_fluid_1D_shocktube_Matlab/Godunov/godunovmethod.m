for k=1:N+2
w(:,k) = [q(1,k) , q(2,k)/q(1,k), (gamma-1)*q(3,k)-1/2*q(2,k)^2/q(1,k)];
end


for k = [2 1 3:N+2] %Spatial loop

if k == 1
% Boundary conditions
% soft boundary condition
if soft == 1
dw(:,1) = zeros(3,1);
else
%Hard boundary condition
dw(1,k) = -dw(1,k+1);
dw(2,k) = dw(2,k+1);
dw(3,k) = -dw(3,k+1);
end

elseif k == N+2
% Boundary conditions
% soft boundary condition
if soft == 1
dw(:,k) = zeros(3,1);

else
%Hard boundary condition
dw(1,k) = -dw(1,k-1);
dw(2,k) = dw(2,k-1);
dw(3,k) = -dw(3,k-1);
end    
    
    
else    % k=2:N+1
%PREDICTOR

dw1 = w(:,k)-w(:,k-1);
dw2 = w(:,k+1)-w(:,k);

if limiter == 1
%Algebraic average
dw(:,k) = (dw1+dw2)/2;

elseif limiter == 2
%double Minmod limiter
if dw1'*dw2 > 0
DML(1,:) = (dw1+dw2)/2;
DML(2,:) = 2*dw1;
DML(3,:) = 2*dw2;
dw(:,k) = minmod(DML)
else
dw(:,k) = zeros(3,1);
end

elseif limiter == 3
%superbee limiter
if dw1'*dw2 > 0
SBL(1,:) = maxmod(dw1,dw2);
SBL(2,:) = minmod(2*dw1,2*dw2);
dw(:,k) = minmod(SBL);
else
dw(:,k) = zeros(3,1);
end

elseif limiter == 4
%koren limiter
if dw1'*dw2 > 0
KL1 = 2*dw1;
KL2 = 1/3*dw1+2/3*dw2;
KL3 = 2*dw2;
kl1mod = sqrt(sum(KL1.^2));
kl2mod = sqrt(sum(KL2.^2));
kl3mod = sqrt(sum(KL3.^2));
if min([kl1mod kl2mod kl3mod]) == kl1mod
    dw(:,k) = KL1;
elseif min([kl1mod kl2mod kl3mod]) == kl2mod
    dw(:,k) = KL2;
elseif min([kl1mod kl2mod kl3mod]) == kl3mod
    dw(:,k) = KL3;
end
else
dw(:,k) = zeros(3,1);
end
end
end

A = [w(2,k), w(1,k), 0; 0, w(2,k), 1/w(1,k); 0 gamma*w(3,k), w(2,k)]; %Ideal gas

ww(:,k) = w(:,k) - dt/(2*dx)*A*dw(:,k);

%CORRECTOR
wwr(:,k) = ww(:,k) - 1/2*dw(:,k);
wwl(:,k) = ww(:,k) + 1/2*dw(:,k);


qr(:,k) = [wwr(1,k), wwr(1,k)*wwr(2,k), 1/(gamma-1)*wwr(3,k)+1/2*wwr(1,k)*wwr(2,k)^2];
ql(:,k) = [wwl(1,k), wwl(1,k)*wwl(2,k), 1/(gamma-1)*wwl(3,k)+1/2*wwl(1,k)*wwl(2,k)^2];
end %end spatial loop