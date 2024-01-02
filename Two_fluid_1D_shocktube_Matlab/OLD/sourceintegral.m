function S=sourceintegral(w,Wf)
% Source term in cell domain

global N
global gamma1 gamma2
global pi1 pi2

if pi1~=0 || pi2~=0
    % Averaging method
    Wf(:,2:2:2*N-2) = (Wf(:,2:2:2*N-2)+Wf(:,3:2:2*N-1))/2;
    Wf(:,3:2:2*N-1) = Wf(:,2:2:2*N-2);
end

ui     = w(2,:);
pi     = w(3,:);
betai  = w(4,:);
alphai = w(5,:);
bi     = alphai.*(1-alphai).*(gamma1*(pi+pi1)-gamma2*(pi+pi2))./((1-alphai)*gamma1.*(pi+pi1)+alphai*gamma2.*(pi+pi2));

uR     = Wf(2,2*(1:N)-1);
pR     = Wf(3,2*(1:N)-1);
betaR  = Wf(4,2*(1:N)-1);
alphaR = Wf(5,2*(1:N)-1);
bR     = alphaR.*(1-alphaR).*(gamma1*(pR+pi1)-gamma2*(pR+pi2))./((1-alphaR)*gamma1.*(pR+pi1)+alphaR*gamma2.*(pR+pi2));

uL     = Wf(2,2*(1:N));
pL     = Wf(3,2*(1:N));
betaL  = Wf(4,2*(1:N));
alphaL = Wf(5,2*(1:N));
bL     = alphaL.*(1-alphaL).*(gamma1*(pL+pi1)-gamma2*(pL+pi2))./((1-alphaL)*gamma1.*(pL+pi1)+alphaL*gamma2.*(pL+pi2));

S = (uR.*pR/3+uR.*pi/6+ui.*pR/6+ui.*pi/3).*(alphai-alphaR)+... %pualpha_x
    (bR.*pR/3+bR.*pi/6+bi.*pR/6+bi.*pi/3).*(ui-uR)+... %pbu_x
    ((alphaR-betaR).*uR/3+(alphaR-betaR).*ui/6+(alphai-betai).*uR/6+(alphai-betai).*ui/3).*(pi-pR)+... %(alpha-beta)up_x
    (ui.*pi/3+ui.*pL/6+uL.*pi/6+uL.*pL/3).*(alphaL-alphai)+... %pualpha_x
    (bi.*pi/3+bi.*pL/6+bL.*pi/6+bL.*pL/3).*(uL-ui)+... %pbu_x
    ((alphai-betai).*ui/3+(alphai-betai).*uL/6+(alphaL-betaL).*ui/6+(alphaL-betaL).*uL/3).*(pL-pi); %(alpha-beta)up_x