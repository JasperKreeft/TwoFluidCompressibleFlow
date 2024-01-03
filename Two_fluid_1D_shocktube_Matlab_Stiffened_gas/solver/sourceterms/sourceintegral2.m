% Source term in cell domain

function [S]=sourceintegral2(N,w,Wf,gamma1,gamma2,pi1,pi2)

S = zeros(1,N);

% Averaging method
Wf(:,2:2:2*N-2) = (Wf(:,2:2:2*N-2)+Wf(:,3:2:2*N-1))/2;
Wf(:,3:2:2*N-1) = Wf(:,2:2:2*N-2);

for i=1:N

    uR     = Wf(2,2*i-1);
    pR     = Wf(3,2*i-1);
    betaR  = Wf(4,2*i-1);
    alphaR = Wf(5,2*i-1);
    bR     = alphaR*(1-alphaR)*(gamma1*(pR+pi1)-gamma2*(pR+pi2))/((1-alphaR)*gamma1*(pR+pi1)+alphaR*gamma2*(pR+pi2));

    uL     = Wf(2,2*i);
    pL     = Wf(3,2*i);
    betaL  = Wf(4,2*i);
    alphaL = Wf(5,2*i);
    bL     = alphaL*(1-alphaL)*(gamma1*(pL+pi1)-gamma2*(pL+pi2))/((1-alphaL)*gamma1*(pL+pi1)+alphaL*gamma2*(pL+pi2));

    S(i) = (uR*pR/3+uR*pL/6+uL*pR/6+uL*pL/3)*(alphaL-alphaR)+... %pualpha_x
           (bR*pR/3+bR*pL/6+bL*pR/6+bL*pL/3)*(uL-uR)+... %pbu_x
           ((alphaR-betaR)*uR/3+(alphaR-betaR)*uL/6+(alphaL-betaL)*uR/6+(alphaL-betaL)*uL/3)*(pL-pR); %(alpha-beta)up_x

end