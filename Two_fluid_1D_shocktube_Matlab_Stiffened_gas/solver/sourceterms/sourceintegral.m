% Source term in cell domain

function [S]=sourceintegral(N,w,Wf,gamma1,gamma2,pi1,pi2)

S = zeros(1,N);

for i=1:N

    ui     = w(2,i);
    Pi     = w(3,i);
    betai  = w(4,i);
    alphai = w(5,i);
    bi     = alphai*(1-alphai)*(gamma1*(Pi+pi1)-gamma2*(Pi+pi2))/((1-alphai)*gamma1*(Pi+pi1)+alphai*gamma2*(Pi+pi2));

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

    S(i) = (uR*pR/3+uR*Pi/6+ui*pR/6+ui*Pi/3)*(alphai-alphaR)+... %pualpha_x
           (bR*pR/3+bR*Pi/6+bi*pR/6+bi*Pi/3)*(ui-uR)+... %pbu_x
           ((alphaR-betaR)*uR/3+(alphaR-betaR)*ui/6+(alphai-betai)*uR/6+(alphai-betai)*ui/3)*(Pi-pR)+... %(alpha-beta)up_x
           (ui*Pi/3+ui*pL/6+uL*Pi/6+uL*pL/3)*(alphaL-alphai)+... %pualpha_x
           (bi*Pi/3+bi*pL/6+bL*Pi/6+bL*pL/3)*(uL-ui)+... %pbu_x
           ((alphai-betai)*ui/3+(alphai-betai)*uL/6+(alphaL-betaL)*ui/6+(alphaL-betaL)*uL/3)*(pL-Pi); %(alpha-beta)up_x
end