function [S]=sourceintegral(W,Wfx,Wfy)

%--------------------------------------------------%
%                                                  %
%    Source in cell domain                         %
%                                                  %
%--------------------------------------------------%

global NX NY
global gamma1 gamma2

S = zeros(NY,NX);

for i=1:NY
    for j=1:NX

        xR     = Wfx(i,2*j-1,1);
        uR     = Wfx(i,2*j-1,4);
        vR     = Wfx(i,2*j-1,5);
        pR     = Wfx(i,2*j-1,6);
        betaR  = Wfx(i,2*j-1,7);
        alphaR = Wfx(i,2*j-1,8);
        bR     = alphaR*(1-alphaR)*(gamma1-gamma2)/((1-alphaR)*gamma1+alphaR*gamma2);

        xL     = Wfx(i,2*j,1);
        uL     = Wfx(i,2*j,4);
        vL     = Wfx(i,2*j,5);
        pL     = Wfx(i,2*j,6);
        betaL  = Wfx(i,2*j,7);
        alphaL = Wfx(i,2*j,8);
        bL     = alphaL*(1-alphaL)*(gamma1-gamma2)/((1-alphaL)*gamma1+alphaL*gamma2);

        yA     = Wfy(2*i-1,j,2);
        uA     = Wfy(2*i-1,j,4);
        vA     = Wfy(2*i-1,j,5);
        pA     = Wfy(2*i-1,j,6);
        betaA  = Wfy(2*i-1,j,7);
        alphaA = Wfy(2*i-1,j,8);
        bA     = alphaA*(1-alphaA)*(gamma1-gamma2)/((1-alphaA)*gamma1+alphaA*gamma2);

        yB     = Wfy(2*i,j,2);
        uB     = Wfy(2*i,j,4);
        vB     = Wfy(2*i,j,5);
        pB     = Wfy(2*i,j,6);
        betaB  = Wfy(2*i,j,7);
        alphaB = Wfy(2*i,j,8);
        bB     = alphaB*(1-alphaB)*(gamma1-gamma2)/((1-alphaB)*gamma1+alphaB*gamma2);

        xi     = W(i,j,1);
        yi     = W(i,j,2);
        ui     = W(i,j,4);
        vi     = W(i,j,5);
        pi     = W(i,j,6);
        betai  = W(i,j,7);
        alphai = W(i,j,8);
        bi     = alphai*(1-alphai)*(gamma1-gamma2)/((1-alphai)*gamma1+alphai*gamma2);

        if (abs(pR-pL)+abs(pA-pB)+abs(uR)+abs(uL)+abs(uA)+abs(uB)+abs(vR)+abs(vL)+abs(vA)+abs(vB))>1e-8

        S(i,j) = (2*pi*ui-pi*uR-pi*uA-pR*ui+4*pR*uR+3*pR*uA-pA*ui+3*pA*uR+4*pA*uA)*(alphaR-alphai)*(yA-yi)+...
                (2*pi*vi-pi*vR-pi*vA-pR*vi+4*pR*vR+3*pR*vA-pA*vi+3*pA*vR+4*pA*vA)*(alphaA-alphai)*(xR-xi)+...
                (2*pi*bi-pi*bR-pi*bA-pR*bi+4*pR*bR+3*pR*bA-pA*bi+3*pA*bR+4*pA*bA)*(uR-ui)*(yA-yi)+...
                (2*pi*bi-pi*bR-pi*bA-pR*bi+4*pR*bR+3*pR*bA-pA*bi+3*pA*bR+4*pA*bA)*(vA-vi)*(xR-xi)+...
                (2*(alphai-betai)*ui-(alphai-betai)*uR-(alphai-betai)*uA-(alphaR-betaR)*ui+4*(alphaR-betaR)*uR+3*(alphaR-betaR)*uA-(alphaA-betaA)*ui+3*(alphaA-betaA)*uR+4*(alphaA-betaA)*uA)*(pR-pi)*(yA-yi)+...
                (2*(alphai-betai)*vi-(alphai-betai)*vR-(alphai-betai)*vA-(alphaR-betaR)*vi+4*(alphaR-betaR)*vR+3*(alphaR-betaR)*vA-(alphaA-betaA)*vi+3*(alphaA-betaA)*vR+4*(alphaA-betaA)*vA)*(pA-pi)*(xR-xi)+...
               -(2*pi*ui-pi*uL-pi*uA-pL*ui+4*pL*uL+3*pL*uA-pA*ui+3*pA*uL+4*pA*uA)*(alphaL-alphai)*(yA-yi)+...
               -(2*pi*vi-pi*vL-pi*vA-pL*vi+4*pL*vL+3*pL*vA-pA*vi+3*pA*vL+4*pA*vA)*(alphaA-alphai)*(xL-xi)+...
               -(2*pi*bi-pi*bL-pi*bA-pL*bi+4*pL*bL+3*pL*bA-pA*bi+3*pA*bL+4*pA*bA)*(uL-ui)*(yA-yi)+...
               -(2*pi*bi-pi*bL-pi*bA-pL*bi+4*pL*bL+3*pL*bA-pA*bi+3*pA*bL+4*pA*bA)*(vA-vi)*(xL-xi)+...
               -(2*(alphai-betai)*ui-(alphai-betai)*uL-(alphai-betai)*uA-(alphaL-betaL)*ui+4*(alphaL-betaL)*uL+3*(alphaL-betaL)*uA-(alphaA-betaA)*ui+3*(alphaA-betaA)*uL+4*(alphaA-betaA)*uA)*(pL-pi)*(yA-yi)+...
               -(2*(alphai-betai)*vi-(alphai-betai)*vL-(alphai-betai)*vA-(alphaL-betaL)*vi+4*(alphaL-betaL)*vL+3*(alphaL-betaL)*vA-(alphaA-betaA)*vi+3*(alphaA-betaA)*vL+4*(alphaA-betaA)*vA)*(pA-pi)*(xL-xi)+...
                (2*pi*ui-pi*uL-pi*uB-pL*ui+4*pL*uL+3*pL*uB-pB*ui+3*pB*uL+4*pB*uB)*(alphaL-alphai)*(yB-yi)+...
                (2*pi*vi-pi*vL-pi*vB-pL*vi+4*pL*vL+3*pL*vB-pB*vi+3*pB*vL+4*pB*vB)*(alphaB-alphai)*(xL-xi)+...
                (2*pi*bi-pi*bL-pi*bB-pL*bi+4*pL*bL+3*pL*bB-pB*bi+3*pB*bL+4*pB*bB)*(uL-ui)*(yB-yi)+...
                (2*pi*bi-pi*bL-pi*bB-pL*bi+4*pL*bL+3*pL*bB-pB*bi+3*pB*bL+4*pB*bB)*(vB-vi)*(xL-xi)+...
                (2*(alphai-betai)*ui-(alphai-betai)*uL-(alphai-betai)*uB-(alphaL-betaL)*ui+4*(alphaL-betaL)*uL+3*(alphaL-betaL)*uB-(alphaB-betaB)*ui+3*(alphaB-betaB)*uL+4*(alphaB-betaB)*uB)*(pL-pi)*(yB-yi)+...
                (2*(alphai-betai)*vi-(alphai-betai)*vL-(alphai-betai)*vB-(alphaL-betaL)*vi+4*(alphaL-betaL)*vL+3*(alphaL-betaL)*vB-(alphaB-betaB)*vi+3*(alphaB-betaB)*vL+4*(alphaB-betaB)*vB)*(pB-pi)*(xL-xi)+...
               -(2*pi*ui-pi*uR-pi*uB-pR*ui+4*pR*uR+3*pR*uB-pB*ui+3*pB*uR+4*pB*uB)*(alphaR-alphai)*(yB-yi)+...
               -(2*pi*vi-pi*vR-pi*vB-pR*vi+4*pR*vR+3*pR*vB-pB*vi+3*pB*vR+4*pB*vB)*(alphaB-alphai)*(xR-xi)+...
               -(2*pi*bi-pi*bR-pi*bB-pR*bi+4*pR*bR+3*pR*bB-pB*bi+3*pB*bR+4*pB*bB)*(uR-ui)*(yB-yi)+...
               -(2*pi*bi-pi*bR-pi*bB-pR*bi+4*pR*bR+3*pR*bB-pB*bi+3*pB*bR+4*pB*bB)*(vB-vi)*(xR-xi)+...
               -(2*(alphai-betai)*ui-(alphai-betai)*uR-(alphai-betai)*uB-(alphaR-betaR)*ui+4*(alphaR-betaR)*uR+3*(alphaR-betaR)*uB-(alphaB-betaB)*ui+3*(alphaB-betaB)*uR+4*(alphaB-betaB)*uB)*(pR-pi)*(yB-yi)+...
               -(2*(alphai-betai)*vi-(alphai-betai)*vR-(alphai-betai)*vB-(alphaR-betaR)*vi+4*(alphaR-betaR)*vR+3*(alphaR-betaR)*vB-(alphaB-betaB)*vi+3*(alphaB-betaB)*vR+4*(alphaB-betaB)*vB)*(pB-pi)*(xR-xi);

        end

    end
end

S = S/12;



