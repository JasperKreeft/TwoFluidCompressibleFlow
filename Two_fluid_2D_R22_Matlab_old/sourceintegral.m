function [S]=sourceintegral(NX,NY,WRK3,Wfx,Wfy,gamma1,gamma2)

%--------------------------------------------------%
%                                                  %
%    Source in cell domain                         %
%                                                  %
%--------------------------------------------------%

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

        xi     = WRK3(i,j,1);
        yi     = WRK3(i,j,2);
        ui     = WRK3(i,j,4);
        vi     = WRK3(i,j,5);
        pi     = WRK3(i,j,6);
        betai  = WRK3(i,j,7);
        alphai = WRK3(i,j,8);
        bi     = alphai*(1-alphai)*(gamma1-gamma2)/((1-alphai)*gamma1+alphai*gamma2);

        if (abs(pR-pL)<1D-6 && abs(pA-pB)<1D-6 && abs(uR)<1D-6 && abs(uL)<1D-6 && abs(uA)<1D-6 && abs(uB)<1D-6 && abs(vR)<1D-6 && abs(vL)<1D-6 && abs(vA)<1D-6 && abs(vB)<1D-6)

        S(i,j) = 0D0;

        else

        S(i,j) = 1D0/12D0*(2D0*pi*ui-pi*uR-pi*uA-pR*ui+4D0*pR*uR+3D0*pR*uA-pA*ui+3D0*pA*uR+4D0*pA*uA)*(alphaR-alphai)*(yA-yi)+...
                1D0/12D0*(2D0*pi*vi-pi*vR-pi*vA-pR*vi+4D0*pR*vR+3D0*pR*vA-pA*vi+3D0*pA*vR+4D0*pA*vA)*(alphaA-alphai)*(xR-xi)+...
                1D0/12D0*(2D0*pi*bi-pi*bR-pi*bA-pR*bi+4D0*pR*bR+3D0*pR*bA-pA*bi+3D0*pA*bR+4D0*pA*bA)*(uR-ui)*(yA-yi)+...
                1D0/12D0*(2D0*pi*bi-pi*bR-pi*bA-pR*bi+4D0*pR*bR+3D0*pR*bA-pA*bi+3D0*pA*bR+4D0*pA*bA)*(vA-vi)*(xR-xi)+...
                1D0/12D0*(2D0*(alphai-betai)*ui-(alphai-betai)*uR-(alphai-betai)*uA-(alphaR-betaR)*ui+4D0*(alphaR-betaR)*uR+3D0*(alphaR-betaR)*uA-(alphaA-betaA)*ui+3D0*(alphaA-betaA)*uR+4D0*(alphaA-betaA)*uA)*(pR-pi)*(yA-yi)+...
                1D0/12D0*(2D0*(alphai-betai)*vi-(alphai-betai)*vR-(alphai-betai)*vA-(alphaR-betaR)*vi+4D0*(alphaR-betaR)*vR+3D0*(alphaR-betaR)*vA-(alphaA-betaA)*vi+3D0*(alphaA-betaA)*vR+4D0*(alphaA-betaA)*vA)*(pA-pi)*(xR-xi)+...
               -1D0/12D0*(2D0*pi*ui-pi*uL-pi*uA-pL*ui+4D0*pL*uL+3D0*pL*uA-pA*ui+3D0*pA*uL+4D0*pA*uA)*(alphaL-alphai)*(yA-yi)+...
               -1D0/12D0*(2D0*pi*vi-pi*vL-pi*vA-pL*vi+4D0*pL*vL+3D0*pL*vA-pA*vi+3D0*pA*vL+4D0*pA*vA)*(alphaA-alphai)*(xL-xi)+...
               -1D0/12D0*(2D0*pi*bi-pi*bL-pi*bA-pL*bi+4D0*pL*bL+3D0*pL*bA-pA*bi+3D0*pA*bL+4D0*pA*bA)*(uL-ui)*(yA-yi)+...
               -1D0/12D0*(2D0*pi*bi-pi*bL-pi*bA-pL*bi+4D0*pL*bL+3D0*pL*bA-pA*bi+3D0*pA*bL+4D0*pA*bA)*(vA-vi)*(xL-xi)+...
               -1D0/12D0*(2D0*(alphai-betai)*ui-(alphai-betai)*uL-(alphai-betai)*uA-(alphaL-betaL)*ui+4D0*(alphaL-betaL)*uL+3D0*(alphaL-betaL)*uA-(alphaA-betaA)*ui+3D0*(alphaA-betaA)*uL+4D0*(alphaA-betaA)*uA)*(pL-pi)*(yA-yi)+...
               -1D0/12D0*(2D0*(alphai-betai)*vi-(alphai-betai)*vL-(alphai-betai)*vA-(alphaL-betaL)*vi+4D0*(alphaL-betaL)*vL+3D0*(alphaL-betaL)*vA-(alphaA-betaA)*vi+3D0*(alphaA-betaA)*vL+4D0*(alphaA-betaA)*vA)*(pA-pi)*(xL-xi)+...
                1D0/12D0*(2D0*pi*ui-pi*uL-pi*uB-pL*ui+4D0*pL*uL+3D0*pL*uB-pB*ui+3D0*pB*uL+4D0*pB*uB)*(alphaL-alphai)*(yB-yi)+...
                1D0/12D0*(2D0*pi*vi-pi*vL-pi*vB-pL*vi+4D0*pL*vL+3D0*pL*vB-pB*vi+3D0*pB*vL+4D0*pB*vB)*(alphaB-alphai)*(xL-xi)+...
                1D0/12D0*(2D0*pi*bi-pi*bL-pi*bB-pL*bi+4D0*pL*bL+3D0*pL*bB-pB*bi+3D0*pB*bL+4D0*pB*bB)*(uL-ui)*(yB-yi)+...
                1D0/12D0*(2D0*pi*bi-pi*bL-pi*bB-pL*bi+4D0*pL*bL+3D0*pL*bB-pB*bi+3D0*pB*bL+4D0*pB*bB)*(vB-vi)*(xL-xi)+...
                1D0/12D0*(2D0*(alphai-betai)*ui-(alphai-betai)*uL-(alphai-betai)*uB-(alphaL-betaL)*ui+4D0*(alphaL-betaL)*uL+3D0*(alphaL-betaL)*uB-(alphaB-betaB)*ui+3D0*(alphaB-betaB)*uL+4D0*(alphaB-betaB)*uB)*(pL-pi)*(yB-yi)+...
                1D0/12D0*(2D0*(alphai-betai)*vi-(alphai-betai)*vL-(alphai-betai)*vB-(alphaL-betaL)*vi+4D0*(alphaL-betaL)*vL+3D0*(alphaL-betaL)*vB-(alphaB-betaB)*vi+3D0*(alphaB-betaB)*vL+4D0*(alphaB-betaB)*vB)*(pB-pi)*(xL-xi)+...
               -1D0/12D0*(2D0*pi*ui-pi*uR-pi*uB-pR*ui+4D0*pR*uR+3D0*pR*uB-pB*ui+3D0*pB*uR+4D0*pB*uB)*(alphaR-alphai)*(yB-yi)+...
               -1D0/12D0*(2D0*pi*vi-pi*vR-pi*vB-pR*vi+4D0*pR*vR+3D0*pR*vB-pB*vi+3D0*pB*vR+4D0*pB*vB)*(alphaB-alphai)*(xR-xi)+...
               -1D0/12D0*(2D0*pi*bi-pi*bR-pi*bB-pR*bi+4D0*pR*bR+3D0*pR*bB-pB*bi+3D0*pB*bR+4D0*pB*bB)*(uR-ui)*(yB-yi)+...
               -1D0/12D0*(2D0*pi*bi-pi*bR-pi*bB-pR*bi+4D0*pR*bR+3D0*pR*bB-pB*bi+3D0*pB*bR+4D0*pB*bB)*(vB-vi)*(xR-xi)+...
               -1D0/12D0*(2D0*(alphai-betai)*ui-(alphai-betai)*uR-(alphai-betai)*uB-(alphaR-betaR)*ui+4D0*(alphaR-betaR)*uR+3D0*(alphaR-betaR)*uB-(alphaB-betaB)*ui+3D0*(alphaB-betaB)*uR+4D0*(alphaB-betaB)*uB)*(pR-pi)*(yB-yi)+...
               -1D0/12D0*(2D0*(alphai-betai)*vi-(alphai-betai)*vR-(alphai-betai)*vB-(alphaR-betaR)*vi+4D0*(alphaR-betaR)*vR+3D0*(alphaR-betaR)*vB-(alphaB-betaB)*vi+3D0*(alphaB-betaB)*vR+4D0*(alphaB-betaB)*vB)*(pB-pi)*(xR-xi);

        end

    end
end