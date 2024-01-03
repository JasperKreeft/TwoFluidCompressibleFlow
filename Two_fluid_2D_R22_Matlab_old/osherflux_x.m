function [F,sL,sR]=osherflux_x(NX,NY,FBC,Wfx,gamma1,gamma2)

%---------------------------------------%
%                                       %
% Osher's flux for two-fluid flow       %
%                                       %
% author: Jasper Kreeft                 %
%                                       %
%---------------------------------------%

%--------------------------------------%
% osher scheme                         %
%--------------------------------------%

% INTEGER, INTENT(IN) :: NX,NY
% DOUBLE PRECISION, INTENT(IN) :: Wfx(NY,2*NX,8)
% DOUBLE PRECISION, INTENT(IN) :: gamma1,gamma2
% DOUBLE PRECISION, INTENT(IN) :: FBC(NY,NX+1,6)
% DOUBLE PRECISION, INTENT(OUT) :: F(NY,NX+1,6)
% DOUBLE PRECISION, INTENT(OUT) :: sR(NY,NX),sL(NY,NX)
% 
% DOUBLE PRECISION :: rhol,ul,vl,pl,betal,alphal,al,gammal,rhor,ur,vr,pr,betar,alphar,ar,gammar
% DOUBLE PRECISION :: a2,a3,rho2,rho3,pf,uf,vf2,vf3,beta2,beta3,alpha2,alpha3
% DOUBLE PRECISION :: rhosl,usl,vsl,psl,betasl,alphasl,asl,rhosr,usr,vsr,psr,betasr,alphasr,asr
% DOUBLE PRECISION :: sint12,sint1sl,sintsl2,sintCD,sint34,sint3sr,sintsr4
% DOUBLE PRECISION :: lambda0,lambda1,lambda2,lambda3,lambda4,lambda
% DOUBLE PRECISION :: Fl(6),F2(6),F3(6),Fr(6),Fsl(6),Fsr(6)
% DOUBLE PRECISION, PARAMETER :: eps=1D-15
% DOUBLE PRECISION :: diffw,maxdiffw
% INTEGER :: i,j,k

F = FBC;

sL = zeros(NY,NX);
sR = zeros(NY,NX);

for i=1:NY
for j=2:NX

maxdiffw = 0D0;
for k=3:6
diffw = abs(Wfx(i,2*(j-1),k)-Wfx(i,2*j-1,k));
maxdiffw = max(maxdiffw,diffw);
end

if maxdiffw<1e-6

F(i,j,1) = Wfx(i,2*(j-1),3)*Wfx(i,2*(j-1),4);
F(i,j,2) = Wfx(i,2*(j-1),3)*Wfx(i,2*(j-1),4)^(2D0)+Wfx(i,2*(j-1),6);
F(i,j,3) = Wfx(i,2*(j-1),3)*Wfx(i,2*(j-1),4)*Wfx(i,2*(j-1),5);
F(i,j,4) = (Wfx(i,2*(j-1),8)/(gamma1-1D0)+(1D0-Wfx(i,2*(j-1),8))/(gamma2-1D0)+1D0)*Wfx(i,2*(j-1),4)*Wfx(i,2*(j-1),6)+0.5D0*Wfx(i,2*(j-1),3)*Wfx(i,2*(j-1),4)*(Wfx(i,2*(j-1),4)^(2D0)+Wfx(i,2*(j-1),5)^(2D0));
F(i,j,5) = Wfx(i,2*(j-1),7)*Wfx(i,2*(j-1),3)*Wfx(i,2*(j-1),4);
F(i,j,6) = gamma1/(gamma1-1D0)*Wfx(i,2*(j-1),4)*Wfx(i,2*(j-1),6)*Wfx(i,2*(j-1),8)+1D0/2D0*Wfx(i,2*(j-1),3)*Wfx(i,2*(j-1),4)*(Wfx(i,2*(j-1),4)^(2D0)+Wfx(i,2*(j-1),5)^(2D0))*Wfx(i,2*(j-1),7);

sL(i,j)  = 0D0;
sR(i,j)  = 0D0;

elseif (abs(Wfx(i,2*(j-1),6)-Wfx(i,2*j-1,6))<1D-6 && abs(Wfx(i,2*(j-1),4))<1D-6 && abs(Wfx(i,2*j-1,4))<1D-6)

F(i,j,1) = 0D0;
F(i,j,2) = Wfx(i,2*(j-1),6);
F(i,j,3) = 0D0;
F(i,j,4) = 0D0;
F(i,j,5) = 0D0;
F(i,j,6) = 0D0;

sL(i,j)  = 0D0;
sR(i,j)  = 0D0;

else

rhol   = Wfx(i,2*(j-1),3);
ul     = Wfx(i,2*(j-1),4);
vl     = Wfx(i,2*(j-1),5);
pl     = Wfx(i,2*(j-1),6);
betal  = Wfx(i,2*(j-1),7);
alphal = Wfx(i,2*(j-1),8);
gammal = 1D0/(alphal/(gamma1-1D0)+(1D0-alphal)/(gamma2-1D0))+1D0;
al     = sqrt(gammal*pl/rhol);

rhor   = Wfx(i,2*j-1,3);
ur     = Wfx(i,2*j-1,4);
vr     = Wfx(i,2*j-1,5);
pr     = Wfx(i,2*j-1,6);
betar  = Wfx(i,2*j-1,7);
alphar = Wfx(i,2*j-1,8);
gammar = 1D0/(alphar/(gamma1-1D0)+(1D0-alphar)/(gamma2-1D0))+1D0;
ar     = sqrt(gammar*pr/rhor);

% Final state values
[rho2,rho3,uf,vf2,vf3,pf,beta2,beta3,alpha2,alpha3]=oshersolver_x(rhol,ul,vl,pl,betal,alphal,al,rhor,ur,vr,pr,betar,alphar,ar,gamma1,gamma2);

% speed of sound
a2 = sqrt((1D0/(alpha2/(gamma1-1D0)+(1D0-alpha2)/(gamma2-1D0))+1D0)*pf/rho2);
a3 = sqrt((1D0/(alpha3/(gamma1-1D0)+(1D0-alpha3)/(gamma2-1D0))+1D0)*pf/rho3);

% Sonic point values
[rhosl,usl,vsl,psl,betasl,alphasl]=osher_sonic_L_x(ur,pr,ar,gammar,rhol,ul,vl,pl,betal,alphal,al,gammal,gamma1,gamma2);
[rhosr,usr,vsr,psr,betasr,alphasr]=osher_sonic_R_x(ul,pl,al,gammal,rhor,ur,vr,pr,betar,alphar,ar,gammar,gamma1,gamma2);

Fl(1) = rhol*ul;
Fl(2) = rhol*ul^(2D0)+pl;
Fl(3) = rhol*ul*vl;
Fl(4) = (alphal/(gamma1-1D0)+(1D0-alphal)/(gamma2-1D0)+1D0)*ul*pl+1D0/2D0*rhol*ul*(ul^(2D0)+vl^(2D0));
Fl(5) = betal*rhol*ul;
Fl(6) = gamma1/(gamma1-1D0)*ul*pl*alphal+1D0/2D0*rhol*ul*(ul^(2D0)+vl^(2D0))*betal;

F2(1) = rho2*uf;
F2(2) = rho2*uf^(2D0)+pf;
F2(3) = rho2*uf*vf2;
F2(4) = (alpha2/(gamma1-1D0)+(1D0-alpha2)/(gamma2-1D0)+1D0)*uf*pf+1D0/2D0*rho2*uf*(uf^(2D0)+vf2^(2D0));
F2(5) = beta2*rho2*uf;
F2(6) = gamma1/(gamma1-1D0)*uf*pf*alpha2+1D0/2D0*rho2*uf*(uf^(2D0)+vf2^(2D0))*beta2;

F3(1) = rho3*uf;
F3(2) = rho3*uf^(2D0)+pf;
F3(3) = rho3*uf*vf3;
F3(4) = (alpha3/(gamma1-1D0)+(1D0-alpha3)/(gamma2-1D0)+1D0)*uf*pf+1D0/2D0*rho3*uf*(uf^(2D0)+vf3^(2D0));
F3(5) = beta3*rho3*uf;
F3(6) = gamma1/(gamma1-1D0)*uf*pf*alpha3+1D0/2D0*rho3*uf*(uf^(2D0)+vf3^(2D0))*beta3;

Fr(1) = rhor*ur;
Fr(2) = rhor*ur^(2D0)+pr;
Fr(3) = rhor*ur*vr;
Fr(4) = (alphar/(gamma1-1D0)+(1D0-alphar)/(gamma2-1D0)+1D0)*ur*pr+1D0/2D0*rhor*ur*(ur^(2D0)+vr^(2D0));
Fr(5) = betar*rhor*ur;
Fr(6) = gamma1/(gamma1-1D0)*ur*pr*alphar+1D0/2D0*rhor*ur*(ur^(2D0)+vr^(2D0))*betar;

Fsl(1) = rhosl*usl;
Fsl(2) = rhosl*usl^(2D0)+psl;
Fsl(3) = rhosl*usl*vsl;
Fsl(4) = (alphasl/(gamma1-1D0)+(1D0-alphasl)/(gamma2-1D0)+1D0)*usl*psl+1D0/2D0*rhosl*usl*(usl^(2D0)+vsl^(2D0));
Fsl(5) = betasl*rhosl*usl;
Fsl(6) = gamma1/(gamma1-1D0)*usl*psl*alphasl+1D0/2D0*rhosl*usl*(usl^(2D0)+vsl^(2D0))*betasl;

Fsr(1) = rhosr*usr;
Fsr(2) = rhosr*usr^(2D0)+psr;
Fsr(3) = rhosr*usr*vsr;
Fsr(4) = (alphasr/(gamma1-1D0)+(1D0-alphasr)/(gamma2-1D0)+1D0)*usr*psr+1D0/2D0*rhosr*usr*(usr^(2D0)+vsr^(2D0));
Fsr(5) = betasr*rhosr*usr;
Fsr(6) = gamma1/(gamma1-1D0)*usr*psr*alphasr+1D0/2D0*rhosr*usr*(usr^(2D0)+vsr^(2D0))*betasr;

lambda0 = uf;
lambda1 = ul-al;
lambda2 = uf-a2;
lambda3 = uf+a3;
lambda4 = ur+ar;
% Lambda  = max(abs(lambda0),abs(lambda1),abs(lambda2),abs(lambda3),abs(lambda4),abs(Lambda)))

if (lambda0>=0D0 && lambda2>=0D0)
    if (lambda1>=0D0 && lambda4>=0D0)
        F(i,j,:) = Fl;
        CALL sub_Sint12_x(rhol,ul,pl,betal,alphal,uf,gamma1,gamma2,sint12)
        CALL sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
        CALL sub_Sint34_x(rho3,uf,pf,beta3,alpha3,ur,gamma1,gamma2,sint34)
        sL(i,j) = 0D0;
        sR(i,j) = sint12+sintCD+sint34;
    elseif (lambda1>=0D0 && lambda4<=0D0)
        F(i,j,:) = Fl+Fr-Fsr;
        CALL sub_Sintsr4_x(rhosr,usr,psr,betasr,alphasr,ur,gamma1,gamma2,sintsr4)
        CALL sub_Sint12_x(rhol,ul,pl,betal,alphal,uf,gamma1,gamma2,sint12)
        CALL sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
        CALL sub_Sint3sr_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,sint3sr)
        sL(i,j) = sintsr4;
        sR(i,j) = sint12+sintCD+sint3sr;
    elseif (lambda1<=0D0 && lambda4>=0D0)
        F(i,j,:) = Fsl;
        CALL sub_Sint1sl_x(rhol,ul,pl,betal,alphal,usl,gamma1,gamma2,sint1sl)
        CALL sub_Sintsl2_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,sintsl2)
        CALL sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
        CALL sub_Sint34_x(rho3,uf,pf,beta3,alpha3,ur,gamma1,gamma2,sint34)
        sL(i,j) = sint1sl;
        sR(i,j) = sintsl2+sintCD+sint34;
    elseif (lambda1<=0D0 && lambda4<=0D0)
        F(i,j,:) = Fsl-Fsr+Fr;
        CALL sub_Sint1sl_x(rhol,ul,pl,betal,alphal,usl,gamma1,gamma2,sint1sl)
        CALL sub_Sintsr4_x(rhosr,usr,psr,betasr,alphasr,ur,gamma1,gamma2,sintsr4)
        CALL sub_Sintsl2_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,sintsl2)
        CALL sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
        CALL sub_Sint3sr_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,sint3sr)
        sL(i,j) = sint1sl+sintsr4;
        sR(i,j) = sintsl2+sintCD+sint3sr;
    end

elseif  (lambda0>=0D0 && lambda2<=0D0)
    if (lambda1>=0D0 && lambda4>=0D0)
        F(i,j,:) = Fl-Fsl+F2;
        CALL sub_Sintsl2_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,sintsl2)
        CALL sub_Sint1sl_x(rhol,ul,pl,betal,alphal,usl,gamma1,gamma2,sint1sl)
        CALL sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
        CALL sub_Sint34_x(rho3,uf,pf,beta3,alpha3,ur,gamma1,gamma2,sint34)
        sL(i,j) = sintsl2;
        sR(i,j) = sint1sl+sintCD+sint34;
    elseif (lambda1>=0D0 && lambda4<=0D0)
        F(i,j,:) = Fl-Fsl+F2-Fsr+Fr;
        CALL sub_Sintsl2_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,sintsl2)
        CALL sub_Sintsr4_x(rhosr,usr,psr,betasr,alphasr,ur,gamma1,gamma2,sintsr4)
        CALL sub_Sint1sl_x(rhol,ul,pl,betal,alphal,usl,gamma1,gamma2,sint1sl)
        CALL sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
        CALL sub_Sint3sr_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,sint3sr)
        sL(i,j) = sintsl2+sintsr4;
        sR(i,j) = sint1sl+sintCD+sint3sr;
    elseif (lambda1<=0D0 && lambda4>=0D0)
        F(i,j,:) = F2;
        CALL sub_Sint12_x(rhol,ul,pl,betal,alphal,uf,gamma1,gamma2,sint12)
        CALL sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
        CALL sub_Sint34_x(rho3,uf,pf,beta3,alpha3,ur,gamma1,gamma2,sint34)
        sL(i,j) = sint12;
        sR(i,j) = sintCD+sint34;
    elseif (lambda1<=0D0 && lambda4<=0D0)
        F(i,j,:) = Fr+F2-Fsr;
        CALL sub_Sint12_x(rhol,ul,pl,betal,alphal,uf,gamma1,gamma2,sint12)
        CALL sub_Sintsr4_x(rhosr,usr,psr,betasr,alphasr,ur,gamma1,gamma2,sintsr4)
        CALL sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
        CALL sub_Sint3sr_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,sint3sr)
        sL(i,j) = sint12+sintsr4;
        sR(i,j) = sintCD+sint3sr;
    end

elseif  (lambda0<=0D0 && lambda3>=0D0)
    if (lambda1>=0D0 && lambda4>=0D0)
        F(i,j,:) = Fl-Fsl+F3;
        CALL sub_Sintsl2_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,sintsl2)
        CALL sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
        CALL sub_Sint1sl_x(rhol,ul,pl,betal,alphal,usl,gamma1,gamma2,sint1sl)
        CALL sub_Sint34_x(rho3,uf,pf,beta3,alpha3,ur,gamma1,gamma2,sint34)
        sL(i,j) = sintsl2+sintCD;
        sR(i,j) = sint1sl+sint34;
    elseif (lambda1>=0D0 && lambda4<=0D0)
        F(i,j,:) = Fl-Fsl+F3-Fsr+Fr;
        CALL sub_Sintsl2_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,sintsl2)
        CALL sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
        CALL sub_Sintsr4_x(rhosr,usr,psr,betasr,alphasr,ur,gamma1,gamma2,sintsr4)
        CALL sub_Sint1sl_x(rhol,ul,pl,betal,alphal,usl,gamma1,gamma2,sint1sl)
        CALL sub_Sint3sr_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,sint3sr)
        sL(i,j) = sintsl2+sintCD+sintsr4;
        sR(i,j) = sint1sl+sint3sr;
    elseif (lambda1<=0D0 && lambda4>=0D0)
        F(i,j,:) = F3;
        CALL sub_Sint12_x(rhol,ul,pl,betal,alphal,uf,gamma1,gamma2,sint12)
        CALL sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
        CALL sub_Sint34_x(rho3,uf,pf,beta3,alpha3,ur,gamma1,gamma2,sint34)
        sL(i,j) = sint12+sintCD;
        sR(i,j) = sint34;
    elseif (lambda1<=0D0 && lambda4<=0D0)
        F(i,j,:) = Fr+F3-Fsr;
        CALL sub_Sint12_x(rhol,ul,pl,betal,alphal,uf,gamma1,gamma2,sint12)
        CALL sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
        CALL sub_Sintsr4_x(rhosr,usr,psr,betasr,alphasr,ur,gamma1,gamma2,sintsr4)
        CALL sub_Sint3sr_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,sint3sr)
        sL(i,j) = sint12+sintCD+sintsr4;
        sR(i,j) = sint3sr;
    end

elseif  (lambda0<=0D0 && lambda3<=0D0)
    if (lambda1>=0D0 && lambda4>=0D0)
        F(i,j,:) = Fl-Fsl+Fsr;
        CALL sub_Sintsl2_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,sintsl2)
        CALL sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
        CALL sub_Sint3sr_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,sint3sr)
        CALL sub_Sint1sl_x(rhol,ul,pl,betal,alphal,usl,gamma1,gamma2,sint1sl)
        CALL sub_Sintsr4_x(rhosr,usr,psr,betasr,alphasr,ur,gamma1,gamma2,sintsr4)
        sL(i,j) = sintsl2+sintCD+sint3sr;
        sR(i,j) = sint1sl+sintsr4;
    elseif (lambda1>=0D0 && lambda4<=0D0)
        F(i,j,:) = Fl-Fsl+Fr;
        CALL sub_Sintsl2_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,sintsl2)
        CALL sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
        CALL sub_Sint34_x(rho3,uf,pf,beta3,alpha3,ur,gamma1,gamma2,sint34)
        CALL sub_Sint1sl_x(rhol,ul,pl,betal,alphal,usl,gamma1,gamma2,sint1sl)
        sL(i,j) = sintsl2+sintCD+sint34;
        sR(i,j) = sint1sl;
    elseif (lambda1<=0D0 && lambda4>=0D0)
        F(i,j,:) = Fsr;
        CALL sub_Sint12_x(rhol,ul,pl,betal,alphal,uf,gamma1,gamma2,sint12)
        CALL sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
        CALL sub_Sint3sr_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,sint3sr)
        CALL sub_Sintsr4_x(rhosr,usr,psr,betasr,alphasr,ur,gamma1,gamma2,sintsr4)
        sL(i,j) = sint12+sintCD+sint3sr;
        sR(i,j) = sintsr4;
    elseif (lambda1<=0D0 && lambda4<=0D0)
        F(i,j,:) = Fr;
        CALL sub_Sint12_x(rhol,ul,pl,betal,alphal,uf,gamma1,gamma2,sint12)
        CALL sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
        CALL sub_Sint34_x(rho3,uf,pf,beta3,alpha3,ur,gamma1,gamma2,sint34)
        sL(i,j) = sint12+sintCD+sint34;
        sR(i,j) = 0D0;
    end

end %flux choice

end %if L==R

end %for j loop
end %for i loop


end % osherflux_x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rho2,rho3,uf,vf2,vf3,pf,beta2,beta3,alpha2,alpha3]=oshersolver_x(rhol,ul,vl,pl,betal,alphal,al,rhor,ur,vr,pr,betar,alphar,ar,gamma1,gamma2)
%-----------------------------------%
%                                   %
% flux solver                       %
%                                   %
%-----------------------------------%

% forUBLE PRECISION, INTENT(IN)  :: rhol,ul,vl,pl,betal,alphal,al,rhor,ur,vr,pr,betar,alphar,ar,gamma1,gamma2
% forUBLE PRECISION, INTENT(OUT) :: rho2,rho3,uf,vf2,vf3,pf,beta2,beta3,alpha2,alpha3
% INTEGER, PARAMETER :: NRK4 = 300
% INTEGER :: m,k
% forUBLE PRECISION :: diffp,dp,dpold
% forUBLE PRECISION :: rho,u,v,p,beta,alpha,uguess,dul,dur,gamma
% forUBLE PRECISION :: pf2,pf3
% forUBLE PRECISION :: ufnew,ufold,pfold
% forUBLE PRECISION :: rho_k1,rho_k2,rho_k3,rho_k4,p_k1,p_k2,p_k3,p_k4,alpha_k1,alpha_k2,alpha_k3,alpha_k4

%-----------------------------------------------------------------------%
%                       Single-fluid initial guess                      %
%-----------------------------------------------------------------------%
gamma  = ((1D0/(alphal/(gamma1-1D0)+(1D0-alphal)/(gamma2-1D0))+1D0)+(1D0/(alphar/(gamma1-1D0)+(1D0-alphar)/(gamma2-1D0))+1D0))/2D0;
uguess = (ul/al*pl^((gamma-1D0)/(2D0*gamma))+ur/ar*pr^((gamma-1D0)/(2D0*gamma))+2D0/(gamma-1D0)*(pl^((gamma-1D0)/(2D0*gamma))-pr^((gamma-1D0)/(2D0*gamma))))/(pl^((gamma-1D0)/(2D0*gamma))/al+pr^((gamma-1D0)/(2D0*gamma))/ar);

%-----------------------------------------------------------------------%
%                       Runga-Kutta fourth order                        %
%-----------------------------------------------------------------------%
diffp = 1D0;
m = 0;

while diffp>1D-10
m=m+1;

%-----------------------------------------------------------------------%
%                               LEFT WAVE                               %
%-----------------------------------------------------------------------%
if m==1
uf = uguess;
elseif m==2
uf = 0.8D0*uguess;
end

rho   = rhol;
v     = vl;
p     = pl;
beta  = betal;
alpha = alphal;
dul   = (uf - ul)/NRK4;

for k=1:NRK4
    CALL RK4_k1_L_x(rho,p,alpha,gamma1,gamma2,dul,p_k1,rho_k1,alpha_k1)
    CALL RK4_k2_L_x(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dul,p_k2,rho_k2,alpha_k2)
    CALL RK4_k3_L_x(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dul,p_k3,rho_k3,alpha_k3)
    CALL RK4_k4_L_x(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dul,p_k4,rho_k4,alpha_k4)
    rho   = rho+(rho_k1+2D0*rho_k2+2D0*rho_k3+rho_k4)/6D0;
    p     = p+(p_k1+2D0*p_k2+2D0*p_k3+p_k4)/6D0;
    alpha = alpha+(alpha_k1+2D0*alpha_k2+2D0*alpha_k3+alpha_k4)/6D0;
end

rho2   = rho;
vf2    = v;
pf2    = p;
beta2  = beta;
alpha2 = alpha;

%-----------------------------------------------------------------------%
%                              RIGHT WAVE                               %
%-----------------------------------------------------------------------%
rho   = rhor;
v     = vr;
p     = pr;
beta  = betar;
alpha = alphar;
dur   = (uf - ur)/NRK4;

for k=1:NRK4
    CALL RK4_k1_R_x(rho,p,alpha,gamma1,gamma2,dur,p_k1,rho_k1,alpha_k1)
    CALL RK4_k2_R_x(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dur,p_k2,rho_k2,alpha_k2)
    CALL RK4_k3_R_x(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dur,p_k3,rho_k3,alpha_k3)
    CALL RK4_k4_R_x(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dur,p_k4,rho_k4,alpha_k4)
    rho   = rho+(rho_k1+2D0*rho_k2+2D0*rho_k3+rho_k4)/6D0;
    p     = p+(p_k1+2D0*p_k2+2D0*p_k3+p_k4)/6D0;
    alpha = alpha+(alpha_k1+2D0*alpha_k2+2D0*alpha_k3+alpha_k4)/6D0;
end

rho3   = rho;
vf3    = v;
pf3    = p;
beta3  = beta;
alpha3 = alpha;

%----------------------------------------------------

dp    = pf3 - pf2;
diffp = abs(dp);

if m==1
ufold = uf;
dpold = dp;
else
ufnew = uf - dp*(uf - ufold)/(dp - dpold);
ufold = uf;
uf    = ufnew;
dpold = dp;
end

end %while loop

pf = pf2;

end % oshersolver_x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rhosl,usl,vsl,psl,betasl,alphasl]=osher_sonic_L_x(ur,pr,ar,gammar,rhol,ul,vl,pl,betal,alphal,al,gammal,gamma1,gamma2,rhosl,usl,vsl,psl,betasl,alphasl)

%------------------------------------%
%                                    %
% Solver for left sonic point        %
%                                    %
%------------------------------------%

% forUBLE PRECISION, INTENT(IN)  :: ur,pr,ar,gammar,rhol,ul,vl,pl,betal,alphal,al,gammal,gamma1,gamma2
% forUBLE PRECISION, INTENT(out) :: rhosl,usl,vsl,psl,betasl,alphasl
% INTEGER, PARAMETER :: NRK4   = 300
% INTEGER :: m,k
% forUBLE PRECISION :: diffu,dua
% forUBLE PRECISION :: uslnew,uslold,asl,duaold
% forUBLE PRECISION :: rho,u,v,p,beta,alpha,dul
% forUBLE PRECISION :: rho_k1,rho_k2,rho_k3,rho_k4,p_k1,p_k2,p_k3,p_k4,alpha_k1,alpha_k2,alpha_k3,alpha_k4
% forUBLE PRECISION :: gamma,uguess
% forUBLE PRECISION, PARAMETER :: eps=1D-15


%Runga-Kutta Solver
gamma  = (gammar+gammal)/2D0;
uguess = (ul/al*pl^((gamma-1D0)/(2D0*gamma))+ur/ar*pr^((gamma-1D0)/(2D0*gamma))+2D0/(gamma-1D0)*(pl^((gamma-1D0)/(2D0*gamma))-pl^((gamma-1D0)/(2D0*gamma))))/(pr^((gamma-1D0)/(2D0*gamma))/ar+pr^((gamma-1D0)/(2D0*gamma))/ar);

%-----------------------------------------------------------------------%
%                         runga-kutta fourth-order                      %
%-----------------------------------------------------------------------%
diffu = 1D0;
m = 0;

while diffu>1D-10
m=m+1;

%-----------------------------------------------------------------------%
%                               LEFT WAVE                               %
%-----------------------------------------------------------------------%
if m==1
usl = uguess;
elseif m==2
usl = 0.8D0*uguess;
end

rho   = rhol;
v     = vl;
p     = pl;
beta  = betal;
alpha = alphal;
dul   = (usl - ul)/NRK4;

for k=1:NRK4
    CALL RK4_k1_L_x(rho,p,alpha,gamma1,gamma2,dul,p_k1,rho_k1,alpha_k1)
    CALL RK4_k2_L_x(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dul,p_k2,rho_k2,alpha_k2)
    CALL RK4_k3_L_x(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dul,p_k3,rho_k3,alpha_k3)
    CALL RK4_k4_L_x(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dul,p_k4,rho_k4,alpha_k4)
    p = p + (p_k1 + 2D0*p_k2 + 2D0*p_k3 + p_k4)/6D0;
    rho = rho + (rho_k1 + 2D0*rho_k2 + 2D0*rho_k3 + rho_k4)/6D0;
    alpha = alpha + (alpha_k1 + 2D0*alpha_k2 + 2D0*alpha_k3 + alpha_k4)/6D0;
end

rhosl   = rho;
vsl     = v;
psl     = p;
alphasl = alpha;
betasl  = beta;
gamma   = 1D0/(alpha/(gamma1-1D0)+(1D0-alpha)/(gamma2-1D0))+1D0;
asl     = sqrt(gamma*psl/rhosl);

%------------------------------------------------------------------------

dua     = usl - asl;
diffu   = abs(dua);

if m==1
uslold = usl;
duaold = dua;
else
uslnew = usl - dua*(usl-uslold+eps)/(dua-duaold+eps);
uslold = usl;
usl    = uslnew;
duaold = dua;
end

end %while loop


end % osher_sonic_L_x



function [rhosr,usr,vsr,psr,betasr,alphasr]=osher_sonic_R_x(ul,pl,al,gammal,rhor,ur,vr,pr,betar,alphar,ar,gammar,gamma1,gamma2)

% forUBLE PRECISION, INTENT(IN)  :: ul,pl,al,gammal,rhor,ur,vr,pr,betar,alphar,ar,gammar,gamma1,gamma2
% forUBLE PRECISION, INTENT(OUT) :: rhosr,usr,vsr,psr,betasr,alphasr
% INTEGER, PARAMETER :: NRK4 = 500
% INTEGER :: m,k
% forUBLE PRECISION :: diffu,dua
% forUBLE PRECISION :: usrnew,usrold,asr,duaold
% forUBLE PRECISION :: rho,u,v,p,beta,alpha,dur
% forUBLE PRECISION :: rho_k1,rho_k2,rho_k3,rho_k4,p_k1,p_k2,p_k3,p_k4,alpha_k1,alpha_k2,alpha_k3,alpha_k4
% forUBLE PRECISION :: gamma,uguess
% forUBLE PRECISION, PARAMETER :: eps = 1D-15


%Runga-Kutta Solver
gamma  = (gammar+gammal)/2D0;
uguess = (ul/al*pl^((gamma-1D0)/(2D0*gamma))+ur/ar*pr^((gamma-1D0)/(2D0*gamma))+2D0/(gamma-1)*(pl^((gamma-1D0)/(2D0*gamma))-pl^((gamma-1D0)/(2D0*gamma))))/(pr^((gamma-1D0)/(2D0*gamma))/ar+pr^((gamma-1D0)/(2D0*gamma))/ar);


%-----------------------------------------------------------------------%
%                       Runga-Kutta fourth order                        %
%-----------------------------------------------------------------------%
diffu = 1D0;
m = 0;

while diffu>1D-10
m=m+1;

%-----------------------------------------------------------------------%
%                               RIGHT WAVE                              %
%-----------------------------------------------------------------------%
if m==1
usr = uguess;
elseif m==2
usr = 0.8D0*uguess;
end

rho   = rhor;
v     = vr;
p     = pr;
beta  = betar;
alpha = alphar;
dur   = (usr - ur)/NRK4;

for k=1:NRK4
    CALL RK4_k1_R_x(rho,p,alpha,gamma1,gamma2,dur,p_k1,rho_k1,alpha_k1)
    CALL RK4_k2_R_x(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dur,p_k2,rho_k2,alpha_k2)
    CALL RK4_k3_R_x(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dur,p_k3,rho_k3,alpha_k3)
    CALL RK4_k4_R_x(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dur,p_k4,rho_k4,alpha_k4)
    p = p + (p_k1 + 2D0*p_k2 + 2D0*p_k3 + p_k4)/6D0;
    rho = rho + (rho_k1 + 2D0*rho_k2 + 2D0*rho_k3 + rho_k4)/6D0;
    alpha = alpha + (alpha_k1 + 2D0*alpha_k2 + 2D0*alpha_k3 + alpha_k4)/6D0;
end

rhosr   = rho;
vsr     = v;
psr     = p;
betasr  = beta;
alphasr = alpha;
gamma   = 1D0/(alpha/(gamma1-1D0)+(1D0-alpha)/(gamma2-1D0))+1D0;
asr     = sqrt(gamma*psr/rhosr);

dua     = usr + asr;
diffu   = abs(dua);

if m==1
usrold = usr;
duaold = dua;
else
usrnew = usr - dua*(usr-usrold+eps)/(dua-duaold+eps);
usrold = usr;
usr    = usrnew;
duaold = dua;
end

end %while loop


end % osher_sonic_R_x



SUBROUTINE RK4_k1_L_x(rho,p,alpha,gamma1,gamma2,dul,p_k1,rho_k1,alpha_k1)
%------------------------------------------------%
%                                                %
% RK4 k1 step for left wave                      %
%                                                %
%------------------------------------------------%

IMPLICIT NONE

forUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,gamma1,gamma2,dul
forUBLE PRECISION, INTENT(OUT) :: p_k1,rho_k1,alpha_k1
forUBLE PRECISION :: gamma,a,b

rho_k1   = rho
p_k1     = p
alpha_k1 = alpha
gamma    = 1D0/(alpha_k1/(gamma1-1D0)+(1D0-alpha_k1)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k1/rho_k1)
b        = alpha_k1*(1D0-alpha_k1)*(gamma1-gamma2)/((1D0-alpha_k1)*gamma1+alpha_k1*gamma2)
p_k1     = -dul*rho_k1*a
rho_k1   = -dul*rho_k1/a
alpha_k1 = -dul*b/a

end SUBROUTINE RK4_k1_L_x


SUBROUTINE RK4_k1_R_x(rho,p,alpha,gamma1,gamma2,dur,p_k1,rho_k1,alpha_k1)
%------------------------------------------------%
%                                                %
% RK4 k1 step for right wave                     %
%                                                %
%------------------------------------------------%

IMPLICIT NONE

forUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,gamma1,gamma2,dur
forUBLE PRECISION, INTENT(OUT) :: p_k1,rho_k1,alpha_k1
forUBLE PRECISION :: gamma,a,b

p_k1     = p
alpha_k1 = alpha
rho_k1   = rho
gamma    = 1D0/(alpha_k1/(gamma1-1D0)+(1D0-alpha_k1)/(gamma2-1))+1D0
a        = sqrt(gamma*p_k1/rho_k1)
b        = alpha_k1*(1D0-alpha_k1)*(gamma1-gamma2)/((1D0-alpha_k1)*gamma1+alpha_k1*gamma2)
p_k1     = dur*rho_k1*a
rho_k1   = dur*rho_k1/a
alpha_k1 = dur*b/a

end SUBROUTINE RK4_k1_R_x


SUBROUTINE RK4_k2_L_x(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dul,p_k2,rho_k2,alpha_k2)
%------------------------------------------------%
%                                                %
% RK4 k2 step for left wave                      %
%                                                %
%------------------------------------------------%

IMPLICIT NONE

forUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dul
forUBLE PRECISION, INTENT(OUT) :: p_k2,rho_k2,alpha_k2
forUBLE PRECISION :: gamma,a,b

rho_k2   = rho+rho_k1/2D0
p_k2     = p+p_k1/2D0
alpha_k2 = alpha+alpha_k1/2D0
gamma    = 1D0/(alpha_k2/(gamma1-1D0)+(1D0-alpha_k2)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k2/rho_k2)
b        = alpha_k2*(1-alpha_k2)*(gamma1-gamma2)/((1D0-alpha_k2)*gamma1+alpha_k2*gamma2)
p_k2     = -dul*rho_k2*a
rho_k2   = -dul*rho_k2/a
alpha_k2 = -dul*b/a

end SUBROUTINE RK4_k2_L_x


SUBROUTINE RK4_k2_R_x(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dur,p_k2,rho_k2,alpha_k2)
%------------------------------------------------%
%                                                %
% RK4 k2 step for right wave                     %
%                                                %
%------------------------------------------------%

IMPLICIT NONE

forUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dur
forUBLE PRECISION, INTENT(OUT) :: p_k2,rho_k2,alpha_k2
forUBLE PRECISION :: gamma,a,b

p_k2     = p+p_k1/2D0
alpha_k2 = alpha+alpha_k1/2D0
rho_k2   = rho+rho_k1/2D0
gamma    = 1D0/(alpha_k2/(gamma1-1D0)+(1D0-alpha_k2)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k2/rho_k2)
b        = alpha_k2*(1D0-alpha_k2)*(gamma1-gamma2)/((1D0-alpha_k2)*gamma1+alpha_k2*gamma2)
p_k2     = dur*rho_k2*a
rho_k2   = dur*rho_k2/a
alpha_k2 = dur*b/a

end SUBROUTINE RK4_k2_R_x


SUBROUTINE RK4_k3_L_x(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dul,p_k3,rho_k3,alpha_k3)
%------------------------------------------------%
%                                                %
% RK4 k3 step for left wave                      %
%                                                %
%------------------------------------------------%

IMPLICIT NONE

forUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dul
forUBLE PRECISION, INTENT(OUT) :: p_k3,rho_k3,alpha_k3
forUBLE PRECISION :: gamma,a,b

p_k3     = p+p_k2/2D0
alpha_k3 = alpha+alpha_k2/2D0
rho_k3   = rho+rho_k2/2D0
gamma    = 1D0/(alpha_k3/(gamma1-1D0)+(1D0-alpha_k3)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k3/rho_k3)
b        = alpha_k3*(1D0-alpha_k3)*(gamma1-gamma2)/((1D0-alpha_k3)*gamma1+alpha_k3*gamma2)
p_k3     = -dul*rho_k3*a
rho_k3   = -dul*rho_k3/a
alpha_k3 = -dul*b/a

end SUBROUTINE RK4_k3_L_x


SUBROUTINE RK4_k3_R_x(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dur,p_k3,rho_k3,alpha_k3)
%------------------------------------------------%
%                                                %
% RK4 k3 step for right wave                     %
%                                                %
%------------------------------------------------%

IMPLICIT NONE

forUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dur
forUBLE PRECISION, INTENT(OUT) :: p_k3,rho_k3,alpha_k3
forUBLE PRECISION :: gamma,a,b

p_k3     = p+p_k2/2D0
alpha_k3 = alpha+alpha_k2/2D0
rho_k3   = rho+rho_k2/2D0
gamma    = 1D0/(alpha_k3/(gamma1-1D0)+(1D0-alpha_k3)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k3/rho_k3)
b        = alpha_k3*(1D0-alpha_k3)*(gamma1-gamma2)/((1D0-alpha_k3)*gamma1+alpha_k3*gamma2)
p_k3     = dur*rho_k3*a
rho_k3   = dur*rho_k3/a
alpha_k3 = dur*b/a

end SUBROUTINE RK4_k3_R_x




SUBROUTINE RK4_k4_L_x(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dul,p_k4,rho_k4,alpha_k4)
%------------------------------------------------%
%                                                %
% RK4 k4 step for left wave                      %
%                                                %
%------------------------------------------------%

IMPLICIT NONE

forUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dul
forUBLE PRECISION, INTENT(OUT) :: p_k4,rho_k4,alpha_k4
forUBLE PRECISION :: gamma,a,b

p_k4     = p+p_k3
alpha_k4 = alpha+alpha_k3
rho_k4   = rho+rho_k3
gamma    = 1D0/(alpha_k4/(gamma1-1D0)+(1D0-alpha_k4)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k4/rho_k4)
b        = alpha_k4*(1D0-alpha_k4)*(gamma1-gamma2)/((1D0-alpha_k4)*gamma1+alpha_k4*gamma2)
p_k4     = -dul*rho_k4*a
rho_k4   = -dul*rho_k4/a
alpha_k4 = -dul*b/a

end SUBROUTINE RK4_k4_L_x



SUBROUTINE RK4_k4_R_x(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dur,p_k4,rho_k4,alpha_k4)
%------------------------------------------------%
%                                                %
% RK4 k4 step for right wave                     %
%                                                %
%------------------------------------------------%

IMPLICIT NONE

forUBLE PRECISION, INTENT(IN)  :: rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dur
forUBLE PRECISION, INTENT(OUT) :: p_k4,rho_k4,alpha_k4
forUBLE PRECISION :: gamma,a,b

p_k4     = p+p_k3
alpha_k4 = alpha+alpha_k3
rho_k4   = rho+rho_k3
gamma    = 1D0/(alpha_k4/(gamma1-1D0)+(1D0-alpha_k4)/(gamma2-1D0))+1D0
a        = sqrt(gamma*p_k4/rho_k4)
b        = alpha_k4*(1D0-alpha_k4)*(gamma1-gamma2)/((1D0-alpha_k4)*gamma1+alpha_k4*gamma2)
p_k4     = dur*rho_k4*a
rho_k4   = dur*rho_k4/a
alpha_k4 = dur*b/a

end SUBROUTINE RK4_k4_R_x




SUBROUTINE sub_Sint12_x(rhol,ul,pl,betal,alphal,uf,gamma1,gamma2,sint12)
%--------------------------------------------------%
%                                                  %
%    Source in isentropic wave 1 to 2              %
%                                                  %
%--------------------------------------------------%

IMPLICIT NONE

forUBLE PRECISION, INTENT(IN)  :: rhol,ul,pl,betal,alphal,uf,gamma1,gamma2
forUBLE PRECISION, INTENT(OUT) :: sint12
forUBLE PRECISION :: dul,rho,u,p,beta,alpha,gamma,a,bS
forUBLE PRECISION :: rho_k1,rho_k2,rho_k3,rho_k4,p_k1,p_k2,p_k3,p_k4,alpha_k1,alpha_k2,alpha_k3,alpha_k4
INTEGER, PARAMETER :: NRK4 = 500
INTEGER :: k

sint12 = 0D0
rho   = rhol
u     = ul
p     = pl
beta  = betal
alpha = alphal
dul   = (uf - ul)/NRK4

for k=1,NRK4
    CALL RK4_k1_L_x(rho,p,alpha,gamma1,gamma2,dul,p_k1,rho_k1,alpha_k1)
    CALL RK4_k2_L_x(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dul,p_k2,rho_k2,alpha_k2)
    CALL RK4_k3_L_x(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dul,p_k3,rho_k3,alpha_k3)
    CALL RK4_k4_L_x(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dul,p_k4,rho_k4,alpha_k4)
    p = p + (p_k1 + 2D0*p_k2 + 2D0*p_k3 + p_k4)/6D0
    rho = rho + (rho_k1 + 2D0*rho_k2 + 2D0*rho_k3 + rho_k4)/6D0
    alpha = alpha + (alpha_k1 + 2D0*alpha_k2 + 2D0*alpha_k3 + alpha_k4)/6D0
    gamma  = 1D0/(alpha/(gamma1-1D0)+(1D0-alpha)/(gamma2-1D0))+1D0
    a = sqrt(gamma*p/rho)
    bS = alpha*(1D0-alpha)*(gamma1-gamma2)/((1D0-alpha)*gamma1+alpha*gamma2)
    sint12 = sint12 + (bS*p*(1D0-u/a)-(alpha-beta)*rho*a*u)*dul
end for

end SUBROUTINE sub_Sint12_x



SUBROUTINE sub_Sint1sl_x(rhol,ul,pl,betal,alphal,usl,gamma1,gamma2,sint1sl)
%--------------------------------------------------%
%                                                  %
%    Source in isentropic wave 1 to sonic line     %
%                                                  %
%--------------------------------------------------%

IMPLICIT NONE

forUBLE PRECISION, INTENT(IN)  :: rhol,ul,pl,betal,alphal,usl,gamma1,gamma2
forUBLE PRECISION, INTENT(OUT) :: sint1sl
forUBLE PRECISION :: dul,rho,u,p,beta,alpha,gamma,a,bS
forUBLE PRECISION :: rho_k1,rho_k2,rho_k3,rho_k4,p_k1,p_k2,p_k3,p_k4,alpha_k1,alpha_k2,alpha_k3,alpha_k4
INTEGER, PARAMETER :: NRK4 = 500
INTEGER :: k

sint1sl = 0D0
rho   = rhol
u     = ul
p     = pl
beta  = betal
alpha = alphal
dul   = (usl - ul)/NRK4

for k=1,NRK4
    CALL RK4_k1_L_x(rho,p,alpha,gamma1,gamma2,dul,p_k1,rho_k1,alpha_k1)
    CALL RK4_k2_L_x(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dul,p_k2,rho_k2,alpha_k2)
    CALL RK4_k3_L_x(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dul,p_k3,rho_k3,alpha_k3)
    CALL RK4_k4_L_x(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dul,p_k4,rho_k4,alpha_k4)
    p = p + (p_k1 + 2D0*p_k2 + 2D0*p_k3 + p_k4)/6D0
    rho = rho + (rho_k1 + 2D0*rho_k2 + 2D0*rho_k3 + rho_k4)/6D0
    alpha = alpha + (alpha_k1 + 2D0*alpha_k2 + 2D0*alpha_k3 + alpha_k4)/6D0
    gamma  = 1D0/(alpha/(gamma1-1D0)+(1D0-alpha)/(gamma2-1D0))+1D0
    a = sqrt(gamma*p/rho)
    bS = alpha*(1D0-alpha)*(gamma1-gamma2)/((1D0-alpha)*gamma1+alpha*gamma2)
    sint1sl = sint1sl + (bS*p*(1D0-u/a)-(alpha-beta)*rho*a*u)*dul
end for

end SUBROUTINE sub_Sint1sl_x


SUBROUTINE sub_Sintsl2_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,sintsl2)
%--------------------------------------------------%
%                                                  %
%    Source in isentropic wave sonic line to 2     %
%                                                  %
%--------------------------------------------------%

IMPLICIT NONE

forUBLE PRECISION, INTENT(IN)  :: rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2
forUBLE PRECISION, INTENT(OUT) :: sintsl2
forUBLE PRECISION :: dul,rho,u,p,beta,alpha,gamma,a,bS
forUBLE PRECISION :: rho_k1,rho_k2,rho_k3,rho_k4,p_k1,p_k2,p_k3,p_k4,alpha_k1,alpha_k2,alpha_k3,alpha_k4
INTEGER, PARAMETER :: NRK4 = 500
INTEGER :: k

sintsl2 = 0D0
rho   = rhosl
u     = usl
p     = psl
beta  = betasl
alpha = alphasl
dul   = (uf - usl)/NRK4

for k=1,NRK4
    CALL RK4_k1_L_x(rho,p,alpha,gamma1,gamma2,dul,p_k1,rho_k1,alpha_k1)
    CALL RK4_k2_L_x(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dul,p_k2,rho_k2,alpha_k2)
    CALL RK4_k3_L_x(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dul,p_k3,rho_k3,alpha_k3)
    CALL RK4_k4_L_x(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dul,p_k4,rho_k4,alpha_k4)
    p = p + (p_k1 + 2D0*p_k2 + 2D0*p_k3 + p_k4)/6D0
    rho = rho + (rho_k1 + 2D0*rho_k2 + 2D0*rho_k3 + rho_k4)/6D0
    alpha = alpha + (alpha_k1 + 2D0*alpha_k2 + 2D0*alpha_k3 + alpha_k4)/6D0
    gamma  = 1D0/(alpha/(gamma1-1D0)+(1D0-alpha)/(gamma2-1D0))+1D0
    a = sqrt(gamma*p/rho)
    bS = alpha*(1D0-alpha)*(gamma1-gamma2)/((1D0-alpha)*gamma1+alpha*gamma2)
    sintsl2 = sintsl2 + (bS*p*(1D0-u/a)-(alpha-beta)*rho*a*u)*dul
end for

end SUBROUTINE sub_Sintsl2_x



SUBROUTINE sub_SintCD_x(uf,pf,alpha3,alpha2,sintCD)
%--------------------------------------------------%
%                                                  %
%    Source over contact discontinuity             %
%                                                  %
%--------------------------------------------------%

IMPLICIT NONE

forUBLE PRECISION, INTENT(IN)  :: uf,pf,alpha3,alpha2
forUBLE PRECISION, INTENT(OUT) :: sintCD

sintCD = uf*pf*(alpha3-alpha2)

end SUBROUTINE sub_sintCD_x



SUBROUTINE sub_Sint34_x(rho3,uf,pf,beta3,alpha3,ur,gamma1,gamma2,sint34)
%--------------------------------------------------%
%                                                  %
%    Source in isentropic wave 3 to 4              %
%                                                  %
%--------------------------------------------------%

IMPLICIT NONE

forUBLE PRECISION, INTENT(IN)  :: rho3,uf,pf,beta3,alpha3,ur,gamma1,gamma2
forUBLE PRECISION, INTENT(OUT) :: sint34
forUBLE PRECISION :: dur,rho,u,p,beta,alpha,gamma,a,bS
forUBLE PRECISION :: rho_k1,rho_k2,rho_k3,rho_k4,p_k1,p_k2,p_k3,p_k4,alpha_k1,alpha_k2,alpha_k3,alpha_k4
INTEGER, PARAMETER :: NRK4 = 500
INTEGER :: k

sint34 = 0D0
rho   = rho3
u     = uf
p     = pf
beta  = beta3
alpha = alpha3
dur   = (ur - uf)/NRK4

for k=1,NRK4
    CALL RK4_k1_R_x(rho,p,alpha,gamma1,gamma2,dur,p_k1,rho_k1,alpha_k1)
    CALL RK4_k2_R_x(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dur,p_k2,rho_k2,alpha_k2)
    CALL RK4_k3_R_x(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dur,p_k3,rho_k3,alpha_k3)
    CALL RK4_k4_R_x(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dur,p_k4,rho_k4,alpha_k4)
    p = p + (p_k1 + 2D0*p_k2 + 2D0*p_k3 + p_k4)/6D0
    rho = rho + (rho_k1 + 2D0*rho_k2 + 2D0*rho_k3 + rho_k4)/6D0
    alpha = alpha + (alpha_k1 + 2D0*alpha_k2 + 2D0*alpha_k3 + alpha_k4)/6D0
    gamma  = 1D0/(alpha/(gamma1-1D0)+(1D0-alpha)/(gamma2-1D0))+1D0
    a = sqrt(gamma*p/rho)
    bS = alpha*(1D0-alpha)*(gamma1-gamma2)/((1D0-alpha)*gamma1+alpha*gamma2)
    sint34 = sint34 + (bS*p*(1D0+u/a)+(alpha-beta)*rho*a*u)*dur
end for

end SUBROUTINE sub_Sint34_x




SUBROUTINE sub_Sint3sr_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,sint3sr)
%--------------------------------------------------%
%                                                  %
%    Source in isentropic wave 3 to sonic line     %
%                                                  %
%--------------------------------------------------%

IMPLICIT NONE

forUBLE PRECISION, INTENT(IN)  :: rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2
forUBLE PRECISION, INTENT(OUT) :: sint3sr
forUBLE PRECISION :: dur,rho,u,p,beta,alpha,gamma,a,bS
forUBLE PRECISION :: rho_k1,rho_k2,rho_k3,rho_k4,p_k1,p_k2,p_k3,p_k4,alpha_k1,alpha_k2,alpha_k3,alpha_k4
INTEGER, PARAMETER :: NRK4 = 500
INTEGER :: k

sint3sr = 0D0
rho   = rho3
u     = uf
p     = pf
beta  = beta3
alpha = alpha3
dur   = (usr - uf)/NRK4

for k=1,NRK4
    CALL RK4_k1_R_x(rho,p,alpha,gamma1,gamma2,dur,p_k1,rho_k1,alpha_k1)
    CALL RK4_k2_R_x(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dur,p_k2,rho_k2,alpha_k2)
    CALL RK4_k3_R_x(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dur,p_k3,rho_k3,alpha_k3)
    CALL RK4_k4_R_x(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dur,p_k4,rho_k4,alpha_k4)
    p = p + (p_k1 + 2D0*p_k2 + 2D0*p_k3 + p_k4)/6D0
    rho = rho + (rho_k1 + 2D0*rho_k2 + 2D0*rho_k3 + rho_k4)/6D0
    alpha = alpha + (alpha_k1 + 2D0*alpha_k2 + 2D0*alpha_k3 + alpha_k4)/6D0
    gamma  = 1D0/(alpha/(gamma1-1D0)+(1D0-alpha)/(gamma2-1D0))+1D0
    a = sqrt(gamma*p/rho)
    bS = alpha*(1D0-alpha)*(gamma1-gamma2)/((1D0-alpha)*gamma1+alpha*gamma2)
    sint3sr = sint3sr + (bS*p*(1D0+u/a)+(alpha-beta)*rho*a*u)*dur
end for

end SUBROUTINE sub_Sint3sr_x




SUBROUTINE sub_Sintsr4_x(rhosr,usr,psr,betasr,alphasr,ur,gamma1,gamma2,sintsr4)
%--------------------------------------------------%
%                                                  %
%    Source in isentropic wave sonic line to 4     %
%                                                  %
%--------------------------------------------------%

IMPLICIT NONE

forUBLE PRECISION, INTENT(IN)  :: rhosr,usr,psr,betasr,alphasr,ur,gamma1,gamma2
forUBLE PRECISION, INTENT(OUT) :: sintsr4
forUBLE PRECISION :: dur,rho,u,p,beta,alpha,gamma,a,bS
forUBLE PRECISION :: rho_k1,rho_k2,rho_k3,rho_k4,p_k1,p_k2,p_k3,p_k4,alpha_k1,alpha_k2,alpha_k3,alpha_k4
INTEGER, PARAMETER :: NRK4 = 500
INTEGER :: k

sintsr4 = 0D0
rho   = rhosr
u     = usr
p     = psr
alpha = alphasr
beta  = betasr
dur   = (ur - usr)/NRK4

for k=1,NRK4
    CALL RK4_k1_R_x(rho,p,alpha,gamma1,gamma2,dur,p_k1,rho_k1,alpha_k1)
    CALL RK4_k2_R_x(rho,p,alpha,rho_k1,p_k1,alpha_k1,gamma1,gamma2,dur,p_k2,rho_k2,alpha_k2)
    CALL RK4_k3_R_x(rho,p,alpha,rho_k2,p_k2,alpha_k2,gamma1,gamma2,dur,p_k3,rho_k3,alpha_k3)
    CALL RK4_k4_R_x(rho,p,alpha,rho_k3,p_k3,alpha_k3,gamma1,gamma2,dur,p_k4,rho_k4,alpha_k4)
    p = p + (p_k1 + 2D0*p_k2 + 2D0*p_k3 + p_k4)/6D0
    rho = rho + (rho_k1 + 2D0*rho_k2 + 2D0*rho_k3 + rho_k4)/6D0
    alpha = alpha + (alpha_k1 + 2D0*alpha_k2 + 2D0*alpha_k3 + alpha_k4)/6D0
    gamma  = 1D0/(alpha/(gamma1-1D0)+(1D0-alpha)/(gamma2-1D0))+1D0
    a = sqrt(gamma*p/rho)
    bS = alpha*(1D0-alpha)*(gamma1-gamma2)/((1D0-alpha)*gamma1+alpha*gamma2)
    sintsr4 = sintsr4 + (bS*p*(1D0+u/a)+(alpha-beta)*rho*a*u)*dur
end for

end SUBROUTINE sub_Sintsr4_x


end MODULE mod_flux_x
