function [F,sL,sR]=linosherflux_x(NX,NY,FBC,Wfx,gamma1,gamma2)

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

F = FBC;

sL = zeros(NY,NX+1);
sR = zeros(NY,NX+1);

for i=1:NY
for j=2:NX

maxdiffw = 0e0;
for k=3:6
diffw = abs(Wfx(i,2*(j-1),k)-Wfx(i,2*j-1,k));
maxdiffw = max(maxdiffw,diffw);
end

if maxdiffw<1e-6

F(i,j,1) = Wfx(i,2*(j-1),3)*Wfx(i,2*(j-1),4);
F(i,j,2) = Wfx(i,2*(j-1),3)*Wfx(i,2*(j-1),4)^(2e0)+Wfx(i,2*(j-1),6);
F(i,j,3) = Wfx(i,2*(j-1),3)*Wfx(i,2*(j-1),4)*Wfx(i,2*(j-1),5);
F(i,j,4) = (Wfx(i,2*(j-1),8)/(gamma1-1e0)+(1e0-Wfx(i,2*(j-1),8))/(gamma2-1e0)+1e0)*Wfx(i,2*(j-1),4)*Wfx(i,2*(j-1),6)+0.5e0*Wfx(i,2*(j-1),3)*Wfx(i,2*(j-1),4)*(Wfx(i,2*(j-1),4)^(2e0)+Wfx(i,2*(j-1),5)^(2e0));
F(i,j,5) = Wfx(i,2*(j-1),7)*Wfx(i,2*(j-1),3)*Wfx(i,2*(j-1),4);
F(i,j,6) = gamma1/(gamma1-1e0)*Wfx(i,2*(j-1),4)*Wfx(i,2*(j-1),6)*Wfx(i,2*(j-1),8)+1e0/2e0*Wfx(i,2*(j-1),3)*Wfx(i,2*(j-1),4)*(Wfx(i,2*(j-1),4)^(2e0)+Wfx(i,2*(j-1),5)^(2e0))*Wfx(i,2*(j-1),7);

sL(i,j)  = 0e0;
sR(i,j)  = 0e0;

elseif (abs(Wfx(i,2*(j-1),6)-Wfx(i,2*j-1,6))<1e-6 && abs(Wfx(i,2*(j-1),4))<1e-6 && abs(Wfx(i,2*j-1,4))<1e-6)

F(i,j,1) = 0e0;
F(i,j,2) = Wfx(i,2*(j-1),6);
F(i,j,3) = 0e0;
F(i,j,4) = 0e0;
F(i,j,5) = 0e0;
F(i,j,6) = 0e0;

sL(i,j)  = 0e0;
sR(i,j)  = 0e0;

else

rho1   = Wfx(i,2*(j-1),3);
u1     = Wfx(i,2*(j-1),4);
v1     = Wfx(i,2*(j-1),5);
p1     = Wfx(i,2*(j-1),6);
beta1  = Wfx(i,2*(j-1),7);
alpha1 = Wfx(i,2*(j-1),8);
a1     = sqrt(1.0e0/(alpha1/gamma1+(1.0e0-alpha1)/gamma2)*p1/rho1);

rho4   = Wfx(i,2*j-1,3);
u4     = Wfx(i,2*j-1,4);
v4     = Wfx(i,2*j-1,5);
p4     = Wfx(i,2*j-1,6);
beta4  = Wfx(i,2*j-1,7);
alpha4 = Wfx(i,2*j-1,8);
a4     = sqrt(1.0e0/(alpha4/gamma1+(1.0e0-alpha4)/gamma2)*p4/rho4);

% Final state values
[rho2,rho3,uf,vf2,vf3,pf,beta2,beta3,alpha2,alpha3]=linoshersolver_x(rho1,u1,v1,p1,beta1,alpha1,a1,rho4,u4,v4,p4,beta4,alpha4,a4,gamma1,gamma2);

% speed of sound
a2 = sqrt(1e0/(alpha2/gamma1+(1e0-alpha2)/gamma2)*pf/rho2);
a3 = sqrt(1e0/(alpha3/gamma1+(1e0-alpha3)/gamma2)*pf/rho3);

% Sonic point values
[rhosl,usl,vsl,psl,betasl,alphasl]=linosher_sonic_L_x(rho1,u1,v1,p1,beta1,alpha1,a1,gamma1,gamma2);
[rhosr,usr,vsr,psr,betasr,alphasr]=linosher_sonic_R_x(rho4,u4,v4,p4,beta4,alpha4,a4,gamma1,gamma2);

Fl(1) = rho1*u1;
Fl(2) = rho1*u1^(2e0)+p1;
Fl(3) = rho1*u1*v1;
Fl(4) = (alpha1/(gamma1-1e0)+(1e0-alpha1)/(gamma2-1e0)+1e0)*u1*p1+1e0/2e0*rho1*u1*(u1^(2e0)+v1^(2e0));
Fl(5) = beta1*rho1*u1;
Fl(6) = gamma1/(gamma1-1e0)*u1*p1*alpha1+1e0/2e0*rho1*u1*(u1^(2e0)+v1^(2e0))*beta1;

F2(1) = rho2*uf;
F2(2) = rho2*uf^(2e0)+pf;
F2(3) = rho2*uf*vf2;
F2(4) = (alpha2/(gamma1-1e0)+(1e0-alpha2)/(gamma2-1e0)+1e0)*uf*pf+1e0/2e0*rho2*uf*(uf^(2e0)+vf2^(2e0));
F2(5) = beta2*rho2*uf;
F2(6) = gamma1/(gamma1-1e0)*uf*pf*alpha2+1e0/2e0*rho2*uf*(uf^(2e0)+vf2^(2e0))*beta2;

F3(1) = rho3*uf;
F3(2) = rho3*uf^(2e0)+pf;
F3(3) = rho3*uf*vf3;
F3(4) = (alpha3/(gamma1-1e0)+(1e0-alpha3)/(gamma2-1e0)+1e0)*uf*pf+1e0/2e0*rho3*uf*(uf^(2e0)+vf3^(2e0));
F3(5) = beta3*rho3*uf;
F3(6) = gamma1/(gamma1-1e0)*uf*pf*alpha3+1e0/2e0*rho3*uf*(uf^(2e0)+vf3^(2e0))*beta3;

Fr(1) = rho4*u4;
Fr(2) = rho4*u4^(2e0)+p4;
Fr(3) = rho4*u4*v4;
Fr(4) = (alpha4/(gamma1-1e0)+(1e0-alpha4)/(gamma2-1e0)+1e0)*u4*p4+1e0/2e0*rho4*u4*(u4^(2e0)+v4^(2e0));
Fr(5) = beta4*rho4*u4;
Fr(6) = gamma1/(gamma1-1e0)*u4*p4*alpha4+1e0/2e0*rho4*u4*(u4^(2e0)+v4^(2e0))*beta4;

Fsl(1) = rhosl*usl;
Fsl(2) = rhosl*usl^(2e0)+psl;
Fsl(3) = rhosl*usl*vsl;
Fsl(4) = (alphasl/(gamma1-1e0)+(1e0-alphasl)/(gamma2-1e0)+1e0)*usl*psl+1e0/2e0*rhosl*usl*(usl^(2e0)+vsl^(2e0));
Fsl(5) = betasl*rhosl*usl;
Fsl(6) = gamma1/(gamma1-1e0)*usl*psl*alphasl+1e0/2e0*rhosl*usl*(usl^(2e0)+vsl^(2e0))*betasl;

Fsr(1) = rhosr*usr;
Fsr(2) = rhosr*usr^(2e0)+psr;
Fsr(3) = rhosr*usr*vsr;
Fsr(4) = (alphasr/(gamma1-1e0)+(1e0-alphasr)/(gamma2-1e0)+1e0)*usr*psr+1e0/2e0*rhosr*usr*(usr^(2e0)+vsr^(2e0));
Fsr(5) = betasr*rhosr*usr;
Fsr(6) = gamma1/(gamma1-1e0)*usr*psr*alphasr+1e0/2e0*rhosr*usr*(usr^(2e0)+vsr^(2e0))*betasr;

lambda0 = uf;
lambda1 = u1-a1;
lambda2 = uf-a2;
lambda3 = uf+a3;
lambda4 = u4+a4;
% Lambda  = max(abs(lambda0),abs(lambda1),abs(lambda2),abs(lambda3),abs(lambda4),abs(Lambda)))

if (lambda0>=0e0 && lambda2>=0e0)
    if (lambda1>=0e0 && lambda4>=0e0)
        F(i,j,:) = Fl;
        [sint12] = sub_Sint_lin_x(rho1,u1,p1,beta1,alpha1,uf,gamma1,gamma2,-1);
        [sintCD] = sub_SintCD_x(uf,pf,alpha3,alpha2);
        [sint34] = sub_Sint_lin_x(rho3,uf,pf,beta3,alpha3,u4,gamma1,gamma2,+1);
        sL(i,j) = 0e0;
        sR(i,j) = sint12+sintCD+sint34;
    elseif (lambda1>=0e0 && lambda4<=0e0)
        F(i,j,:) = Fl+Fr-Fsr;
        [sintsr4] = sub_Sint_lin_x(rhosr,usr,psr,betasr,alphasr,u4,gamma1,gamma2,+1);
        [sint12]  = sub_Sint_lin_x(rho1,u1,p1,beta1,alpha1,uf,gamma1,gamma2,-1);
        [sintCD]  = sub_SintCD_x(uf,pf,alpha3,alpha2);
        [sint3sr] = sub_Sint_lin_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,+1);
        sL(i,j) = sintsr4;
        sR(i,j) = sint12+sintCD+sint3sr;
    elseif (lambda1<=0e0 && lambda4>=0e0)
        F(i,j,:) = Fsl;
        [sint1sl] = sub_Sint_lin_x(rho1,u1,p1,beta1,alpha1,usl,gamma1,gamma2,-1);
        [sintsl2] = sub_Sint_lin_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,-1);
        [sintCD]  = sub_SintCD_x(uf,pf,alpha3,alpha2);
        [sint34]  = sub_Sint_lin_x(rho3,uf,pf,beta3,alpha3,u4,gamma1,gamma2,+1);
        sL(i,j) = sint1sl;
        sR(i,j) = sintsl2+sintCD+sint34;
    elseif (lambda1<=0e0 && lambda4<=0e0)
        F(i,j,:) = Fsl-Fsr+Fr;
        [sint1sl] = sub_Sint_lin_x(rho1,u1,p1,beta1,alpha1,usl,gamma1,gamma2,-1);
        [sintsr4] = sub_Sint_lin_x(rhosr,usr,psr,betasr,alphasr,u4,gamma1,gamma2,+1);
        [sintsl2] = sub_Sint_lin_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,-1);
        [sintCD]  = sub_SintCD_x(uf,pf,alpha3,alpha2);
        [sint3sr] = sub_Sint_lin_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,+1);
        sL(i,j) = sint1sl+sintsr4;
        sR(i,j) = sintsl2+sintCD+sint3sr;
    end

elseif  (lambda0>=0e0 && lambda2<=0e0)
    if (lambda1>=0e0 && lambda4>=0e0)
        F(i,j,:) = Fl-Fsl+F2;
        [sintsl2] = sub_Sint_lin_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,-1);
        [sint1sl] = sub_Sint_lin_x(rho1,u1,p1,beta1,alpha1,usl,gamma1,gamma2,-1);
        [sintCD]  = sub_SintCD_x(uf,pf,alpha3,alpha2);
        [sint34]  = sub_Sint_lin_x(rho3,uf,pf,beta3,alpha3,u4,gamma1,gamma2,+1);
        sL(i,j) = sintsl2;
        sR(i,j) = sint1sl+sintCD+sint34;
    elseif (lambda1>=0e0 && lambda4<=0e0)
        F(i,j,:) = Fl-Fsl+F2-Fsr+Fr;
        [sintsl2] = sub_Sint_lin_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,-1);
        [sintsr4] = sub_Sint_lin_x(rhosr,usr,psr,betasr,alphasr,u4,gamma1,gamma2,+1);
        [sint1sl] = sub_Sint_lin_x(rho1,u1,p1,beta1,alpha1,usl,gamma1,gamma2,-1);
        [sintCD]  = sub_SintCD_x(uf,pf,alpha3,alpha2);
        [sint3sr] = sub_Sint_lin_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,+1);
        sL(i,j) = sintsl2+sintsr4;
        sR(i,j) = sint1sl+sintCD+sint3sr;
    elseif (lambda1<=0e0 && lambda4>=0e0)
        F(i,j,:) = F2;
        [sint12] = sub_Sint_lin_x(rho1,u1,p1,beta1,alpha1,uf,gamma1,gamma2,-1);
        [sintCD] = sub_SintCD_x(uf,pf,alpha3,alpha2);
        [sint34] = sub_Sint_lin_x(rho3,uf,pf,beta3,alpha3,u4,gamma1,gamma2,+1);
        sL(i,j) = sint12;
        sR(i,j) = sintCD+sint34;
    elseif (lambda1<=0e0 && lambda4<=0e0)
        F(i,j,:) = Fr+F2-Fsr;
        [sint12]  = sub_Sint_lin_x(rho1,u1,p1,beta1,alpha1,uf,gamma1,gamma2,-1);
        [sintsr4] = sub_Sint_lin_x(rhosr,usr,psr,betasr,alphasr,u4,gamma1,gamma2,+1);
        [sintCD]  = sub_SintCD_x(uf,pf,alpha3,alpha2);
        [sint3sr] = sub_Sint_lin_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,+1);
        sL(i,j) = sint12+sintsr4;
        sR(i,j) = sintCD+sint3sr;
    end

elseif  (lambda0<=0e0 && lambda3>=0e0)
    if (lambda1>=0e0 && lambda4>=0e0)
        F(i,j,:) = Fl-Fsl+F3;
        [sintsl2] = sub_Sint_lin_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,-1);
        [sintCD]  = sub_SintCD_x(uf,pf,alpha3,alpha2);
        [sint1sl] = sub_Sint_lin_x(rho1,u1,p1,beta1,alpha1,usl,gamma1,gamma2,-1);
        [sint34]  = sub_Sint_lin_x(rho3,uf,pf,beta3,alpha3,u4,gamma1,gamma2,+1);
        sL(i,j) = sintsl2+sintCD;
        sR(i,j) = sint1sl+sint34;
    elseif (lambda1>=0e0 && lambda4<=0e0)
        F(i,j,:) = Fl-Fsl+F3-Fsr+Fr;
        [sintsl2] = sub_Sint_lin_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,-1);
        [sintCD]  = sub_SintCD_x(uf,pf,alpha3,alpha2);
        [sintsr4] = sub_Sint_lin_x(rhosr,usr,psr,betasr,alphasr,u4,gamma1,gamma2,+1);
        [sint1sl] = sub_Sint_lin_x(rho1,u1,p1,beta1,alpha1,usl,gamma1,gamma2,-1);
        [sint3sr] = sub_Sint_lin_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,+1);
        sL(i,j) = sintsl2+sintCD+sintsr4;
        sR(i,j) = sint1sl+sint3sr;
    elseif (lambda1<=0e0 && lambda4>=0e0)
        F(i,j,:) = F3;
        [sint12] = sub_Sint_lin_x(rho1,u1,p1,beta1,alpha1,uf,gamma1,gamma2,-1);
        [sintCD] = sub_SintCD_x(uf,pf,alpha3,alpha2);
        [sint34] = sub_Sint_lin_x(rho3,uf,pf,beta3,alpha3,u4,gamma1,gamma2,+1);
        sL(i,j) = sint12+sintCD;
        sR(i,j) = sint34;
    elseif (lambda1<=0e0 && lambda4<=0e0)
        F(i,j,:) = Fr+F3-Fsr;
        [sint12]  = sub_Sint_lin_x(rho1,u1,p1,beta1,alpha1,uf,gamma1,gamma2,-1);
        [sintCD]  = sub_SintCD_x(uf,pf,alpha3,alpha2);
        [sintsr4] = sub_Sint_lin_x(rhosr,usr,psr,betasr,alphasr,u4,gamma1,gamma2,+1);
        [sint3sr] = sub_Sint_lin_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,+1);
        sL(i,j) = sint12+sintCD+sintsr4;
        sR(i,j) = sint3sr;
    end

elseif  (lambda0<=0e0 && lambda3<=0e0)
    if (lambda1>=0e0 && lambda4>=0e0)
        F(i,j,:) = Fl-Fsl+Fsr;
        [sintsl2] = sub_Sint_lin_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,-1);
        [sintCD]  = sub_SintCD_x(uf,pf,alpha3,alpha2);
        [sint3sr] = sub_Sint_lin_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,+1);
        [sint1sl] = sub_Sint_lin_x(rho1,u1,p1,beta1,alpha1,usl,gamma1,gamma2,-1);
        [sintsr4] = sub_Sint_lin_x(rhosr,usr,psr,betasr,alphasr,u4,gamma1,gamma2,+1);
        sL(i,j) = sintsl2+sintCD+sint3sr;
        sR(i,j) = sint1sl+sintsr4;
    elseif (lambda1>=0e0 && lambda4<=0e0)
        F(i,j,:) = Fl-Fsl+Fr;
        [sintsl2] = sub_Sint_lin_x(rhosl,usl,psl,betasl,alphasl,uf,gamma1,gamma2,-1);
        [sintCD]  = sub_SintCD_x(uf,pf,alpha3,alpha2);
        [sint34]  = sub_Sint_lin_x(rho3,uf,pf,beta3,alpha3,u4,gamma1,gamma2,+1);
        [sint1sl] = sub_Sint_lin_x(rho1,u1,p1,beta1,alpha1,usl,gamma1,gamma2,-1);
        sL(i,j) = sintsl2+sintCD+sint34;
        sR(i,j) = sint1sl;
    elseif (lambda1<=0e0 && lambda4>=0e0)
        F(i,j,:) = Fsr;
        [sint12]  = sub_Sint_lin_x(rho1,u1,p1,beta1,alpha1,uf,gamma1,gamma2,-1);
        [sintCD]  = sub_SintCD_x(uf,pf,alpha3,alpha2);
        [sint3sr] = sub_Sint_lin_x(rho3,uf,pf,beta3,alpha3,usr,gamma1,gamma2,+1);
        [sintsr4] = sub_Sint_lin_x(rhosr,usr,psr,betasr,alphasr,u4,gamma1,gamma2,+1);
        sL(i,j) = sint12+sintCD+sint3sr;
        sR(i,j) = sintsr4;
    elseif (lambda1<=0e0 && lambda4<=0e0)
        F(i,j,:) = Fr;
        [sint12] = sub_Sint_lin_x(rho1,u1,p1,beta1,alpha1,uf,gamma1,gamma2,-1);
        [sintCD] = sub_SintCD_x(uf,pf,alpha3,alpha2);
        [sint34] = sub_Sint_lin_x(rho3,uf,pf,beta3,alpha3,u4,gamma1,gamma2,+1);
        sL(i,j) = sint12+sintCD+sint34;
        sR(i,j) = 0e0;
    end

end %flux choice

end %if L==R

end %for j loop
end %for i loop


end % linosherflux_x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rho2,rho3,uf,vf2,vf3,pf,beta2,beta3,alpha2,alpha3]=linoshersolver_x(rho1,u1,v1,p1,beta1,alpha1,a1,rho4,u4,v4,p4,beta4,alpha4,a4,gamma1,gamma2)

%Final state values
        uf    = (rho1*a1*u1+rho4*a4*u4+(p1-p4))/(rho1*a1+rho4*a4);
        pf    = (rho4*a4*p1+rho1*a1*p4+rho1*a1*rho4*a4*(u1-u4))/(rho1*a1+rho4*a4);
        rho2   = rho1+(pf-p1)/a1^2;
        rho3   = rho4+(pf-p4)/a4^2;
        vf2    = v1;
        vf3    = v4;
        beta2  = beta1;
        beta3  = beta4;
        alpha2 = alpha1+(pf-p1)*(1/gamma2-1/gamma1)*alpha1*(1-alpha1)/p1;
        alpha3 = alpha4+(pf-p4)*(1/gamma2-1/gamma1)*alpha4*(1-alpha4)/p4;
        
end % linoshersolver_x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rhosl,usl,vsl,psl,betasl,alphasl]=linosher_sonic_L_x(rho1,u1,v1,p1,beta1,alpha1,a1,gamma1,gamma2)

%------------------------------------%
%                                    %
% Solver for left sonic point        %
%                                    %
%------------------------------------%

uguess = a1;

%-----------------------------------------------------------------------%
%                         runga-kutta fou4th-order                      %
%-----------------------------------------------------------------------%
diffu = 1e0;
m = 0;

while diffu>1e-10
m=m+1;

%-----------------------------------------------------------------------%
%                               LEFT WAVE                               %
%-----------------------------------------------------------------------%
if m==1
usl = uguess;
elseif m==2
usl = 0.8e0*uguess;
end

psl     = p1-rho1*a1*(usl-u1);
rhosl   = rho1-rho1/a1*(usl-u1);
b1      = alpha1*(1-alpha1)*(gamma1-gamma2)/((1-alpha1)*gamma1+alpha1*gamma2);
alphasl = alpha1-b1/a1*(usl-u1);
asl     = sqrt(1/(alphasl/gamma1+(1-alphasl)/gamma2)*psl/rhosl);

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

psl     = p1-rho1*a1*(usl-u1);
rhosl   = rho1-rho1/a1*(usl-u1);
vsl     = v1;
betasl  = beta1;
b1      = alpha1*(1-alpha1)*(gamma1-gamma2)/((1-alpha1)*gamma1+alpha1*gamma2);
alphasl = alpha1-b1/a1*(usl-u1);

end % linosher_sonic_L_x



function [rhosr,usr,vsr,psr,betasr,alphasr]=linosher_sonic_R_x(rho4,u4,v4,p4,beta4,alpha4,a4,gamma1,gamma2)

%Runga-Kutta Solver
uguess = -a4;


%-----------------------------------------------------------------------%
%                       Runga-Kutta fourth order                        %
%-----------------------------------------------------------------------%
diffu = 1e0;
m = 0;

while diffu>1e-10
m=m+1;

%-----------------------------------------------------------------------%
%                               RIGHT WAVE                              %
%-----------------------------------------------------------------------%
if m==1
usr = uguess;
elseif m==2
usr = 0.8e0*uguess;
end

psr     = p4+rho4*a4*(usr-u4);
rhosr   = rho4+rho4/a4*(usr-u4);
b4      = alpha4*(1-alpha4)*(gamma1-gamma2)/((1-alpha4)*gamma1+alpha4*gamma2);
alphasr = alpha4+b4/a4*(usr-u4);
asr     = sqrt(1/(alphasr/gamma1+(1-alphasr)/gamma2)*psr/rhosr);

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

psr     = p4+rho4*a4*(usr-u4);
rhosr   = rho4+rho4/a4*(usr-u4);
vsr     = v4;
betasr  = beta4;
b4      = alpha4*(1-alpha4)*(gamma1-gamma2)/((1-alpha4)*gamma1+alpha4*gamma2);
alphasr = alpha4+b4/a4*(usr-u4);

end % linosher_sonic_R_x



function [sint]=sub_Sint_lin_x(rho,u,p,beta,alpha,u_end,gamma1,gamma2,lr)
%--------------------------------------------------%
%                                                  %
%    Source in isentropic wave                     %
%                                                  %
%--------------------------------------------------%

% lr is a sign function for left running wave (-) or right running wave (+)

du   = (u_end - u);
a    = sqrt(1e0/(alpha/gamma1+(1e0-alpha)/gamma2)*p/rho);
bS   = alpha*(1e0-alpha)*(gamma1-gamma2)/((1e0-alpha)*gamma1+alpha*gamma2);
sint = (bS*p*(1e0+lr*u/a)+lr*(alpha-beta)*rho*a*u)*du;

end % sub_Sint_lin_x


function [sintCD]=sub_SintCD_x(uf,pf,alpha3,alpha2)
%--------------------------------------------------%
%                                                  %
%    Sou4ce over contact discontinuity             %
%                                                  %
%--------------------------------------------------%

sintCD = uf*pf*(alpha3-alpha2);

end % sub_sintCD_x