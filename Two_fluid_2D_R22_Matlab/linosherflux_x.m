function [F,sL,sR]=linosherflux_x(F,Wfx)

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

global NX NY
global gamma1 gamma2

sL = zeros(NY,NX+1);
sR = zeros(NY,NX+1);

for i=1:NY
for j=2:NX

rho1   = Wfx(i,2*(j-1),3);
u1     = Wfx(i,2*(j-1),4);
v1     = Wfx(i,2*(j-1),5);
p1     = Wfx(i,2*(j-1),6);
beta1  = Wfx(i,2*(j-1),7);
alpha1 = Wfx(i,2*(j-1),8);
a1     = sqrt(1/(alpha1/gamma1+(1-alpha1)/gamma2)*p1/rho1);

rho4   = Wfx(i,2*j-1,3);
u4     = Wfx(i,2*j-1,4);
v4     = Wfx(i,2*j-1,5);
p4     = Wfx(i,2*j-1,6);
beta4  = Wfx(i,2*j-1,7);
alpha4 = Wfx(i,2*j-1,8);
a4     = sqrt(1/(alpha4/gamma1+(1-alpha4)/gamma2)*p4/rho4);


if abs(rho1-rho4)+abs(u1-u4)+abs(v1-v4)+abs(p1-p4)<1e-8

F(i,j,1) = rho1*u1;
F(i,j,2) = rho1*u1^2+p1;
F(i,j,3) = rho1*u1*v1;
F(i,j,4) = (alpha1/(gamma1-1)+(1-alpha1)/(gamma2-1)+1)*u1*p1+1/2*rho1*u1*(u1^2+v1^2);
F(i,j,5) = beta1*rho1*u1;
F(i,j,6) = gamma1/(gamma1-1)*u1*p1*alpha1+1/2*rho1*u1*(u1^2+v1^2)*beta1;

elseif ( abs(p1-p4) + abs(u1)<1e-6 + abs(u4) )<1e-8

F(i,j,2) = p1;

else

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

%speed of sound
a2 = sqrt(1/(alpha2/gamma1+(1-alpha2)/gamma2)*pf/rho2);
a3 = sqrt(1/(alpha3/gamma1+(1-alpha3)/gamma2)*pf/rho3);



% Sonic point values
[rhosl,usl,vsl,psl,betasl,alphasl]=linosher_sonic_L_x(rho1,u1,v1,p1,beta1,alpha1,a1);
[rhosr,usr,vsr,psr,betasr,alphasr]=linosher_sonic_R_x(rho4,u4,v4,p4,beta4,alpha4,a4);

Fl  = [rho1*u1  ;rho1*u1^2+p1   ;rho1*u1*v1   ;(alpha1/(gamma1-1)+(1-alpha1)/(gamma2-1)+1)*u1*p1+1/2*rho1*u1*(u1^2+v1^2)        ;beta1*rho1*u1   ;gamma1/(gamma1-1)*u1*p1*alpha1+1/2*rho1*u1*(u1^2+v1^2)*beta1        ];
F2  = [rho2*uf  ;rho2*uf^2+pf   ;rho2*uf*vf2  ;(alpha2/(gamma1-1)+(1-alpha2)/(gamma2-1)+1)*uf*pf+1/2*rho2*uf*(uf^2+vf2^2)       ;beta2*rho2*uf   ;gamma1/(gamma1-1)*uf*pf*alpha2+1/2*rho2*uf*(uf^2+vf2^2)*beta2       ];
F3  = [rho3*uf  ;rho3*uf^2+pf   ;rho3*uf*vf3  ;(alpha3/(gamma1-1)+(1-alpha3)/(gamma2-1)+1)*uf*pf+1/2*rho3*uf*(uf^2+vf3^2)       ;beta3*rho3*uf   ;gamma1/(gamma1-1)*uf*pf*alpha3+1/2*rho3*uf*(uf^2+vf3^2)*beta3       ];
Fr  = [rho4*u4  ;rho4*u4^2+p4   ;rho4*u4*v4   ;(alpha4/(gamma1-1)+(1-alpha4)/(gamma2-1)+1)*u4*p4+1/2*rho4*u4*(u4^2+v4^2)        ;beta4*rho4*u4   ;gamma1/(gamma1-1)*u4*p4*alpha4+1/2*rho4*u4*(u4^2+v4^2)*beta4        ];
Fsl = [rhosl*usl;rhosl*usl^2+psl;rhosl*usl*vsl;(alphasl/(gamma1-1)+(1-alphasl)/(gamma2-1)+1)*usl*psl+1/2*rhosl*usl*(usl^2+vsl^2);betasl*rhosl*usl;gamma1/(gamma1-1)*usl*psl*alphasl+1/2*rhosl*usl*(usl^2+vsl^2)*betasl];
Fsr = [rhosr*usr;rhosr*usr^2+psr;rhosr*usr*vsr;(alphasr/(gamma1-1)+(1-alphasr)/(gamma2-1)+1)*usr*psr+1/2*rhosr*usr*(usr^2+vsr^2);betasr*rhosr*usr;gamma1/(gamma1-1)*usr*psr*alphasr+1/2*rhosr*usr*(usr^2+vsr^2)*betasr];

lambda0 = uf;
lambda1 = u1-a1;
lambda2 = uf-a2;
lambda3 = uf+a3;
lambda4 = u4+a4;
% Lambda  = max(abs(lambda0),abs(lambda1),abs(lambda2),abs(lambda3),abs(lambda4),abs(Lambda)))

if lambda0>=0 && lambda2>=0
    if lambda1>=0 && lambda4>=0
        F(i,j,:) = Fl;
        sint12 = Sint_lin(rho1,u1,p1,beta1,alpha1,uf,-1);
        sintCD = SintCD(uf,pf,alpha3,alpha2);
        sint34 = Sint_lin(rho3,uf,pf,beta3,alpha3,u4,+1);
        sL(i,j) = 0;
        sR(i,j) = sint12+sintCD+sint34;
    elseif lambda1>=0 && lambda4<=0
        F(i,j,:) = Fl+Fr-Fsr;
        sintsr4 = Sint_lin(rhosr,usr,psr,betasr,alphasr,u4,+1);
        sint12  = Sint_lin(rho1,u1,p1,beta1,alpha1,uf,-1);
        sintCD  = SintCD(uf,pf,alpha3,alpha2);
        sint3sr = Sint_lin(rho3,uf,pf,beta3,alpha3,usr,+1);
        sL(i,j) = sintsr4;
        sR(i,j) = sint12+sintCD+sint3sr;
    elseif lambda1<=0 && lambda4>=0
        F(i,j,:) = Fsl;
        sint1sl = Sint_lin(rho1,u1,p1,beta1,alpha1,usl,-1);
        sintsl2 = Sint_lin(rhosl,usl,psl,betasl,alphasl,uf,-1);
        sintCD  = SintCD(uf,pf,alpha3,alpha2);
        sint34  = Sint_lin(rho3,uf,pf,beta3,alpha3,u4,+1);
        sL(i,j) = sint1sl;
        sR(i,j) = sintsl2+sintCD+sint34;
    elseif lambda1<=0 && lambda4<=0
        F(i,j,:) = Fsl-Fsr+Fr;
        sint1sl = Sint_lin(rho1,u1,p1,beta1,alpha1,usl,-1);
        sintsr4 = Sint_lin(rhosr,usr,psr,betasr,alphasr,u4,+1);
        sintsl2 = Sint_lin(rhosl,usl,psl,betasl,alphasl,uf,-1);
        sintCD  = SintCD(uf,pf,alpha3,alpha2);
        sint3sr = Sint_lin(rho3,uf,pf,beta3,alpha3,usr,+1);
        sL(i,j) = sint1sl+sintsr4;
        sR(i,j) = sintsl2+sintCD+sint3sr;
    end

elseif lambda0>=0 && lambda2<=0
    if lambda1>=0 && lambda4>=0
        F(i,j,:) = Fl-Fsl+F2;
        sintsl2 = Sint_lin(rhosl,usl,psl,betasl,alphasl,uf,-1);
        sint1sl = Sint_lin(rho1,u1,p1,beta1,alpha1,usl,-1);
        sintCD  = SintCD(uf,pf,alpha3,alpha2);
        sint34  = Sint_lin(rho3,uf,pf,beta3,alpha3,u4,+1);
        sL(i,j) = sintsl2;
        sR(i,j) = sint1sl+sintCD+sint34;
    elseif lambda1>=0 && lambda4<=0
        F(i,j,:) = Fl-Fsl+F2-Fsr+Fr;
        sintsl2 = Sint_lin(rhosl,usl,psl,betasl,alphasl,uf,-1);
        sintsr4 = Sint_lin(rhosr,usr,psr,betasr,alphasr,u4,+1);
        sint1sl = Sint_lin(rho1,u1,p1,beta1,alpha1,usl,-1);
        sintCD  = SintCD(uf,pf,alpha3,alpha2);
        sint3sr = Sint_lin(rho3,uf,pf,beta3,alpha3,usr,+1);
        sL(i,j) = sintsl2+sintsr4;
        sR(i,j) = sint1sl+sintCD+sint3sr;
    elseif lambda1<=0 && lambda4>=0
        F(i,j,:) = F2;
        sint12 = Sint_lin(rho1,u1,p1,beta1,alpha1,uf,-1);
        sintCD = SintCD(uf,pf,alpha3,alpha2);
        sint34 = Sint_lin(rho3,uf,pf,beta3,alpha3,u4,+1);
        sL(i,j) = sint12;
        sR(i,j) = sintCD+sint34;
    elseif lambda1<=0 && lambda4<=0
        F(i,j,:) = Fr+F2-Fsr;
        sint12  = Sint_lin(rho1,u1,p1,beta1,alpha1,uf,-1);
        sintsr4 = Sint_lin(rhosr,usr,psr,betasr,alphasr,u4,+1);
        sintCD  = SintCD(uf,pf,alpha3,alpha2);
        sint3sr = Sint_lin(rho3,uf,pf,beta3,alpha3,usr,+1);
        sL(i,j) = sint12+sintsr4;
        sR(i,j) = sintCD+sint3sr;
    end

elseif  (lambda0<=0 && lambda3>=0)
    if (lambda1>=0 && lambda4>=0)
        F(i,j,:) = Fl-Fsl+F3;
        sintsl2 = Sint_lin(rhosl,usl,psl,betasl,alphasl,uf,-1);
        sintCD  = SintCD(uf,pf,alpha3,alpha2);
        sint1sl = Sint_lin(rho1,u1,p1,beta1,alpha1,usl,-1);
        sint34  = Sint_lin(rho3,uf,pf,beta3,alpha3,u4,+1);
        sL(i,j) = sintsl2+sintCD;
        sR(i,j) = sint1sl+sint34;
    elseif (lambda1>=0 && lambda4<=0)
        F(i,j,:) = Fl-Fsl+F3-Fsr+Fr;
        sintsl2 = Sint_lin(rhosl,usl,psl,betasl,alphasl,uf,-1);
        sintCD  = SintCD(uf,pf,alpha3,alpha2);
        sintsr4 = Sint_lin(rhosr,usr,psr,betasr,alphasr,u4,+1);
        sint1sl = Sint_lin(rho1,u1,p1,beta1,alpha1,usl,-1);
        sint3sr = Sint_lin(rho3,uf,pf,beta3,alpha3,usr,+1);
        sL(i,j) = sintsl2+sintCD+sintsr4;
        sR(i,j) = sint1sl+sint3sr;
    elseif (lambda1<=0 && lambda4>=0)
        F(i,j,:) = F3;
        sint12 = Sint_lin(rho1,u1,p1,beta1,alpha1,uf,-1);
        sintCD = SintCD(uf,pf,alpha3,alpha2);
        sint34 = Sint_lin(rho3,uf,pf,beta3,alpha3,u4,+1);
        sL(i,j) = sint12+sintCD;
        sR(i,j) = sint34;
    elseif (lambda1<=0 && lambda4<=0)
        F(i,j,:) = Fr+F3-Fsr;
        sint12  = Sint_lin(rho1,u1,p1,beta1,alpha1,uf,-1);
        sintCD  = SintCD(uf,pf,alpha3,alpha2);
        sintsr4 = Sint_lin(rhosr,usr,psr,betasr,alphasr,u4,+1);
        sint3sr = Sint_lin(rho3,uf,pf,beta3,alpha3,usr,+1);
        sL(i,j) = sint12+sintCD+sintsr4;
        sR(i,j) = sint3sr;
    end

elseif  (lambda0<=0 && lambda3<=0)
    if (lambda1>=0 && lambda4>=0)
        F(i,j,:) = Fl-Fsl+Fsr;
        sintsl2 = Sint_lin(rhosl,usl,psl,betasl,alphasl,uf,-1);
        sintCD  = SintCD(uf,pf,alpha3,alpha2);
        sint3sr = Sint_lin(rho3,uf,pf,beta3,alpha3,usr,+1);
        sint1sl = Sint_lin(rho1,u1,p1,beta1,alpha1,usl,-1);
        sintsr4 = Sint_lin(rhosr,usr,psr,betasr,alphasr,u4,+1);
        sL(i,j) = sintsl2+sintCD+sint3sr;
        sR(i,j) = sint1sl+sintsr4;
    elseif (lambda1>=0 && lambda4<=0)
        F(i,j,:) = Fl-Fsl+Fr;
        sintsl2 = Sint_lin(rhosl,usl,psl,betasl,alphasl,uf,-1);
        sintCD  = SintCD(uf,pf,alpha3,alpha2);
        sint34  = Sint_lin(rho3,uf,pf,beta3,alpha3,u4,+1);
        sint1sl = Sint_lin(rho1,u1,p1,beta1,alpha1,usl,-1);
        sL(i,j) = sintsl2+sintCD+sint34;
        sR(i,j) = sint1sl;
    elseif (lambda1<=0 && lambda4>=0)
        F(i,j,:) = Fsr;
        sint12  = Sint_lin(rho1,u1,p1,beta1,alpha1,uf,-1);
        sintCD  = SintCD(uf,pf,alpha3,alpha2);
        sint3sr = Sint_lin(rho3,uf,pf,beta3,alpha3,usr,+1);
        sintsr4 = Sint_lin(rhosr,usr,psr,betasr,alphasr,u4,+1);
        sL(i,j) = sint12+sintCD+sint3sr;
        sR(i,j) = sintsr4;
    elseif (lambda1<=0 && lambda4<=0)
        F(i,j,:) = Fr;
        sint12 = Sint_lin(rho1,u1,p1,beta1,alpha1,uf,-1);
        sintCD = SintCD(uf,pf,alpha3,alpha2);
        sint34 = Sint_lin(rho3,uf,pf,beta3,alpha3,u4,+1);
        sL(i,j) = sint12+sintCD+sint34;
        sR(i,j) = 0;
    end

end %flux choice

end %if L==R

end %for j loop
end %for i loop


end % linosherflux_x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rhosl,usl,vsl,psl,betasl,alphasl]=linosher_sonic_L_x(rho1,u1,v1,p1,beta1,alpha1,a1)

global gamma1 gamma2

%------------------------------------%
%                                    %
% Solver for left sonic point        %
%                                    %
%------------------------------------%

uguess = a1;

%-----------------------------------------------------------------------%
%                         runga-kutta fourth-order                      %
%-----------------------------------------------------------------------%
diffu = 1;
m = 0;

while diffu>1e-8
m=m+1;

%-----------------------------------------------------------------------%
%                               LEFT WAVE                               %
%-----------------------------------------------------------------------%
if m==1
usl = uguess;
elseif m==2
usl = 0.8*uguess;
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



function [rhosr,usr,vsr,psr,betasr,alphasr]=linosher_sonic_R_x(rho4,u4,v4,p4,beta4,alpha4,a4)

global gamma1 gamma2

%Runga-Kutta Solver
uguess = -a4;


%-----------------------------------------------------------------------%
%                       Runga-Kutta fourth order                        %
%-----------------------------------------------------------------------%
diffu = 1;
m = 0;

while diffu>1e-8
m=m+1;
if m>100; keyboard; end
%-----------------------------------------------------------------------%
%                               RIGHT WAVE                              %
%-----------------------------------------------------------------------%
if m==1
usr = uguess;
elseif m==2
usr = 0.8*uguess;
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



function sint=Sint_lin(rho,u,p,beta,alpha,u_end,lr)

global gamma1 gamma2

%--------------------------------------------------%
%                                                  %
%    Source in isentropic wave                     %
%                                                  %
%--------------------------------------------------%

% lr is a sign function for left running wave (-) or right running wave (+)

du   = (u_end - u);
a    = sqrt(1/(alpha/gamma1+(1-alpha)/gamma2)*p/rho);
bS   = alpha*(1-alpha)*(gamma1-gamma2)/((1-alpha)*gamma1+alpha*gamma2);
sint = (bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u)*du;

end % Sint_lin


function sintCD=SintCD(uf,pf,alpha3,alpha2)

%--------------------------------------------------%
%                                                  %
%    Source over contact discontinuity             %
%                                                  %
%--------------------------------------------------%

sintCD = uf*pf*(alpha3-alpha2);

end % sub_sintCD_x