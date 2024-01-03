function [G,sA,sB]=linosherflux_y(G,Wfy)

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

sA = zeros(NY+1,NX);
sB = zeros(NY+1,NX);

for i=2:NY
for j=1:NX

rho1   = Wfy(2*(i-1),j,3);
u1     = Wfy(2*(i-1),j,4);
v1     = Wfy(2*(i-1),j,5);
p1     = Wfy(2*(i-1),j,6);
beta1  = Wfy(2*(i-1),j,7);
alpha1 = Wfy(2*(i-1),j,8);
a1     = sqrt(1/(alpha1/gamma1+(1-alpha1)/gamma2)*p1/rho1);

rho4   = Wfy(2*i-1,j,3);
u4     = Wfy(2*i-1,j,4);
v4     = Wfy(2*i-1,j,5);
p4     = Wfy(2*i-1,j,6);
beta4  = Wfy(2*i-1,j,7);
alpha4 = Wfy(2*i-1,j,8);
a4     = sqrt(1/(alpha4/gamma1+(1-alpha4)/gamma2)*p4/rho4);


if  abs(rho1-rho4)+abs(u1-u4)+abs(v1-v4)+abs(p1-p4)<1e-8

G(i,j,1) = rho1*v1;
G(i,j,2) = rho1*u1*v1;
G(i,j,3) = rho1*v1^2+p1;
G(i,j,4) = (alpha1/(gamma1-1)+(1-alpha1)/(gamma2-1)+1)*v1*p1+1/2*rho1*v1*(u1^2+v1^2);
G(i,j,5) = beta1*rho1*v1;
G(i,j,6) = gamma1/(gamma1-1)*v1*p1*alpha1+1/2*rho1*v1*(u1^2+v1^2)*beta1;

elseif ( abs(p1-p4) + abs(v1) + abs(v4) ) <1e-8

G(i,j,3) = p1;

else

% Final state values
[rho2,rho3,uf2,uf3,vf,pf,beta2,beta3,alpha2,alpha3]=linoshersolver_y(rho1,u1,v1,p1,beta1,alpha1,a1,rho4,u4,v4,p4,beta4,alpha4,a4);

% speed of sound
a2 = sqrt(1/(alpha2/gamma1+(1-alpha2)/gamma2)*pf/rho2);
a3 = sqrt(1/(alpha3/gamma1+(1-alpha3)/gamma2)*pf/rho3);

% Sonic point values
[rhosl,usl,vsl,psl,betasl,alphasl] = linosher_sonic_L_y(rho1,u1,v1,p1,beta1,alpha1,a1);
[rhosr,usr,vsr,psr,betasr,alphasr] = linosher_sonic_R_y(rho4,u4,v4,p4,beta4,alpha4,a4);

Gl(1) = rho1*v1;
Gl(2) = rho1*v1*u1;
Gl(3) = rho1*v1^2+p1;
Gl(4) = (alpha1/(gamma1-1)+(1-alpha1)/(gamma2-1)+1)*v1*p1+1/2*rho1*v1*(u1^2+v1^2);
Gl(5) = beta1*rho1*v1;
Gl(6) = gamma1/(gamma1-1)*v1*p1*alpha1+1/2*rho1*v1*(u1^2+v1^2)*beta1;

G2(1) = rho2*vf;
G2(2) = rho2*uf2*vf;
G2(3) = rho2*vf^2+pf;
G2(4) = (alpha2/(gamma1-1)+(1-alpha2)/(gamma2-1)+1)*vf*pf+1/2*rho2*vf*(uf2^2+vf^2);
G2(5) = beta2*rho2*vf;
G2(6) = gamma1/(gamma1-1)*vf*pf*alpha2+1/2*rho2*vf*(uf2^2+vf^2)*beta2;

G3(1) = rho3*vf;
G3(2) = rho3*uf3*vf;
G3(3) = rho3*vf^2+pf;
G3(4) = (alpha3/(gamma1-1)+(1-alpha3)/(gamma2-1)+1)*vf*pf+1/2*rho3*vf*(uf3^2+vf^2);
G3(5) = beta3*rho3*vf;
G3(6) = gamma1/(gamma1-1)*vf*pf*alpha3+1/2*rho3*vf*(uf3^2+vf^2)*beta3;

Gr(1) = rho4*v4;
Gr(2) = rho4*u4*v4;
Gr(3) = rho4*v4^2+p4;
Gr(4) = (alpha4/(gamma1-1)+(1-alpha4)/(gamma2-1)+1)*v4*p4+1/2*rho4*v4*(u4^2+v4^2);
Gr(5) = beta4*rho4*v4;
Gr(6) = gamma1/(gamma1-1)*v4*p4*alpha4+1/2*rho4*v4*(u4^2+v4^2)*beta4;

Gsl(1) = rhosl*vsl;
Gsl(2) = rhosl*usl*vsl;
Gsl(3) = rhosl*vsl^2+psl;
Gsl(4) = (alphasl/(gamma1-1)+(1-alphasl)/(gamma2-1)+1)*vsl*psl+1/2*rhosl*vsl*(usl^2+vsl^2);
Gsl(5) = betasl*rhosl*vsl;
Gsl(6) = gamma1/(gamma1-1)*vsl*psl*alphasl+1/2*rhosl*vsl*(usl^2+vsl^2)*betasl;

Gsr(1) = rhosr*vsr;
Gsr(2) = rhosr*usr*vsr;
Gsr(3) = rhosr*vsr^2+psr;
Gsr(4) = (alphasr/(gamma1-1)+(1-alphasr)/(gamma2-1)+1)*vsr*psr+1/2*rhosr*vsr*(usr^2+vsr^2);
Gsr(5) = betasr*rhosr*vsr;
Gsr(6) = gamma1/(gamma1-1)*vsr*psr*alphasr+1/2*rhosr*vsr*(usr^2+vsr^2)*betasr;

lambda0 = vf;
lambda1 = v1-a1;
lambda2 = vf-a2;
lambda3 = vf+a3;
lambda4 = v4+a4;
% Lambda  = max(abslambda0),abslambda1),abslambda2),abslambda3),abslambda4),abs(Lambda)));


if lambda0>=0 && lambda2>=0
    if lambda1>=0 && lambda4>=0
        G(i,j,:) = Gl;
        sint12 = sub_Sint_lin_y(rho1,v1,p1,beta1,alpha1,vf,-1);
        sintCD = sub_SintCD_y(vf,pf,alpha3,alpha2);
        sint34 = sub_Sint_lin_y(rho3,vf,pf,beta3,alpha3,v4,+1);
        sA(i,j) = 0;
        sB(i,j) = sint12+sintCD+sint34;
    elseif lambda1>=0 && lambda4<=0
        G(i,j,:) = Gl+Gr-Gsr;
        sintsr4 = sub_Sint_lin_y(rhosr,vsr,psr,betasr,alphasr,v4,+1);
        sint12  = sub_Sint_lin_y(rho1,v1,p1,beta1,alpha1,vf,-1);
        sintCD  = sub_SintCD_y(vf,pf,alpha3,alpha2);
        sint3sr = sub_Sint_lin_y(rho3,vf,pf,beta3,alpha3,vsr,+1);
        sA(i,j) = sintsr4;
        sB(i,j) = sint12+sintCD+sint3sr;
    elseif lambda1<=0 && lambda4>=0
        G(i,j,:) = Gsl;
        sint1sl = sub_Sint_lin_y(rho1,v1,p1,beta1,alpha1,vsl,-1);
        sintsl2 = sub_Sint_lin_y(rhosl,vsl,psl,betasl,alphasl,vf,-1);
        sintCD  = sub_SintCD_y(vf,pf,alpha3,alpha2);
        sint34  = sub_Sint_lin_y(rho3,vf,pf,beta3,alpha3,v4,+1);
        sA(i,j) = sint1sl;
        sB(i,j) = sintsl2+sintCD+sint34;
    elseif lambda1<=0 && lambda4<=0
        G(i,j,:) = Gsl-Gsr+Gr;
        sint1sl = sub_Sint_lin_y(rho1,v1,p1,beta1,alpha1,vsl,-1);
        sintsr4 = sub_Sint_lin_y(rhosr,vsr,psr,betasr,alphasr,v4,+1);
        sintsl2 = sub_Sint_lin_y(rhosl,vsl,psl,betasl,alphasl,vf,-1);
        sintCD  = sub_SintCD_y(vf,pf,alpha3,alpha2);
        sint3sr = sub_Sint_lin_y(rho3,vf,pf,beta3,alpha3,vsr,+1);
        sA(i,j) = sint1sl+sintsr4;
        sB(i,j) = sintsl2+sintCD+sint3sr;
    end

elseif  lambda0>=0 && lambda2<=0
    if lambda1>=0 && lambda4>=0
        G(i,j,:) = Gl-Gsl+G2;
        sintsl2 = sub_Sint_lin_y(rhosl,vsl,psl,betasl,alphasl,vf,-1);
        sint1sl = sub_Sint_lin_y(rho1,v1,p1,beta1,alpha1,vsl,-1);
        sintCD  = sub_SintCD_y(vf,pf,alpha3,alpha2);
        sint34  = sub_Sint_lin_y(rho3,vf,pf,beta3,alpha3,v4,+1);
        sA(i,j) = sintsl2;
        sB(i,j) = sint1sl+sintCD+sint34;
    elseif lambda1>=0 && lambda4<=0
        G(i,j,:) = Gl-Gsl+G2-Gsr+Gr;
        sintsl2 = sub_Sint_lin_y(rhosl,vsl,psl,betasl,alphasl,vf,-1);
        sintsr4 = sub_Sint_lin_y(rhosr,vsr,psr,betasr,alphasr,v4,+1);
        sint1sl = sub_Sint_lin_y(rho1,v1,p1,beta1,alpha1,vsl,-1);
        sintCD  = sub_SintCD_y(vf,pf,alpha3,alpha2);
        sint3sr = sub_Sint_lin_y(rho3,vf,pf,beta3,alpha3,vsr,-1);
        sA(i,j) = sintsl2+sintsr4;
        sB(i,j) = sint1sl+sintCD+sint3sr;
    elseif lambda1<=0 && lambda4>=0
        G(i,j,:) = G2;
        sint12 = sub_Sint_lin_y(rho1,v1,p1,beta1,alpha1,vf,-1);
        sintCD = sub_SintCD_y(vf,pf,alpha3,alpha2);
        sint34 = sub_Sint_lin_y(rho3,vf,pf,beta3,alpha3,v4,+1);
        sA(i,j) = sint12;
        sB(i,j) = sintCD+sint34;
    elseif lambda1<=0 && lambda4<=0
        G(i,j,:) = Gr+G2-Gsr;
        sint12  = sub_Sint_lin_y(rho1,v1,p1,beta1,alpha1,vf,-1);
        sintsr4 = sub_Sint_lin_y(rhosr,vsr,psr,betasr,alphasr,v4,+1);
        sintCD  = sub_SintCD_y(vf,pf,alpha3,alpha2);
        sint3sr = sub_Sint_lin_y(rho3,vf,pf,beta3,alpha3,vsr,+1);
        sA(i,j) = sint12+sintsr4;
        sB(i,j) = sintCD+sint3sr;
    end

elseif  lambda0<=0 && lambda3>=0
    if lambda1>=0 && lambda4>=0
        G(i,j,:) = Gl-Gsl+G3;
        sintsl2 = sub_Sint_lin_y(rhosl,vsl,psl,betasl,alphasl,vf,-1);
        sintCD  = sub_SintCD_y(vf,pf,alpha3,alpha2);
        sint1sl = sub_Sint_lin_y(rho1,v1,p1,beta1,alpha1,vsl,-1);
        sint34  = sub_Sint_lin_y(rho3,vf,pf,beta3,alpha3,v4,+1);
        sA(i,j) = sintsl2+sintCD;
        sB(i,j) = sint1sl+sint34;
    elseif lambda1>=0 && lambda4<=0
        G(i,j,:) = Gl-Gsl+G3-Gsr+Gr;
        sintsl2 = sub_Sint_lin_y(rhosl,vsl,psl,betasl,alphasl,vf,-1);
        sintCD  = sub_SintCD_y(vf,pf,alpha3,alpha2);
        sintsr4 = sub_Sint_lin_y(rhosr,vsr,psr,betasr,alphasr,v4,+1);
        sint1sl = sub_Sint_lin_y(rho1,v1,p1,beta1,alpha1,vsl,-1);
        sint3sr = sub_Sint_lin_y(rho3,vf,pf,beta3,alpha3,vsr,+1);
        sA(i,j) = sintsl2+sintCD+sintsr4;
        sB(i,j) = sint1sl+sint3sr;
    elseif lambda1<=0 && lambda4>=0
        G(i,j,:) = G3;
        sint12 = sub_Sint_lin_y(rho1,v1,p1,beta1,alpha1,vf,-1);
        sintCD = sub_SintCD_y(vf,pf,alpha3,alpha2);
        sint34 = sub_Sint_lin_y(rho3,vf,pf,beta3,alpha3,v4,+1);
        sA(i,j) = sint12+sintCD;
        sB(i,j) = sint34;
    elseif lambda1<=0 && lambda4<=0
        G(i,j,:) = Gr+G3-Gsr;
        sint12  = sub_Sint_lin_y(rho1,v1,p1,beta1,alpha1,vf,-1);
        sintCD  = sub_SintCD_y(vf,pf,alpha3,alpha2);
        sintsr4 = sub_Sint_lin_y(rhosr,vsr,psr,betasr,alphasr,v4,+1);
        sint3sr = sub_Sint_lin_y(rho3,vf,pf,beta3,alpha3,vsr,+1);
        sA(i,j) = sint12+sintCD+sintsr4;
        sB(i,j) = sint3sr;
    end

elseif  lambda0<=0 && lambda3<=0
    if lambda1>=0 && lambda4>=0
        G(i,j,:) = Gl-Gsl+Gsr;
        sintsl2 = sub_Sint_lin_y(rhosl,vsl,psl,betasl,alphasl,vf,-1);
        sintCD  = sub_SintCD_y(vf,pf,alpha3,alpha2);
        sint3sr = sub_Sint_lin_y(rho3,vf,pf,beta3,alpha3,vsr,+1);
        sint1sl = sub_Sint_lin_y(rho1,v1,p1,beta1,alpha1,vsl,-1);
        sintsr4 = sub_Sint_lin_y(rhosr,vsr,psr,betasr,alphasr,v4,+1);
        sA(i,j) = sintsl2+sintCD+sint3sr;
        sB(i,j) = sint1sl+sintsr4;
    elseif lambda1>=0 && lambda4<=0
        G(i,j,:) = Gl-Gsl+Gr;
        sintsl2 = sub_Sint_lin_y(rhosl,vsl,psl,betasl,alphasl,vf,-1);
        sintCD  = sub_SintCD_y(vf,pf,alpha3,alpha2);
        sint34  = sub_Sint_lin_y(rho3,vf,pf,beta3,alpha3,v4,+1);
        sint1sl = sub_Sint_lin_y(rho1,v1,p1,beta1,alpha1,vsl,-1);
        sA(i,j) = sintsl2+sintCD+sint34;
        sB(i,j) = sint1sl;
    elseif lambda1<=0 && lambda4>=0
        G(i,j,:) = Gsr;
        sint12  = sub_Sint_lin_y(rho1,v1,p1,beta1,alpha1,vf,-1);
        sintCD  = sub_SintCD_y(vf,pf,alpha3,alpha2);
        sint3sr = sub_Sint_lin_y(rho3,vf,pf,beta3,alpha3,vsr,+1);
        sintsr4 = sub_Sint_lin_y(rhosr,vsr,psr,betasr,alphasr,v4,+1);
        sA(i,j) = sint12+sintCD+sint3sr;
        sB(i,j) = sintsr4;
    elseif lambda1<=0 && lambda4<=0
        G(i,j,:) = Gr;
        sint12 = sub_Sint_lin_y(rho1,v1,p1,beta1,alpha1,vf,-1);
        sintCD = sub_SintCD_y(vf,pf,alpha3,alpha2);
        sint34 = sub_Sint_lin_y(rho3,vf,pf,beta3,alpha3,v4,+1);
        sA(i,j) = sint12+sintCD+sint34;
        sB(i,j) = 0;
    end;

end %flux choice

end %if L==R

end %for j loop
end %for i loop

end % linosherflux_y


function [rho2,rho3,uf2,uf3,vf,pf,beta2,beta3,alpha2,alpha3]=linoshersolver_y(rho1,u1,v1,p1,beta1,alpha1,a1,rho4,u4,v4,p4,beta4,alpha4,a4)
%-----------------------------------%
%                                   %
% flux solver                       %
%                                   %
%-----------------------------------%

global gamma1 gamma2

% Fina1 state va1ues
vf     = (rho1*a1*v1+rho4*a4*v4+(p1-p4))/(rho1*a1+rho4*a4);
pf     = (rho4*a4*p1+rho1*a1*p4+rho1*a1*rho4*a4*(v1-v4))/(rho1*a1+rho4*a4);
rho2   = rho1+(pf-p1)/a1^2;
rho3   = rho4+(pf-p4)/a4^2;
uf2    = u1;
uf3    = u4;
beta2  = beta1;
beta3  = beta4;
alpha2 = alpha1+(pf-p1)*(1/gamma2-1/gamma1)*alpha1*(1-alpha1)/p1;
alpha3 = alpha4+(pf-p4)*(1/gamma2-1/gamma1)*alpha4*(1-alpha4)/p4;

end % linoshersolver_y


function [rhosl,usl,vsl,psl,betasl,alphasl]=linosher_sonic_L_y(rho1,u1,v1,p1,beta1,alpha1,a1)

global gamma1 gamma2

%------------------------------------%
%                                    %
% Solver for left sonic point        %
%                                    %
%------------------------------------%

vguess = a1;

diffv = 1.0;
m = 0;

while diffv>1e-8
    m=m+1;

    %-----------------------------------------------------------------------%
    %                               LEFT WAVE                               %
    %-----------------------------------------------------------------------%
    if (m==1)
        vsl = vguess;
    elseif (m==2)
        vsl = 0.8*vguess;
    end

    psl     = p1-rho1*a1*(vsl-v1);
    rhosl   = rho1-rho1/a1*(vsl-v1);
    b1      = alpha1*(1-alpha1)*(gamma1-gamma2)/((1-alpha1)*gamma1+alpha1*gamma2);
    alphasl = alpha1-b1/a1*(vsl-v1);
    asl     = sqrt(1/(alphasl/gamma1+(1-alphasl)/gamma2)*psl/rhosl);

    %------------------------------------------------------------------------

    dva     = vsl - asl;
    diffv   = abs(dva);

    if (m==1)
        vslold = vsl;
        dvaold = dva;
    else
        vslnew = vsl - dva*(vsl-vslold+eps)/(dva-dvaold+eps);
        vslold = vsl;
        vsl    = vslnew;
        dvaold = dva;
    end

end %while loop

psl     = p1-rho1*a1*(vsl-v1);
rhosl   = rho1-rho1/a1*(vsl-v1);
usl     = u1;
betasl  = beta1;
b1      = alpha1*(1-alpha1)*(gamma1-gamma2)/((1-alpha1)*gamma1+alpha1*gamma2);
alphasl = alpha1-b1/a1*(vsl-v1);

end % linosher_sonic_L_y


function [rhosr,usr,vsr,psr,betasr,alphasr]=linosher_sonic_R_y(rho4,u4,v4,p4,beta4,alpha4,a4)

global gamma1 gamma2

vguess = -a4;

diffv = 1.0;
m = 0;

while (diffv>1e-8)
    m=m+1;
if m>100; keyboard; end
    %-----------------------------------------------------------------------%
    %                               LEFT WAVE                               %
    %-----------------------------------------------------------------------%
    if (m==1)
        vsr = vguess;
    elseif (m==2)
        vsr = 0.8*vguess;
    end

    psr     = p4+rho4*a4*(vsr-v4);
    rhosr   = rho4+rho4/a4*(vsr-v4);
    b4      = alpha4*(1-alpha4)*(gamma1-gamma2)/((1-alpha4)*gamma1+alpha4*gamma2);
    alphasr = alpha4+b4/a4*(vsr-v4);
    asr     = sqrt(1/(alphasr/gamma1+(1-alphasr)/gamma2)*psr/rhosr);
    %-----------------------------------------------------------------------%

    dva     = vsr + asr;
    diffv   = abs(dva);

    if (m==1)
        vsrold = vsr;
        dvaold = dva;
    else
        vsrnew = vsr - dva*(vsr-vsrold+eps)/(dva-dvaold+eps);
        vsrold = vsr;
        vsr    = vsrnew;
        dvaold = dva;
    end

end %while loop

psr     = p4+rho4*a4*(vsr-v4);
rhosr   = rho4+rho4/a4*(vsr-v4);
usr     = u4;
betasr  = beta4;
b4      = alpha4*(1-alpha4)*(gamma1-gamma2)/((1-alpha4)*gamma1+alpha4*gamma2);
alphasr = alpha4+b4/a4*(vsr-v4);

end % linosher_sonic_R_y


function [sint]=sub_Sint_lin_y(rho,v,p,beta,alpha,v_end,lr)
%--------------------------------------------------%
%                                                  %
%    Source in isentropic wave                     %
%                                                  %
%--------------------------------------------------%

global gamma1 gamma2

% lr is a sign function for left running wave (-) or right running wave (+)

dv   = (v_end - v);
a    = sqrt(1/(alpha/gamma1+(1-alpha)/gamma2)*p/rho);
bS   = alpha*(1-alpha)*(gamma1-gamma2)/((1-alpha)*gamma1+alpha*gamma2);
sint = (bS*p*(1+lr*v/a)+lr*(alpha-beta)*rho*a*v)*dv;

end % sub_Sint_lin_y


function sintCD=sub_SintCD_y(vf,pf,alpha3,alpha2)
%--------------------------------------------------%
%                                                  %
%    Source over contact discontinuity             %
%                                                  %
%--------------------------------------------------%

sintCD = vf*pf*(alpha3-alpha2);

end % sub_sintCD_y