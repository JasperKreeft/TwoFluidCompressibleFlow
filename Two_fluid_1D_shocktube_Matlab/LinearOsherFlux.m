%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              %
% Written by: Jasper Kreeft  (2007)            %
% Updated by: Jasper Kreeft  (2015)            %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,sL,sR,Lambda] = LinearOsherFlux(Wf,Lambda)

global N
global gamma1 gamma2
global pi1 pi2

epss = 1e-12;

F  = zeros(5,N+1);
sL = zeros(1,N); % waarden zijn gedefinieerd voor de cel index en niet voor de flux index
sR = zeros(1,N);

% i = j - 1/2
for i  = 2:N

    rho1   = Wf(1,2*i-2);
    u1     = Wf(2,2*i-2);
    p1     = Wf(3,2*i-2);
    beta1  = Wf(4,2*i-2);
    alpha1 = Wf(5,2*i-2);
%     a1     = sqrt(1/(alpha1/gamma1+(1-alpha1)/gamma2)*p1/rho1);
    a1     = sqrt(1/(alpha1/(gamma1*(p1+pi1))+(1-alpha1)/(gamma2*(p1+pi2)))/rho1);

    rho4   = Wf(1,2*i-1);
    u4     = Wf(2,2*i-1);
    p4     = Wf(3,2*i-1);
    beta4  = Wf(4,2*i-1);
    alpha4 = Wf(5,2*i-1);
%     a4     = sqrt(1/(alpha4/gamma1+(1-alpha4)/gamma2)*p4/rho4);
    a4     = sqrt(1/(alpha4/(gamma1*(p4+pi1))+(1-alpha4)/(gamma2*(p4+pi2)))/rho4);

    % Check for negative densities and pressures
    if rho1<=0
        disp('rho1 negative')
        rho1 = epss;
    end

    if p1<=0
        disp('p1 negative')
        p1 = epss;
    end

    if rho4<=0
        disp('rho4 negative')
        rho4 = epss;
    end

    if p4<=0
        disp('p4 negative')
        p4 = epss;
    end

    % Equal state in both cells
    if abs(rho1-rho4)+abs(u1-u4)+abs(p1-p4)+abs(beta1-beta4)+abs(alpha1-alpha4)<1e-10

%        F(:,i) = [rho1*u1;
%                  rho1*u1^2+p1;
%                  (alpha1/(gamma1-1)+(1-alpha1)/(gamma2-1)+1)*u1*p1+0.5*rho1*u1^3;
%                  beta1*rho1*u1;
%                  gamma1/(gamma1-1)*u1*p1*alpha1+1/2*rho1*u1^3*beta1 ];

        F(:,i) = [rho1*u1;
                 rho1*u1^2+p1;
                 alpha1*gamma1*(p1+pi1)*u1/(gamma1-1)+(1-alpha1)*gamma2*(p1+pi2)*u1/(gamma2-1)+1/2*rho1*u1^3;
                 rho1*u1*beta1;
                 alpha1*gamma1*(p1+pi1)*u1/(gamma1-1)+1/2*rho1*u1^3*beta1];
             
        Lambda  = max(abs([u1 u1-a1 u4+a4 Lambda]));

    % non-moving interface
    elseif abs(p1-p4)+abs(u1)+abs(u4)<1e-10

        F(2,i) = p1;

    % Two-fluid problem
    else

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate Linear Osher's solution                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        %Final state values
        u23    = (rho1*a1*u1+rho4*a4*u4+(p1-p4))/(rho1*a1+rho4*a4);
        p23    = (rho4*a4*p1+rho1*a1*p4+rho1*a1*rho4*a4*(u1-u4))/(rho1*a1+rho4*a4);
        rho2   = rho1+(p23-p1)/a1^2;
        rho3   = rho4+(p23-p4)/a4^2;
        beta2  = beta1;
        beta3  = beta4;
%         alpha2 = alpha1+(p23-p1)*(1/gamma2-1/gamma1)*alpha1*(1-alpha1)/p1;
%         alpha3 = alpha4+(p23-p4)*(1/gamma2-1/gamma1)*alpha4*(1-alpha4)/p4;
        alpha2 = alpha1+(p23-p1)*alpha1*(1-alpha1)*(1/(gamma2*(p1+pi1))-1/(gamma1*(p1+pi2)));
        alpha3 = alpha4+(p23-p4)*alpha4*(1-alpha4)*(1/(gamma2*(p4+pi1))-1/(gamma1*(p4+pi2)));

        %speed of sound
%         a2 = sqrt(1/(alpha2/gamma1+(1-alpha2)/gamma2)*p23/rho2);
%         a3 = sqrt(1/(alpha3/gamma1+(1-alpha3)/gamma2)*p23/rho3);
        a2 = sqrt(1/(alpha2/(gamma1*(p23+pi1))+(1-alpha2)/(gamma2*(p23+pi2)))/rho2);
        a3 = sqrt(1/(alpha3/(gamma1*(p23+pi1))+(1-alpha3)/(gamma2*(p23+pi2)))/rho3);


        %Sonic point va1ues
        [rhosl,usl,psl,betasl,alphasl]=linosher_sonic_L(rho1,u1,p1,beta1,alpha1,a1);
        [rhosr,usr,psr,betasr,alphasr]=linosher_sonic_R(rho4,u4,p4,beta4,alpha4,a4);

%         F1 = [rho1*u1; rho1*u1^2+p1; (alpha1/(gamma1-1)+(1-alpha1)/(gamma2-1)+1)*u1*p1+0.5*rho1*u1^3; beta1*rho1*u1; gamma1/(gamma1-1)*u1*p1*alpha1+1/2*rho1*u1^3*beta1];
%         F2 = [rho2*u23; rho2*u23^2+p23; (alpha2/(gamma1-1)+(1-alpha2)/(gamma2-1)+1)*u23*p23+0.5*rho2*u23^3; beta2*rho2*u23; gamma1/(gamma1-1)*u23*p23*alpha2+1/2*rho2*u23^3*beta2];
%         F3 = [rho3*u23; rho3*u23^2+p23; (alpha3/(gamma1-1)+(1-alpha3)/(gamma2-1)+1)*u23*p23+0.5*rho3*u23^3; beta3*rho3*u23; gamma1/(gamma1-1)*u23*p23*alpha3+1/2*rho3*u23^3*beta3];
%         F4 = [rho4*u4; rho4*u4^2+p4; (alpha4/(gamma1-1)+(1-alpha4)/(gamma2-1)+1)*u4*p4+0.5*rho4*u4^3; beta4*rho4*u4; gamma1/(gamma1-1)*u4*p4*alpha4+1/2*rho4*u4^3*beta4];
%         Fsl = [rhosl*usl; rhosl*usl^2+psl; (alphasl/(gamma1-1)+(1-alphasl)/(gamma2-1)+1)*usl*psl+0.5*rhosl*usl^3; betasl*rhosl*usl; gamma1/(gamma1-1)*usl*psl*alphasl+1/2*rhosl*usl^3*betasl];
%         Fsr = [rhosr*usr; rhosr*usr^2+psr; (alphasr/(gamma1-1)+(1-alphasr)/(gamma2-1)+1)*usr*psr+0.5*rhosr*usr^3; betasr*rhosr*usr; gamma1/(gamma1-1)*usr*psr*alphasr+1/2*rhosr*usr^3*betasr];
        F1 = [rho1*u1; rho1*u1^2+p1; alpha1*gamma1*(p1+pi1)*u1/(gamma1-1)+(1-alpha1)*gamma2*(p1+pi2)*u1/(gamma2-1)+0.5*rho1*u1^3; beta1*rho1*u1; gamma1/(gamma1-1)*u1*(p1+pi1)*alpha1+1/2*rho1*u1^3*beta1];
        F2 = [rho2*u23; rho2*u23^2+p23; alpha2*gamma1*(p23+pi1)*u23/(gamma1-1)+(1-alpha2)*gamma2*(p23+pi2)*u23/(gamma2-1)+0.5*rho2*u23^3; beta2*rho2*u23; gamma1/(gamma1-1)*u23*(p23+pi1)*alpha2+1/2*rho2*u23^3*beta2];
        F3 = [rho3*u23; rho3*u23^2+p23; alpha3*gamma1*(p23+pi1)*u23/(gamma1-1)+(1-alpha3)*gamma2*(p23+pi2)*u23/(gamma2-1)+0.5*rho3*u23^3; beta3*rho3*u23; gamma1/(gamma1-1)*u23*(p23+pi1)*alpha3+1/2*rho3*u23^3*beta3];
        F4 = [rho4*u4; rho4*u4^2+p4; alphasl*gamma1*(psl+pi1)*usl/(gamma1-1)+(1-alphasl)*gamma2*(psl+pi2)*usl/(gamma2-1)+0.5*rho4*u4^3; beta4*rho4*u4; gamma1/(gamma1-1)*u4*(p4+pi1)*alpha4+1/2*rho4*u4^3*beta4];
        Fsl = [rhosl*usl; rhosl*usl^2+psl; alphasl*gamma1*(psl+pi1)*usl/(gamma1-1)+(1-alphasl)*gamma2*(psl+pi2)*usl/(gamma2-1)+0.5*rhosl*usl^3; betasl*rhosl*usl; gamma1/(gamma1-1)*usl*(psl+pi1)*alphasl+1/2*rhosl*usl^3*betasl];
        Fsr = [rhosr*usr; rhosr*usr^2+psr; alphasr*gamma1*(psr+pi1)*usr/(gamma1-1)+(1-alphasr)*gamma2*(psr+pi2)*usr/(gamma2-1)+0.5*rhosr*usr^3; betasr*rhosr*usr; gamma1/(gamma1-1)*usr*(psr+pi1)*alphasr+1/2*rhosr*usr^3*betasr];
        
        lambda0 = u23;
        lambda1 = u1-a1;
        lambda2 = u23-a2;
        lambda3 = u23+a3;
        lambda4 = u4+a4;
        Lambda  = max(abs([lambda0 lambda1 lambda2 lambda3 lambda4 Lambda]));

        if lambda0 >= 0 && lambda2>=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = F1;
                sint12 = Sint_lin(rho1,u1,p1,beta1,alpha1,u23,-1);
                sintCD = SintCD(u23,p23,alpha3,alpha2);
                sint34 = Sint_lin(rho3,u23,p23,beta3,alpha3,u4,+1);
                sL(i-1) = 0;
                sR(i)   = sint12+sintCD+sint34;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = F1+F4-Fsr;
                sintsr4 = Sint_lin(rhosr,usr,psr,betasr,alphasr,u4,+1);
                sint12  = Sint_lin(rho1,u1,p1,beta1,alpha1,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint3sr = Sint_lin(rho3,u23,p23,beta3,alpha3,usr,+1);
                sL(i-1) = sintsr4;
                sR(i)   = sint12+sintCD+sint3sr;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = Fsl;
                sint1sl = Sint_lin(rho1,u1,p1,beta1,alpha1,usl,-1);
                sintsl2 = Sint_lin(rhosl,usl,psl,betasl,alphasl,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint34  = Sint_lin(rho3,u23,p23,beta3,alpha3,u4,+1);
                sL(i-1) = sint1sl;
                sR(i)   = sintsl2+sintCD+sint34;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = Fsl-Fsr+F4;
                sint1sl = Sint_lin(rho1,u1,p1,beta1,alpha1,usl,-1);
                sintsr4 = Sint_lin(rhosr,usr,psr,betasr,alphasr,u4,+1);
                sintsl2 = Sint_lin(rhosl,usl,psl,betasl,alphasl,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint3sr = Sint_lin(rho3,u23,p23,beta3,alpha3,usr,+1);
                sL(i-1) = sint1sl+sintsr4;
                sR(i)   = sintsl2+sintCD+sint3sr;
            end

        elseif  lambda0 >= 0 && lambda2<=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = F1-Fsl+F2;
                sintsl2 =Sint_lin(rhosl,usl,psl,betasl,alphasl,u23,-1);
                sint1sl =Sint_lin(rho1,u1,p1,beta1,alpha1,usl,-1);
                sintCD  =SintCD(u23,p23,alpha3,alpha2);
                sint34  =Sint_lin(rho3,u23,p23,beta3,alpha3,u4,+1);
                sL(i-1) = sintsl2;
                sR(i)   = sint1sl+sintCD+sint34;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = F1-Fsl+F2-Fsr+F4;
                sintsl2 = Sint_lin(rhosl,usl,psl,betasl,alphasl,u23,-1);
                sintsr4 = Sint_lin(rhosr,usr,psr,betasr,alphasr,u4,+1);
                sint1sl = Sint_lin(rho1,u1,p1,beta1,alpha1,usl,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint3sr = Sint_lin(rho3,u23,p23,beta3,alpha3,usr,+1);
                sL(i-1) = sintsl2+sintsr4;
                sR(i)   = sint1sl+sintCD+sint3sr;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = F2;
                sint12 = Sint_lin(rho1,u1,p1,beta1,alpha1,u23,-1);
                sintCD = SintCD(u23,p23,alpha3,alpha2);
                sint34 = Sint_lin(rho3,u23,p23,beta3,alpha3,u4,+1);
                sL(i-1) = sint12;
                sR(i)   = sintCD+sint34;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = F4+F2-Fsr;
                sint12  = Sint_lin(rho1,u1,p1,beta1,alpha1,u23,-1);
                sintsr4 = Sint_lin(rhosr,usr,psr,betasr,alphasr,u4,+1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint3sr = Sint_lin(rho3,u23,p23,beta3,alpha3,usr,+1);
                sL(i-1) = sint12+sintsr4;
                sR(i)   = sintCD+sint3sr;
            end

        elseif  lambda0 <= 0 && lambda3>=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = F1-Fsl+F3;
                sintsl2 = Sint_lin(rhosl,usl,psl,betasl,alphasl,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint1sl = Sint_lin(rho1,u1,p1,beta1,alpha1,usl,-1);
                sint34  = Sint_lin(rho3,u23,p23,beta3,alpha3,u4,+1);
                sL(i-1) = sintsl2+sintCD;
                sR(i)   = sint1sl+sint34;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = F1-Fsl+F3-Fsr+F4;
                sintsl2 = Sint_lin(rhosl,usl,psl,betasl,alphasl,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sintsr4 = Sint_lin(rhosr,usr,psr,betasr,alphasr,u4,+1);
                sint1sl = Sint_lin(rho1,u1,p1,beta1,alpha1,usl,-1);
                sint3sr = Sint_lin(rho3,u23,p23,beta3,alpha3,usr,+1);
                sL(i-1) = sintsl2+sintCD+sintsr4;
                sR(i)   = sint1sl+sint3sr;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = F3;
                sint12 = Sint_lin(rho1,u1,p1,beta1,alpha1,u23,-1);
                sintCD = SintCD(u23,p23,alpha3,alpha2);
                sint34 = Sint_lin(rho3,u23,p23,beta3,alpha3,u4,+1);
                sL(i-1) = sint12+sintCD;
                sR(i)   = sint34;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = F4+F3-Fsr;
                sint12  = Sint_lin(rho1,u1,p1,beta1,alpha1,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sintsr4 = Sint_lin(rhosr,usr,psr,betasr,alphasr,u4,+1);
                sint3sr = Sint_lin(rho3,u23,p23,beta3,alpha3,usr,+1);
                sL(i-1) = sint12+sintCD+sintsr4;
                sR(i)   = sint3sr;
            end

        elseif  lambda0 <= 0 && lambda3<=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = F1-Fsl+Fsr;
                sintsl2 = Sint_lin(rhosl,usl,psl,betasl,alphasl,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint3sr = Sint_lin(rho3,u23,p23,beta3,alpha3,usr,+1);
                sint1sl = Sint_lin(rho1,u1,p1,beta1,alpha1,usl,-1);
                sintsr4 = Sint_lin(rhosr,usr,psr,betasr,alphasr,u4,+1);
                sL(i-1) = sintsl2+sintCD+sint3sr;
                sR(i)   = sint1sl+sintsr4;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = F1-Fsl+F4;
                sintsl2 = Sint_lin(rhosl,usl,psl,betasl,alphasl,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint34  = Sint_lin(rho3,u23,p23,beta3,alpha3,u4,+1);
                sint1sl = Sint_lin(rho1,u1,p1,beta1,alpha1,usl,-1);
                sL(i-1) = sintsl2+sintCD+sint34;
                sR(i)   = sint1sl;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = Fsr;
                sint12  = Sint_lin(rho1,u1,p1,beta1,alpha1,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint3sr = Sint_lin(rho3,u23,p23,beta3,alpha3,usr,+1);
                sintsr4 = Sint_lin(rhosr,usr,psr,betasr,alphasr,u4,+1);
                sL(i-1) = sint12+sintCD+sint3sr;
                sR(i)   = sintsr4;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = F4;
                sint12 = Sint_lin(rho1,u1,p1,beta1,alpha1,u23,-1);
                sintCD = SintCD(u23,p23,alpha3,alpha2);
                sint34 = Sint_lin(rho3,u23,p23,beta3,alpha3,u4,+1);
                sL(i-1) = sint12+sintCD+sint34;
                sR(i)   = 0;
            end
        end %flux choice

    end %if L==R
    
end %for i loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left sonic point                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rhosl,usl,psl,betasl,alphasl]=linosher_sonic_L(rho1,u1,p1,beta1,alpha1,a1)

global gamma1 gamma2
global pi1 pi2

%-----------------------------------------------------------------------%
%                         Newton-Raphson solver                         %
%-----------------------------------------------------------------------%

uguess = a1;

diffu = 1;
m = 0;

while diffu>1e-6
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
%     b1      = alpha1*(1-alpha1)*(gamma1-gamma2)/((1-alpha1)*gamma1+alpha1*gamma2);
    b1      = alpha1*(1-alpha1)*(gamma1*(p1+pi1)-gamma2*(p1+pi2))/((1-alpha1)*gamma1*(p1+pi1)+alpha1*gamma2*(p1+pi2));
    alphasl = alpha1-b1/a1*(usl-u1);
%     asl     = sqrt(1/(alphasl/gamma1+(1-alphasl)/gamma2)*psl/rhosl);
    asl     = sqrt(1/(alphasl/(gamma1*(psl+pi1))+(1-alpha1)/(gamma2*(psl+pi2)))/rhosl);

%-----------------------------------------------------------------------%

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
betasl  = beta1;
% b1      = alpha1*(1-alpha1)*(gamma1-gamma2)/((1-alpha1)*gamma1+alpha1*gamma2);
b1      = alpha1*(1-alpha1)*(gamma1*(p1+pi1)-gamma2*(p1+pi2))/((1-alpha1)*gamma1*(p1+pi1)+alpha1*gamma2*(p1+pi2));
alphasl = alpha1-b1/a1*(usl-u1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Right sonic point                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rhosr,usr,psr,betasr,alphasr]=linosher_sonic_R(rho4,u4,p4,beta4,alpha4,a4)

global gamma1 gamma2
global pi1 pi2

%-----------------------------------------------------------------------%
%                         Newton-Raphson solver                         %
%-----------------------------------------------------------------------%

uguess = -a4;

diffu = 1;
m = 0;

while diffu>1e-6
    m=m+1;

%-----------------------------------------------------------------------%
%                               LEFT WAVE                               %
%-----------------------------------------------------------------------%
    if m==1
        usr = uguess;
    elseif m==2
        usr = 0.8*uguess;
    end

    psr     = p4+rho4*a4*(usr-u4);
    rhosr   = rho4+rho4/a4*(usr-u4);
%     b4      = alpha4*(1-alpha4)*(gamma1-gamma2)/((1-alpha4)*gamma1+alpha4*gamma2);
    b4      = alpha4*(1-alpha4)*(gamma1*(p4+pi1)-gamma2*(p4+pi2))/((1-alpha4)*gamma1*(p4+pi1)+alpha4*gamma2*(p4+pi2));
    alphasr = alpha4+b4/a4*(usr-u4);
%     asr     = sqrt(1/(alphasr/gamma1+(1-alphasr)/gamma2)*psr/rhosr);
    asr     = sqrt(1/(alphasr/(gamma1*(psr+pi1))+(1-alphasr)/(gamma2*(psr+pi2)))/rhosr);
    %-----------------------------------------------------------------------%

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
betasr  = beta4;
% b4      = alpha4*(1-alpha4)*(gamma1-gamma2)/((1-alpha4)*gamma1+alpha4*gamma2);
b4      = alpha4*(1-alpha4)*(gamma1*(p4+pi1)-gamma2*(p4+pi2))/((1-alpha4)*gamma1*(p4+pi1)+alpha4*gamma2*(p4+pi2));
alphasr = alpha4+b4/a4*(usr-u4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sourse term integral                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sint]=Sint_lin(rho,u,p,beta,alpha,u_end,lr)

global gamma1 gamma2
global pi1 pi2

% lr is a sign function for left running wave (-) or right running wave (+)

du   = (u_end - u);
% a    = sqrt(1/(alpha/gamma1+(1-alpha)/gamma2)*p/rho);
% bS   = alpha*(1-alpha)*(gamma1-gamma2)/((1-alpha)*gamma1+alpha*gamma2);
a    = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
bS   = alpha*(1-alpha)*(gamma1*(p+pi1)-gamma2*(p+pi2))/((1-alpha)*gamma1*(p+pi1)+alpha*gamma2*(p+pi2));
sint = (bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u)*du;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sourse term integral (Interface / Contact discontinuity                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sintCD]=SintCD(uf,pf,alpha3,alpha2)

sintCD = uf*pf*(alpha3-alpha2);

