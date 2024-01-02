%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              %
% Osher's flux solver                          %
%                                              %
% Written by: Jasper Kreeft  (2007)            %
% Updated by: Jasper Kreeft  (2015)            %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,sL,sR,Lambda] = OsherFlux(Wf,Lambda)

global N
global gamma1 gamma2
global pi1 pi2

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
    a1     = sqrt(1/(alpha1/(gamma1*(p1+pi1))+(1-alpha1)/(gamma2*(p1+pi2)))/rho1);

    rho4   = Wf(1,2*i-1);
    u4     = Wf(2,2*i-1);
    p4     = Wf(3,2*i-1);
    beta4  = Wf(4,2*i-1);
    alpha4 = Wf(5,2*i-1);
    a4     = sqrt(1/(alpha4/(gamma1*(p4+pi1))+(1-alpha4)/(gamma2*(p4+pi2)))/rho4);

    % Check for negative densities and pressures
    if rho1<=0
        disp('rho1 negative')
        rho1 = eps;
    end

    if p1<=0
        disp('p1 negative')
        p1 = eps;
    end

    if rho4<=0
        disp('rho4 negative')
        rho4 = eps;
    end

    if p4<=0
        disp('p4 negative')
        p4 = eps;
    end

    % State in beide cellen hetzelfde
    if abs(rho1-rho4)+abs(u1-u4)+abs(p1-p4)+abs(beta1-beta4)+abs(alpha1-alpha4)<1e-10

        F(:,i) = [rho1*u1;
                 rho1*u1^2+p1;
                 alpha1*gamma1/(gamma1-1)*(p1+pi1)*u1+(1-alpha1)*gamma2/(gamma2-1)*(p1+pi2)*u1+1/2*rho1*u1^3;   % (alpha1/(gamma1-1)+(1-alpha1)/(gamma2-1)+1)*u1*p1+
                 beta1*rho1*u1;
                 gamma1/(gamma1-1)*u1*(p1+pi1)*alpha1+1/2*rho1*u1^3*beta1 ];
        
        Lambda  = max(abs([u1 u1-a1 u4+a4 Lambda]));

    % stilstaande interface
    elseif abs(p1-p4)+abs(u1)+abs(u4)<1e-10

        F(2,i) = p1;

    % Single-fluid Osher
    elseif ((alpha1+alpha4<1e-10) || ((1-alpha1)+(1-alpha4)<1e-10)) ...
           && ((beta1+beta4<1e-10) || ((1-beta1)+(1-beta4)<1e-10)) && (pi1==0 && pi2==0)
        
        gamma = 1/2*(1/(alpha1/gamma1+(1-alpha1)/gamma2)+1/(alpha4/gamma1+(1-alpha4)/gamma2));

        % Fina1 state va1ues
        z      = (gamma-1)/(2*gamma);
        p23    = ((a1+a4-(gamma-1)/2*(u4-u1))/(a1/p1^z+a4/p4^z))^(1/z);
        u23    = ((p1/p4)^z*u1/a1+u4/a4+2/(gamma-1)*((p1/p4)^z-1))/(1/a1*(p1/p4)^z+1/a4);
        rho2   = rho1*(p23/p1)^(1/gamma);
        rho3   = rho4*(p23/p4)^(1/gamma);
        alpha2 = alpha1;
        beta2  = beta1;
        alpha3 = alpha4;
        beta3  = beta4;

        % Speed of sound
        a2 = sqrt(gamma*p23/rho2);
        a3 = sqrt(gamma*p23/rho3);
        
        % Sonic point va1ues
        usl     = (gamma-1)/(gamma+1)*u1+2*a1/(gamma+1);
        asl     = usl;
        rhosl   = rho1*(asl/a1)^(2/(gamma-1));
        psl     = p1*(rhosl/rho1)^gamma;
        alphasl = alpha1;
        betasl  = beta1;
        
        usr     = (gamma-1)/(gamma+1)*u4-2*a4/(gamma+1);
        asr     = -usr;
        rhosr   = rho4*(asr/a4)^(2/(gamma-1));
        psr     = p4*(rhosr/rho4)^gamma;
        alphasr = alpha4;
        betasr  = beta4;

        Fl = [rho1*u1;    rho1*u1^2+p1;    (alpha1/(gamma1-1)+(1-alpha1)/(gamma2-1)+1)*u1*p1+0.5*rho1*u1^3;       beta1*rho1*u1;    gamma1/(gamma1-1)*u1*p1*alpha1+1/2*rho1*u1^3*beta1      ];
        F2 = [rho2*u23;   rho2*u23^2+p23;  (alpha2/(gamma1-1)+(1-alpha2)/(gamma2-1)+1)*u23*p23+0.5*rho2*u23^3;    beta2*rho2*u23;   gamma1/(gamma1-1)*u23*p23*alpha2+1/2*rho2*u23^3*beta2   ];
        F3 = [rho3*u23;   rho3*u23^2+p23;  (alpha3/(gamma1-1)+(1-alpha3)/(gamma2-1)+1)*u23*p23+0.5*rho3*u23^3;    beta3*rho3*u23;   gamma1/(gamma1-1)*u23*p23*alpha3+1/2*rho3*u23^3*beta3   ];
        Fr = [rho4*u4;    rho4*u4^2+p4;    (alpha4/(gamma1-1)+(1-alpha4)/(gamma2-1)+1)*u4*p4+0.5*rho4*u4^3;       beta4*rho4*u4;    gamma1/(gamma1-1)*u4*p4*alpha4+1/2*rho4*u4^3*beta4      ];
        Fsl = [rhosl*usl; rhosl*usl^2+psl; (alphasl/(gamma1-1)+(1-alphasl)/(gamma2-1)+1)*usl*psl+0.5*rhosl*usl^3; betasl*rhosl*usl; gamma1/(gamma1-1)*usl*psl*alphasl+1/2*rhosl*usl^3*betasl];
        Fsr = [rhosr*usr; rhosr*usr^2+psr; (alphasr/(gamma1-1)+(1-alphasr)/(gamma2-1)+1)*usr*psr+0.5*rhosr*usr^3; betasr*rhosr*usr; gamma1/(gamma1-1)*usr*psr*alphasr+1/2*rhosr*usr^3*betasr];

        lambda0 = u23;
        lambda1 = u1-a1;
        lambda2 = u23-a2;
        lambda3 = u23+a3;
        lambda4 = u4+a4;
        Lambda  = max(abs([lambda0 lambda1 lambda2 lambda3 lambda4 Lambda]));

        if lambda0 >= 0 && lambda2>=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = Fl;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = Fl+Fr-Fsr;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = Fsl;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = Fsl-Fsr+Fr;
            end

        elseif  lambda0 >= 0 && lambda2<=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = Fl-Fsl+F2;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = Fl-Fsl+F2-Fsr+Fr;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = F2;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = Fr+F2-Fsr;
            end

        elseif  lambda0 <= 0 && lambda3>=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = Fl-Fsl+F3;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = Fl-Fsl+F3-Fsr+Fr;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = F3;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = Fr+F3-Fsr;
            end
    
        elseif  lambda0 <= 0 && lambda3<=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = Fl-Fsl+Fsr;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = Fl-Fsl+Fr;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = Fsr;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = Fr;
            end
        end

    % Two-fluid problem
    else

        %Final state values
        [rho2,rho3,u23,p23,beta2,beta3,alpha2,alpha3]=TFosher(rho1,u1,p1,beta1,alpha1,a1,rho4,u4,p4,beta4,alpha4,a4);

        %speed of sound
        a2 = sqrt(1/(alpha2/(gamma1*(p23+pi1))+(1-alpha2)/(gamma2*(p23+pi2)))/rho2);
        a3 = sqrt(1/(alpha3/(gamma1*(p23+pi1))+(1-alpha3)/(gamma2*(p23+pi2)))/rho3);

        %Sonic point values
        [rhosl,usl,psl,betasl,alphasl]=TFO_sonic_L(u4,p4,a4,rho1,u1,p1,beta1,alpha1,a1);
        [rhosr,usr,psr,betasr,alphasr]=TFO_sonic_R(u1,p1,a1,rho4,u4,p4,beta4,alpha4,a4);

        Fl  = [rho1*u1  ; rho1*u1^2+p1   ; alpha1*gamma1*(p1+pi1)*u1/(gamma1-1)+(1-alpha1)*gamma2*(p1+pi2)*u1/(gamma2-1)+0.5*rho1*u1^3        ; beta1*rho1*u1   ; alpha1*gamma1*(p1 +pi1)*u1/(gamma1-1)+1/2*rho1*u1^3*beta1      ];
        F2  = [rho2*u23 ; rho2*u23^2+p23 ; alpha2*gamma1*(p23+pi1)*u23/(gamma1-1)+(1-alpha2)*gamma2*(p23+pi2)*u23/(gamma2-1)+0.5*rho2*u23^3   ; beta2*rho2*u23  ; alpha2*gamma1*(p23+pi1)*u23/(gamma1-1)+1/2*rho2*u23^3*beta2    ];
        F3  = [rho3*u23 ; rho3*u23^2+p23 ; alpha3*gamma1*(p23+pi1)*u23/(gamma1-1)+(1-alpha3)*gamma2*(p23+pi2)*u23/(gamma2-1)+0.5*rho3*u23^3   ; beta3*rho3*u23  ; alpha3*gamma1*(p23+pi1)*u23/(gamma1-1)+1/2*rho3*u23^3*beta3    ];
        Fr  = [rho4*u4  ; rho4*u4^2+p4   ; alpha4*gamma1*(p4+pi1)*u4/(gamma1-1)+(1-alpha4)*gamma2*(p4+pi2)*u4/(gamma2-1)+0.5*rho4*u4^3        ; beta4*rho4*u4   ; alpha4*gamma1*(p4 +pi1)*u4/(gamma1-1)+1/2*rho4*u4^3*beta4      ];
        Fsl = [rhosl*usl; rhosl*usl^2+psl; alphasl*gamma1*(psl+pi1)*usl/(gamma1-1)+(1-alphasl)*gamma2*(psl+pi2)*usl/(gamma2-1)+0.5*rhosl*usl^3; betasl*rhosl*usl; alphasl*gamma1*(psl+pi1)*usl/(gamma1-1)+1/2*rhosl*usl^3*betasl ];
        Fsr = [rhosr*usr; rhosr*usr^2+psr; alphasr*gamma1*(psr+pi1)*usr/(gamma1-1)+(1-alphasr)*gamma2*(psr+pi2)*usr/(gamma2-1)+0.5*rhosr*usr^3; betasr*rhosr*usr; alphasr*gamma1*(psr+pi1)*usr/(gamma1-1)+1/2*rhosr*usr^3*betasr ];

        lambda0 = u23;
        lambda1 = u1-a1;
        lambda2 = u23-a2;
        lambda3 = u23+a3;
        lambda4 = u4+a4;
        Lambda  = max(abs([lambda0 lambda1 lambda2 lambda3 lambda4 Lambda]));

        if lambda0>=0 && lambda2>=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = Fl;
                sint12 = Sint(rho1,u1,p1,beta1,alpha1,u23,-1);
                sintCD = SintCD(u23,p23,alpha3,alpha2);
                sint34 = Sint(rho4,u4,p4,beta4,alpha4,u23,+1);
                sL(i-1) = 0;
                sR(i)   = sint12+sintCD-sint34;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = Fl+Fr-Fsr;
                sintsr4 = Sint(rho4,u4,p4,beta4,alpha4,usr,+1);
                sint12  = Sint(rho1,u1,p1,beta1,alpha1,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint3sr = Sint(rhosr,usr,psr,betasr,alphasr,u23,+1);
                sL(i-1) = -sintsr4;
                sR(i)   = sint12+sintCD-sint3sr;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = Fsl;
                sint1sl = Sint(rho1,u1,p1,beta1,alpha1,usl,-1);
                sintsl2 = Sint(rhosl,usl,psl,betasl,alphasl,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint34 = Sint(rho4,u4,p4,beta4,alpha4,u23,+1);
                sL(i-1) = sint1sl;
                sR(i)   = sintsl2+sintCD-sint34;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = Fsl-Fsr+Fr;
                sint1sl = Sint(rho1,u1,p1,beta1,alpha1,usl,-1);
                sintsr4 = Sint(rho4,u4,p4,beta4,alpha4,usr,+1);
                sintsl2 = Sint(rhosl,usl,psl,betasl,alphasl,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint3sr = Sint(rhosr,usr,psr,betasr,alphasr,u23,+1);
                sL(i-1) = sint1sl-sintsr4;
                sR(i)   = sintsl2+sintCD-sint3sr;
            end

        elseif lambda0>=0 && lambda2<=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = Fl-Fsl+F2;
                sintsl2 =Sint(rhosl,usl,psl,betasl,alphasl,u23,-1);
                sint1sl =Sint(rho1,u1,p1,beta1,alpha1,usl,-1);
                sintCD  =SintCD(u23,p23,alpha3,alpha2);
                sint34 = Sint(rho4,u4,p4,beta4,alpha4,u23,+1);
                sL(i-1) = sintsl2;
                sR(i)   = sint1sl+sintCD-sint34;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = Fl-Fsl+F2-Fsr+Fr;
                sintsl2 = Sint(rhosl,usl,psl,betasl,alphasl,u23,-1);
                sintsr4 = Sint(rho4,u4,p4,beta4,alpha4,usr,+1);
                sint1sl = Sint(rho1,u1,p1,beta1,alpha1,usl,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint3sr = Sint(rhosr,usr,psr,betasr,alphasr,u23,+1);
                sL(i-1) = sintsl2-sintsr4;
                sR(i)   = sint1sl+sintCD-sint3sr;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = F2;
                sint12 = Sint(rho1,u1,p1,beta1,alpha1,u23,-1);
                sintCD = SintCD(u23,p23,alpha3,alpha2);
                sint34 = Sint(rho4,u4,p4,beta4,alpha4,u23,+1);
                sL(i-1) = sint12;
                sR(i)   = sintCD-sint34;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = Fr+F2-Fsr;
                sint12  = Sint(rho1,u1,p1,beta1,alpha1,u23,-1);
                sintsr4 = Sint(rho4,u4,p4,beta4,alpha4,usr,+1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint3sr = Sint(rhosr,usr,psr,betasr,alphasr,u23,+1);
                sL(i-1) = sint12-sintsr4;
                sR(i)   = sintCD-sint3sr;
            end

        elseif lambda0<=0 && lambda3>=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = Fl-Fsl+F3;
                sintsl2 = Sint(rhosl,usl,psl,betasl,alphasl,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint1sl = Sint(rho1,u1,p1,beta1,alpha1,usl,-1);
                sint34 = Sint(rho4,u4,p4,beta4,alpha4,u23,+1);
                sL(i-1) = sintsl2+sintCD;
                sR(i)   = sint1sl-sint34;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = Fl-Fsl+F3-Fsr+Fr;
                sintsl2 = Sint(rhosl,usl,psl,betasl,alphasl,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sintsr4 = Sint(rho4,u4,p4,beta4,alpha4,usr,+1);
                sint1sl = Sint(rho1,u1,p1,beta1,alpha1,usl,-1);
                sint3sr = Sint(rhosr,usr,psr,betasr,alphasr,u23,+1);
                sL(i-1) = sintsl2+sintCD-sintsr4;
                sR(i)   = sint1sl-sint3sr;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = F3;
                sint12 = Sint(rho1,u1,p1,beta1,alpha1,u23,-1);
                sintCD = SintCD(u23,p23,alpha3,alpha2);
                sint34 = Sint(rho4,u4,p4,beta4,alpha4,u23,+1);
                sL(i-1) = sint12+sintCD;
                sR(i)   = -sint34;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = Fr+F3-Fsr;
                sint12  = Sint(rho1,u1,p1,beta1,alpha1,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sintsr4 = Sint(rho4,u4,p4,beta4,alpha4,usr,+1);
                sint3sr = Sint(rhosr,usr,psr,betasr,alphasr,u23,+1);
                sL(i-1) = sint12+sintCD-sintsr4;
                sR(i)   = -sint3sr;
            end

        elseif lambda0<=0 && lambda3<=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = Fl-Fsl+Fsr;
                sintsl2 = Sint(rhosl,usl,psl,betasl,alphasl,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint3sr = Sint(rhosr,usr,psr,betasr,alphasr,u23,+1);
                sint1sl = Sint(rho1,u1,p1,beta1,alpha1,usl,-1);
                sintsr4 = Sint(rho4,u4,p4,beta4,alpha4,usr,+1);
                sL(i-1) = sintsl2+sintCD-sint3sr;
                sR(i)   = sint1sl-sintsr4;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = Fl-Fsl+Fr;
                sintsl2 = Sint(rhosl,usl,psl,betasl,alphasl,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint34 = Sint(rho4,u4,p4,beta4,alpha4,u23,+1);
                sint1sl = Sint(rho1,u1,p1,beta1,alpha1,usl,-1);
                sL(i-1) = sintsl2+sintCD-sint34;
                sR(i)   = sint1sl;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = Fsr;
                sint12  = Sint(rho1,u1,p1,beta1,alpha1,u23,-1);
                sintCD  = SintCD(u23,p23,alpha3,alpha2);
                sint3sr = Sint(rhosr,usr,psr,betasr,alphasr,u23,+1);
                sintsr4 = Sint(rho4,u4,p4,beta4,alpha4,usr,+1);
                sL(i-1) = sint12+sintCD-sint3sr;
                sR(i)   = -sintsr4;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = Fr;
                sint12 = Sint(rho1,u1,p1,beta1,alpha1,u23,-1);
                sintCD = SintCD(u23,p23,alpha3,alpha2);
                sint34 = Sint(rho4,u4,p4,beta4,alpha4,u23,+1);
                sL(i-1) = sint12+sintCD-sint34;
                sR(i)   = 0;
            end
        end %flux choice

    end %if L==R
end %for i loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two-Fluid Osher                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho2,rho3,u23,p23,beta2,beta3,alpha2,alpha3]=TFosher(rhol,ul,pl,betal,alphal,al,rhor,ur,pr,betar,alphar,ar)

global gamma1 gamma2
global pi1 pi2

NRK4   = 6;

%----------------------------%
% Single-fluid initial guess %
%----------------------------%
gamma  = (1/(alphal/gamma1+(1-alphal)/gamma2)+1/(alphar/gamma1+(1-alphar)/gamma2))/2;
% gamma  = (1/(alphal/(gamma1-1)+(1-alphal)/(gamma2-1))+1/(alphar/(gamma1-1)+(1-alphar)/(gamma2-1)))/2+1;
pil = (alphal*gamma1*pi1/(gamma1-1)+(1-alphal)*gamma2*pi2/(gamma2-1))/(alphal/(gamma1-1)+(1-alphal)/(gamma-1));
pir = (alphar*gamma1*pi1/(gamma1-1)+(1-alphar)*gamma2*pi2/(gamma2-1))/(alphar/(gamma1-1)+(1-alphar)/(gamma-1));
uguess = (ul/al*((pl+pil)/(pr+pir))^((gamma-1)/(2*gamma))+ur/ar+2/(gamma-1)*(((pl+pil)/(pr+pir))^((gamma-1)/(2*gamma))-1))/...
         (((pl+pil)/(pr+pir))^((gamma-1)/(2*gamma))/al+1/ar);
% uguess = (ul/al*pl^((gamma-1)/(2*gamma))+ur/ar*pr^((gamma-1)/(2*gamma))+2/(gamma-1)*(pl^((gamma-1)/(2*gamma))-pr^((gamma-1)/(2*gamma))))/(pl^((gamma-1)/(2*gamma))/al+pr^((gamma-1)/(2*gamma))/ar);

%--------------------------%
% Runga-Kutta fourth order %
%--------------------------%
diffp = 1;
m = 0;

while diffp>1e-6
    m=m+1;

    %--------------------------%
    % Left wave                %
    %--------------------------%
    if m == 1
        uf = 1.0*uguess;
    elseif m==2
        uf = 0.8*uguess;
    end

    rho   = rhol;
    p     = pl;
    beta  = betal;
    alpha = alphal;
    du    = (uf - ul)/NRK4;

    for k=1:NRK4
        [dp_k1,drho_k1,dalpha_k1] = RK4_u(rho          ,p        ,alpha            ,du,-1);
        [dp_k2,drho_k2,dalpha_k2] = RK4_u(rho+drho_k1/2,p+dp_k1/2,alpha+dalpha_k1/2,du,-1);
        [dp_k3,drho_k3,dalpha_k3] = RK4_u(rho+drho_k2/2,p+dp_k2/2,alpha+dalpha_k2/2,du,-1);
        [dp_k4,drho_k4,dalpha_k4] = RK4_u(rho+drho_k3  ,p+dp_k3  ,alpha+dalpha_k3  ,du,-1);
        p     = p     + (dp_k1     + 2*dp_k2     + 2*dp_k3     + dp_k4    )/6;
        rho   = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
        alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
    end

    rho2   = rho;
    pf2    = p;
    beta2  = beta;
    alpha2 = alpha;

    %--------------------------%
    % Right wave               %
    %--------------------------%
    rho   = rhor;
    p     = pr;
    beta  = betar;
    alpha = alphar;
    du    = (uf - ur)/NRK4;

    for k=1:NRK4
        [dp_k1,drho_k1,dalpha_k1] = RK4_u(rho          ,p        ,alpha            ,du,+1);
        [dp_k2,drho_k2,dalpha_k2] = RK4_u(rho+drho_k1/2,p+dp_k1/2,alpha+dalpha_k1/2,du,+1);
        [dp_k3,drho_k3,dalpha_k3] = RK4_u(rho+drho_k2/2,p+dp_k2/2,alpha+dalpha_k2/2,du,+1);
        [dp_k4,drho_k4,dalpha_k4] = RK4_u(rho+drho_k3  ,p+dp_k3  ,alpha+dalpha_k3  ,du,+1);
        p     = p     + (dp_k1     + 2*dp_k2     + 2*dp_k3     + dp_k4    )/6;
        rho   = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
        alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
    end

    rho3   = rho;
    pf3    = p;
    beta3  = beta;
    alpha3 = alpha;

    dp    = pf3 - pf2;
    diffp = abs(dp);

    if m==1
        ufold = uf;
        dpold = dp;
    else
        ufnew = uf - dp*(uf - ufold + eps)/(dp - dpold + eps);
        ufold = uf;
        uf    = ufnew;
        dpold = dp;
    end

end %while loop

u23 = uf;
p23 = pf2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left sonic point                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rhosl,usl,psl,betasl,alphasl]=TFO_sonic_L(u4,p4,a4,rho1,u1,p1,beta1,alpha1,a1)

global gamma1 gamma2
global pi1 pi2

%Runga-Kutta Solver
NRK4   = 6;
gamma  = alpha1*gamma1+(1-alpha1)*gamma2;
pim = (alpha1*gamma1*pi1/(gamma1-1)+(1-alpha1)*gamma2*pi2/(gamma2-1))/(alpha1/(gamma1-1)+(1-alpha1)/(gamma-1));
uguess = (u1/a1*((p1+pim)/(p4+pim))^((gamma-1)/(2*gamma))+u4/a4+2/(gamma-1)*(((p1+pim)/(p4+pim))^((gamma-1)/(2*gamma))-1))/...
         (((p1+pim)/(p4+pim))^((gamma-1)/(2*gamma))/a1+1/a4);
% uguess = (u1/a1*p1^((gamma-1)/(2*gamma))+u4/a4*p4^((gamma-1)/(2*gamma))+...
%          2/(gamma-1)*(p1^((gamma-1)/(2*gamma))-p4^((gamma-1)/(2*gamma))))/(p1^((gamma-1)/(2*gamma))/a1+p4^((gamma-1)/(2*gamma))/a4);

%--------------------------%
% Runga-Kutta fourth order %
%--------------------------%
diffu = 1;
m = 0;

while diffu>1e-6
    m=m+1;

    %--------------------------%
    % Left wave                %
    %--------------------------%
    if m == 1
        usl = 1.0*uguess;
    elseif m == 2
        usl = 0.8*uguess;
    end

    p     = p1;
    rho   = rho1;
    beta  = beta1;
    alpha = alpha1;
    du    = (usl - u1)/NRK4;

    for k=1:NRK4
        [dp_k1,drho_k1,dalpha_k1] = RK4_u(rho          ,p        ,alpha            ,du,-1);
        [dp_k2,drho_k2,dalpha_k2] = RK4_u(rho+drho_k1/2,p+dp_k1/2,alpha+dalpha_k1/2,du,-1);
        [dp_k3,drho_k3,dalpha_k3] = RK4_u(rho+drho_k2/2,p+dp_k2/2,alpha+dalpha_k2/2,du,-1);
        [dp_k4,drho_k4,dalpha_k4] = RK4_u(rho+drho_k3  ,p+dp_k3  ,alpha+dalpha_k3  ,du,-1);
        p     = p     + (dp_k1     + 2*dp_k2     + 2*dp_k3     + dp_k4    )/6;
        rho   = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
        alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
    end

    rhosl   = rho;
    psl     = p;
    alphasl = alpha;
    betasl  = beta;
    asl     = sqrt(1/(alphasl/(gamma1*(psl+pi1))+(1-alphasl)/(gamma2*(psl+pi2)))/rhosl);


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Right sonic point                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rhosr,usr,psr,betasr,alphasr]=TFO_sonic_R(u1,p1,a1,rho4,u4,p4,beta4,alpha4,a4)

global gamma1 gamma2
global pi1 pi2

%Runga-Kutta Solver
NRK4   = 6;
gamma  = alpha4*gamma1+(1-alpha4)*gamma2;
pim = (alpha4*gamma1*pi1/(gamma1-1)+(1-alpha4)*gamma2*pi2/(gamma2-1))/(alpha4/(gamma1-1)+(1-alpha4)/(gamma-1));
uguess = (u1/a1*((p1+pim)/(p4+pim))^((gamma-1)/(2*gamma))+u4/a4+2/(gamma-1)*(((p1+pim)/(p4+pim))^((gamma-1)/(2*gamma))-1))/...
         (((p1+pim)/(p4+pim))^((gamma-1)/(2*gamma))/a1+1/a4);
% uguess = (u1/a1*p1^((gamma-1)/(2*gamma))+u4/a4*p4^((gamma-1)/(2*gamma))+2/(gamma-1)*(p1^((gamma-1)/(2*gamma))...
%          -p4^((gamma-1)/(2*gamma))))/(p1^((gamma-1)/(2*gamma))/a1+p4^((gamma-1)/(2*gamma))/a4);
     
%--------------------------%
% Runga-Kutta fourth order %
%--------------------------%
diffu = 1;
m = 0;

while diffu>1e-6
    m=m+1;

    %--------------------------%
    % Right wave               %
    %--------------------------%
    if m == 1
        usr = 1.0*uguess;
    elseif m == 2
        usr = 0.8*uguess;
    end

    rho   = rho4;
    p     = p4;
    beta  = beta4;
    alpha = alpha4;
    du    = (usr - u4)/NRK4;

    for k=1:NRK4
        [dp_k1,drho_k1,dalpha_k1] = RK4_u(rho          ,p        ,alpha            ,du,+1);
        [dp_k2,drho_k2,dalpha_k2] = RK4_u(rho+drho_k1/2,p+dp_k1/2,alpha+dalpha_k1/2,du,+1);
        [dp_k3,drho_k3,dalpha_k3] = RK4_u(rho+drho_k2/2,p+dp_k2/2,alpha+dalpha_k2/2,du,+1);
        [dp_k4,drho_k4,dalpha_k4] = RK4_u(rho+drho_k3  ,p+dp_k3  ,alpha+dalpha_k3  ,du,+1);
        p     = p     + (dp_k1     + 2*dp_k2     + 2*dp_k3     + dp_k4    )/6;
        rho   = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
        alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
    end

    psr     = p;
    alphasr = alpha;
    betasr  = beta;
    rhosr   = rho;
    asr     = sqrt(1/(alphasr/(gamma1*(psr+pi1))+(1-alphasr)/(gamma2*(psr+pi2)))/rhosr);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sourse term integral                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sint=Sint(rho,u,p,beta,alpha,u_end,lr)

global gamma1 gamma2
global pi1 pi2

% lr is a sign function for left running wave (-) or right running wave (+)

s = [0 0 0 0];

du   = (u_end - u);

a     = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
bS    = (gamma1*(p+pi1)-gamma2*(p+pi2))/(gamma1*(p+pi1)/alpha+gamma2*(p+pi2)/(1-alpha));
s(1) = bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u;

for k=1:3
[dp_k1,drho_k1,dalpha_k1] = RK4_u(rho          ,p        ,alpha            ,du/3,lr);
[dp_k2,drho_k2,dalpha_k2] = RK4_u(rho+drho_k1/2,p+dp_k1/2,alpha+dalpha_k1/2,du/3,lr);
[dp_k3,drho_k3,dalpha_k3] = RK4_u(rho+drho_k2/2,p+dp_k2/2,alpha+dalpha_k2/2,du/3,lr);
[dp_k4,drho_k4,dalpha_k4] = RK4_u(rho+drho_k3  ,p+dp_k3  ,alpha+dalpha_k3  ,du/3,lr);
p      = p     + (dp_k1     + 2*dp_k2     + 2*dp_k3     + dp_k4    )/6;
rho    = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
alpha  = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
u      = u+du/3;
a      = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
bS     = (gamma1*(p+pi1)-gamma2*(p+pi2))/(gamma1*(p+pi1)/alpha+gamma2*(p+pi2)/(1-alpha));
s(k+1) = bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u;
end

sint = du/8*(s(1)+3*s(2)+3*s(3)+s(4));
% Simpson's Rule (third order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sourse term integral (Interface / Contact discontinuity                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sintCD=SintCD(uf,pf,alpha3,alpha2)

sintCD = uf*pf*(alpha3-alpha2);