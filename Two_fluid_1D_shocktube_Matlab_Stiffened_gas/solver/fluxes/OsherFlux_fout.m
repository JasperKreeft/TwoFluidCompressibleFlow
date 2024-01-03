%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Osher's flux solver                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,sL,sR,Lambda] = OsherFlux(N,Wf,gamma1,gamma2,pi1,pi2,Lambda)

% eps = 1e-12;

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
        keyboard
        rho1 = eps;
    end

    if p1<=0
        disp('p1 negative')
        keyboard
        p1 = eps;
    end

    if rho4<=0
        disp('rho4 negative')
        keyboard
        rho4 = eps;
    end

    if p4<=0
        disp('p4 negative')
        keyboard
        p4 = eps;
    end

    % % cavity test
    % if (u1+2/(gamma-1)*a1)<(u4-2/(gamma-1)*a4)
    %     disp('vacuum is created')
    %     pause
    % end

    % State in beide cellen hetzelfde
    if abs(rho1-rho4)<1e-10 && abs(u1-u4)<1e-10 && abs(p1-p4)<1e-10 && abs(beta1-beta4)<1e-10 && abs(alpha1-alpha4)<1e-10
    
       F(:,i) = [rho1*u1;
                  rho1*u1^2+p1;
                  alpha1*gamma1*(p1+pi1)*u1/(gamma1-1)+(1-alpha1)*gamma2*(p1+pi2)*u1/(gamma2-1)+1/2*rho1*u1^3;
                  rho1*u1*beta1;
                  alpha1*gamma1*(p1+pi1)*u1/(gamma1-1)+1/2*rho1*u1^3*beta1 ];

        sL(i-1) = 0;
        sR(i)   = 0;

        Lambda  = max(abs([u1 u1-a1 u4+a4 Lambda]));


    % supersonic tests niet mogelijk. Golf patroon moet bekend zijn voor
    % het berekenen van de source term!!!!

    % stilstaande interface
    elseif abs(p1-p4)<1e-10 && abs(u1)<1e-10 && abs(u4)<1e-10

        F(:,i) = [0;p1;0;0;0];

        sL(i-1) = 0;
        sR(i)   = 0;

    % Two-fluid problem
    else

        %Final state values
        [rho2,rho3,u23,p23,beta2,beta3,alpha2,alpha3]=TFosher(rho1,u1,p1,beta1,alpha1,a1,rho4,u4,p4,beta4,alpha4,a4,gamma1,gamma2,pi1,pi2); %OK
% keyboard
        %speed of sound
        a2 = sqrt(1/(alpha2/(gamma1*(p23+pi1))+(1-alpha2)/(gamma2*(p23+pi2)))/rho2);
        a3 = sqrt(1/(alpha3/(gamma1*(p23+pi1))+(1-alpha3)/(gamma2*(p23+pi2)))/rho3);

        %Sonic point values
        [rhosl,usl,psl,betasl,alphasl]=TFO_sonic_L(u4,p4,a4,rho1,u1,p1,beta1,alpha1,a1,gamma1,gamma2,pi1,pi2);
        [rhosr,usr,psr,betasr,alphasr]=TFO_sonic_R(u1,p1,a1,rho4,u4,p4,beta4,alpha4,a4,gamma1,gamma2,pi1,pi2);

        Fl =  [rho1*u1  ; rho1*u1^2+p1   ; alpha1*gamma1*(p1+pi1)*u1/(gamma1-1)+(1-alpha1)*gamma2*(p1+pi2)*u1/(gamma2-1)+0.5*rho1*u1^3        ; beta1*rho1*u1   ; alpha1*gamma1*(p1 +pi1)*u1/(gamma1-1)+1/2*rho1*u1^3*beta1      ];
        F2 =  [rho2*u23 ; rho2*u23^2+p23 ; alpha2*gamma1*(p23+pi1)*u23/(gamma1-1)+(1-alpha2)*gamma2*(p23+pi2)*u23/(gamma2-1)+0.5*rho2*u23^3   ; beta2*rho2*u23  ; alpha2*gamma1*(p23+pi1)*u23/(gamma1-1)+1/2*rho2*u23^3*beta2    ];
        F3 =  [rho3*u23 ; rho3*u23^2+p23 ; alpha3*gamma1*(p23+pi1)*u23/(gamma1-1)+(1-alpha3)*gamma2*(p23+pi2)*u23/(gamma2-1)+0.5*rho3*u23^3   ; beta3*rho3*u23  ; alpha3*gamma1*(p23+pi1)*u23/(gamma1-1)+1/2*rho3*u23^3*beta3    ];
        Fr =  [rho4*u4  ; rho4*u4^2+p4   ; alpha4*gamma1*(p4+pi1)*u4/(gamma1-1)+(1-alpha4)*gamma2*(p4+pi2)*u4/(gamma2-1)+0.5*rho4*u4^3        ; beta4*rho4*u4   ; alpha4*gamma1*(p4 +pi1)*u4/(gamma1-1)+1/2*rho4*u4^3*beta4      ];
        Fsl = [rhosl*usl; rhosl*usl^2+psl; alphasl*gamma1*(psl+pi1)*usl/(gamma1-1)+(1-alphasl)*gamma2*(psl+pi2)*usl/(gamma2-1)+0.5*rhosl*usl^3; betasl*rhosl*usl; alphasl*gamma1*(psl+pi1)*usl/(gamma1-1)+1/2*rhosl*usl^3*betasl ];
        Fsr = [rhosr*usr; rhosr*usr^2+psr; alphasr*gamma1*(psr+pi1)*usr/(gamma1-1)+(1-alphasr)*gamma2*(psr+pi2)*usr/(gamma2-1)+0.5*rhosr*usr^3; betasr*rhosr*usr; alphasr*gamma1*(psr+pi1)*usr/(gamma1-1)+1/2*rhosr*usr^3*betasr ];

        lambda0 = u23;
        lambda1 = u1-a1;
        lambda2 = u23-a2;
        lambda3 = u23+a3;
        lambda4 = u4+a4;
        Lambda  = max(abs([lambda0 lambda1 lambda2 lambda3 lambda4 Lambda]));

        if lambda0 >= 0 && lambda2>=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = Fl;
                [sint12] = Sint(rho1,u1,p1,beta1,alpha1,u23,gamma1,gamma2,pi1,pi2,-1);
                [sintCD] = SintCD(u23,p23,alpha3,alpha2);
                [sint34] = Sint(rho3,u23,p23,beta3,alpha3,u4,gamma1,gamma2,pi1,pi2,+1);
                sL(i-1) = 0;
                sR(i)   = sint12+sintCD+sint34;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = Fl+Fr-Fsr;
                [sintsr4] = Sint(rhosr,usr,psr,betasr,alphasr,u4,gamma1,gamma2,pi1,pi2,+1);
                [sint12]  = Sint(rho1,u1,p1,beta1,alpha1,u23,gamma1,gamma2,pi1,pi2,-1);
                [sintCD]  = SintCD(u23,p23,alpha3,alpha2);
                [sint3sr] = Sint(rho3,u23,p23,beta3,alpha3,usr,gamma1,gamma2,pi1,pi2,+1);
                sL(i-1) = sintsr4;
                sR(i)   = sint12+sintCD+sint3sr;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = Fsl;
                [sint1sl] = Sint(rho1,u1,p1,beta1,alpha1,usl,gamma1,gamma2,pi1,pi2,-1);
                [sintsl2] = Sint(rhosl,usl,psl,betasl,alphasl,u23,gamma1,gamma2,pi1,pi2,-1);
                [sintCD]  = SintCD(u23,p23,alpha3,alpha2);
                [sint34]  = Sint(rho3,u23,p23,beta3,alpha3,u4,gamma1,gamma2,pi1,pi2,+1);
                sL(i-1) = sint1sl;
                sR(i)   = sintsl2+sintCD+sint34;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = Fsl-Fsr+Fr;
                [sint1sl] = Sint(rho1,u1,p1,beta1,alpha1,usl,gamma1,gamma2,pi1,pi2,-1);
                [sintsr4] = Sint(rhosr,usr,psr,betasr,alphasr,u4,gamma1,gamma2,pi1,pi2,+1);
                [sintsl2] = Sint(rhosl,usl,psl,betasl,alphasl,u23,gamma1,gamma2,pi1,pi2,-1);
                [sintCD]  = SintCD(u23,p23,alpha3,alpha2);
                [sint3sr] = Sint(rho3,u23,p23,beta3,alpha3,usr,gamma1,gamma2,pi1,pi2,+1);
                sL(i-1) = sint1sl+sintsr4;
                sR(i)   = sintsl2+sintCD+sint3sr;
            end

        elseif  lambda0 >= 0 && lambda2<=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = Fl-Fsl+F2;
                [sintsl2] =Sint(rhosl,usl,psl,betasl,alphasl,u23,gamma1,gamma2,pi1,pi2,-1);
                [sint1sl] =Sint(rho1,u1,p1,beta1,alpha1,usl,gamma1,gamma2,pi1,pi2,-1);
                [sintCD]  =SintCD(u23,p23,alpha3,alpha2);
                [sint34]  =Sint(rho3,u23,p23,beta3,alpha3,u4,gamma1,gamma2,pi1,pi2,+1);
                sL(i-1) = sintsl2;
                sR(i)   = sint1sl+sintCD+sint34;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = Fl-Fsl+F2-Fsr+Fr;
                [sintsl2] = Sint(rhosl,usl,psl,betasl,alphasl,u23,gamma1,gamma2,pi1,pi2,-1);
                [sintsr4] = Sint(rhosr,usr,psr,betasr,alphasr,u4,gamma1,gamma2,pi1,pi2,+1);
                [sint1sl] = Sint(rho1,u1,p1,beta1,alpha1,usl,gamma1,gamma2,pi1,pi2,-1);
                [sintCD]  = SintCD(u23,p23,alpha3,alpha2);
                [sint3sr] = Sint(rho3,u23,p23,beta3,alpha3,usr,gamma1,gamma2,pi1,pi2,+1);
                sL(i-1) = sintsl2+sintsr4;
                sR(i)   = sint1sl+sintCD+sint3sr;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = F2;
                [sint12] = Sint(rho1,u1,p1,beta1,alpha1,u23,gamma1,gamma2,pi1,pi2,-1);
                [sintCD] = SintCD(u23,p23,alpha3,alpha2);
                [sint34] = Sint(rho3,u23,p23,beta3,alpha3,u4,gamma1,gamma2,pi1,pi2,+1);
                sL(i-1) = sint12;
                sR(i)   = sintCD+sint34;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = Fr+F2-Fsr;
                [sint12]  = Sint(rho1,u1,p1,beta1,alpha1,u23,gamma1,gamma2,pi1,pi2,-1);
                [sintsr4] = Sint(rhosr,usr,psr,betasr,alphasr,u4,gamma1,gamma2,pi1,pi2,+1);
                [sintCD]  = SintCD(u23,p23,alpha3,alpha2);
                [sint3sr] = Sint(rho3,u23,p23,beta3,alpha3,usr,gamma1,gamma2,pi1,pi2,+1);
                sL(i-1) = sint12+sintsr4;
                sR(i)   = sintCD+sint3sr;
            end

        elseif  lambda0 <= 0 && lambda3>=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = Fl-Fsl+F3;
                [sintsl2] = Sint(rhosl,usl,psl,betasl,alphasl,u23,gamma1,gamma2,pi1,pi2,-1);
                [sintCD]  = SintCD(u23,p23,alpha3,alpha2);
                [sint1sl] = Sint(rho1,u1,p1,beta1,alpha1,usl,gamma1,gamma2,pi1,pi2,-1);
                [sint34]  = Sint(rho3,u23,p23,beta3,alpha3,u4,gamma1,gamma2,pi1,pi2,+1);
                sL(i-1) = sintsl2+sintCD;
                sR(i)   = sint1sl+sint34;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = Fl-Fsl+F3-Fsr+Fr;
                [sintsl2] = Sint(rhosl,usl,psl,betasl,alphasl,u23,gamma1,gamma2,pi1,pi2,-1);
                [sintCD]  = SintCD(u23,p23,alpha3,alpha2);
                [sintsr4] = Sint(rhosr,usr,psr,betasr,alphasr,u4,gamma1,gamma2,pi1,pi2,+1);
                [sint1sl] = Sint(rho1,u1,p1,beta1,alpha1,usl,gamma1,gamma2,pi1,pi2,-1);
                [sint3sr] = Sint(rho3,u23,p23,beta3,alpha3,usr,gamma1,gamma2,pi1,pi2,+1);
                sL(i-1) = sintsl2+sintCD+sintsr4;
                sR(i)   = sint1sl+sint3sr;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = F3;
                [sint12] = Sint(rho1,u1,p1,beta1,alpha1,u23,gamma1,gamma2,pi1,pi2,-1);
                [sintCD] = SintCD(u23,p23,alpha3,alpha2);
                [sint34] = Sint(rho3,u23,p23,beta3,alpha3,u4,gamma1,gamma2,pi1,pi2,+1);
                sL(i-1) = sint12+sintCD;
                sR(i)   = sint34;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = Fr+F3-Fsr;
                [sint12]  = Sint(rho1,u1,p1,beta1,alpha1,u23,gamma1,gamma2,pi1,pi2,-1);
                [sintCD]  = SintCD(u23,p23,alpha3,alpha2);
                [sintsr4] = Sint(rhosr,usr,psr,betasr,alphasr,u4,gamma1,gamma2,pi1,pi2,+1);
                [sint3sr] = Sint(rho3,u23,p23,beta3,alpha3,usr,gamma1,gamma2,pi1,pi2,+1);
                sL(i-1) = sint12+sintCD+sintsr4;
                sR(i)   = sint3sr;
            end

        elseif  lambda0 <= 0 && lambda3<=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = Fl-Fsl+Fsr;
                [sintsl2] = Sint(rhosl,usl,psl,betasl,alphasl,u23,gamma1,gamma2,pi1,pi2,-1);
                [sintCD]  = SintCD(u23,p23,alpha3,alpha2);
                [sint3sr] = Sint(rho3,u23,p23,beta3,alpha3,usr,gamma1,gamma2,pi1,pi2,+1);
                [sint1sl] = Sint(rho1,u1,p1,beta1,alpha1,usl,gamma1,gamma2,pi1,pi2,-1);
                [sintsr4] = Sint(rhosr,usr,psr,betasr,alphasr,u4,gamma1,gamma2,pi1,pi2,+1);
                sL(i-1) = sintsl2+sintCD+sint3sr;
                sR(i)   = sint1sl+sintsr4;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = Fl-Fsl+Fr;
                [sintsl2] = Sint(rhosl,usl,psl,betasl,alphasl,u23,gamma1,gamma2,pi1,pi2,-1);
                [sintCD]  = SintCD(u23,p23,alpha3,alpha2);
                [sint34]  = Sint(rho3,u23,p23,beta3,alpha3,u4,gamma1,gamma2,pi1,pi2,+1);
                [sint1sl] = Sint(rho1,u1,p1,beta1,alpha1,usl,gamma1,gamma2,pi1,pi2,-1);
                sL(i-1) = sintsl2+sintCD+sint34;
                sR(i)   = sint1sl;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = Fsr;
                [sint12]  = Sint(rho1,u1,p1,beta1,alpha1,u23,gamma1,gamma2,pi1,pi2,-1);
                [sintCD]  = SintCD(u23,p23,alpha3,alpha2);
                [sint3sr] = Sint(rho3,u23,p23,beta3,alpha3,usr,gamma1,gamma2,pi1,pi2,+1);
                [sintsr4] = Sint(rhosr,usr,psr,betasr,alphasr,u4,gamma1,gamma2,pi1,pi2,+1);
                sL(i-1) = sint12+sintCD+sint3sr;
                sR(i)   = sintsr4;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = Fr;
                [sint12] = Sint(rho1,u1,p1,beta1,alpha1,u23,gamma1,gamma2,pi1,pi2,-1);
                [sintCD] = SintCD(u23,p23,alpha3,alpha2);
                [sint34] = Sint(rho3,u23,p23,beta3,alpha3,u4,gamma1,gamma2,pi1,pi2,+1);
                sL(i-1) = sint12+sintCD+sint34;
                sR(i)   = 0;
            end
        end %flux choice

    end %if L==R
end %for i loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two-Fluid Osher                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho2,rho3,u23,p23,beta2,beta3,alpha2,alpha3]=TFosher(rhol,ul,pl,betal,alphal,al,rhor,ur,pr,betar,alphar,ar,gamma1,gamma2,pi1,pi2)

NRK4   = 10;

%----------------------------%
% Single-fluid initial guess %
%----------------------------%
gamma  = (1/(alphal/gamma1+(1-alphal)/gamma2)+1/(alphar/gamma1+(1-alphar)/gamma2))/2;
% gamma  = (1/(alphal/(gamma1-1)+(1-alphal)/(gamma2-1))+1/(alphar/(gamma1-1)+(1-alphar)/(gamma2-1)))/2+1;
pil = (alphal*gamma1*pi1/(gamma1-1)+(1-alphal)*gamma2*pi2/(gamma2-1))/(alphal/(gamma1-1)+(1-alphal)/(gamma-1));
pir = (alphar*gamma1*pi1/(gamma1-1)+(1-alphar)*gamma2*pi2/(gamma2-1))/(alphar/(gamma1-1)+(1-alphar)/(gamma-1));
uguess = (ul/al*((pl+pil)/(pr+pir))^((gamma-1)/(2*gamma))+ur/ar+2/(gamma-1)*(((pl+pil)/(pr+pir))^((gamma-1)/(2*gamma))-1))/...
         (((pl+pil)/(pr+pir))^((gamma-1)/(2*gamma))/al+1/ar);
% keyboard
%--------------------------%
%                          %
%--------------------------%
diffp = 1;
m = 0;

while diffp>1e-1
    m=m+1;

    %--------------------------%
    % Left wave                %
    %--------------------------%
    if m == 1
        uf = 1.0*uguess;
    end

    rho   = rhol;
    p     = pl;
    beta  = betal;
    alpha = alphal;
    
    diffrhopalpha = 1;
    while diffrhopalpha>1e-1
        a = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
        b = alpha*(1-alpha)*(gamma1*(p+pi1)-gamma2*(p+pi2))/((1-alpha)*gamma1*(p+pi1)+alpha*gamma2*(p+pi2));
        rhon = rho+rho*uf/a;
        pn = p+rho*a*uf;
        alphan = alpha+b*uf/a;
        diffrhopalpha = max(abs([rhon-rho pn-p alphan-alpha]))
        rho = rhon; p=pn; alpha = alphan;
    end

    rho2   = rho;
    pf2    = p;
    beta2  = beta;
    alpha2 = alpha;
% keyboard
    %--------------------------%
    % Right wave               %
    %--------------------------%
    rho   = rhor;
    p     = pr;
    beta  = betar;
    alpha = alphar;

    diffrhopalpha = 1;
    while diffrhopalpha>1e-1
        a = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
        b = alpha*(1-alpha)*(gamma1*(p+pi1)-gamma2*(p+pi2))/((1-alpha)*gamma1*(p+pi1)+alpha*gamma2*(p+pi2));
        rhon = rho-rho*uf/a;
        pn = p-rho*a*uf;
        alphan = alpha-b*uf/a;
        diffrhopalpha = max(abs([rhon-rho pn-p alphan-alpha]));
        rho = rhon; p=pn; alpha = alphan;
    end

    rho3   = rho;
    pf3    = p;
    beta3  = beta;
    alpha3 = alpha;
% keyboard
    dp    = pf3 - pf2;
    diffp = abs(dp);

    if m==1
        ufold = uf;
        dpold = dp;
%         keyboard
    else
        ufnew = uf - dp*(uf - ufold+eps)/(dp - dpold+eps);
        ufold = uf;
%         keyboard
        uf    = ufnew;
        dpold = dp;
    end

end %while loop
% keyboard
u23 = uf;
p23 = pf2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left sonic point                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rhosl,usl,psl,betasl,alphasl]=TFO_sonic_L(u4,p4,a4,rho1,u1,p1,beta1,alpha1,a1,gamma1,gamma2,pi1,pi2)

%Runga-Kutta Solver
NRK4   = 10;
gamma  = alpha1*gamma1+(1-alpha1)*gamma2;
pim = (alpha1*gamma1*pi1/(gamma1-1)+(1-alpha1)*gamma2*pi2/(gamma2-1))/(alpha1/(gamma1-1)+(1-alpha1)/(gamma-1));
uguess = (u1/a1*((p1+pim)/(p4+pim))^((gamma-1)/(2*gamma))+u4/a4+2/(gamma-1)*(((p1+pim)/(p4+pim))^((gamma-1)/(2*gamma))-1))/...
         (((p1+pim)/(p4+pim))^((gamma-1)/(2*gamma))/a1+1/a4);

%--------------------------%
% Runga-Kutta fourth order %
%--------------------------%
diffu = 1;
m = 0;

while diffu>1e-3
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
    dul   = (usl - u1)/NRK4;

    for k=1:NRK4
        [dp_k1,drho_k1,dalpha_k1]=RK4_p(rho,p,alpha,0,0,0,gamma1,gamma2,pi1,pi2,-dul);
        [dp_k2,drho_k2,dalpha_k2]=RK4_p(rho,p,alpha,drho_k1,dp_k1,dalpha_k1,gamma1,gamma2,pi1,pi2,-dul);
        [dp_k3,drho_k3,dalpha_k3]=RK4_p(rho,p,alpha,drho_k2,dp_k2,dalpha_k2,gamma1,gamma2,pi1,pi2,-dul);
        [dp_k4,drho_k4,dalpha_k4]=RK4_p(rho,p,alpha,drho_k3,dp_k3,dalpha_k3,gamma1,gamma2,pi1,pi2,-dul);
        p = p + (dp_k1 + 2*dp_k2 + 2*dp_k3 + dp_k4)/6;
        rho = rho + (drho_k1 + 2*drho_k2 + 2*drho_k3 + drho_k4)/6;
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
function [rhosr,usr,psr,betasr,alphasr]=TFO_sonic_R(u1,p1,a1,rho4,u4,p4,beta4,alpha4,a4,gamma1,gamma2,pi1,pi2)

%Runga-Kutta Solver
NRK4   = 10;
gamma  = alpha4*gamma1+(1-alpha4)*gamma2;
pim = (alpha4*gamma1*pi1/(gamma1-1)+(1-alpha4)*gamma2*pi2/(gamma2-1))/(alpha4/(gamma1-1)+(1-alpha4)/(gamma-1));
uguess = (u1/a1*((p1+pim)/(p4+pim))^((gamma-1)/(2*gamma))+u4/a4+2/(gamma-1)*(((p1+pim)/(p4+pim))^((gamma-1)/(2*gamma))-1))/...
         (((p1+pim)/(p4+pim))^((gamma-1)/(2*gamma))/a1+1/a4);

%--------------------------%
% Runga-Kutta fourth order %
%--------------------------%
diffu = 1;
m = 0;

while diffu>1e-3
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
    dur   = (usr - u4)/NRK4;

    for k=1:NRK4
        [dp_k1,drho_k1,dalpha_k1]=RK4_p(rho,p,alpha,0,0,0,gamma1,gamma2,pi1,pi2,dur);
        [dp_k2,drho_k2,dalpha_k2]=RK4_p(rho,p,alpha,drho_k1,dp_k1,dalpha_k1,gamma1,gamma2,pi1,pi2,dur);
        [dp_k3,drho_k3,dalpha_k3]=RK4_p(rho,p,alpha,drho_k2,dp_k2,dalpha_k2,gamma1,gamma2,pi1,pi2,dur);
        [dp_k4,drho_k4,dalpha_k4]=RK4_p(rho,p,alpha,drho_k3,dp_k3,dalpha_k3,gamma1,gamma2,pi1,pi2,dur);
        p = p + (dp_k1 + 2*dp_k2 + 2*dp_k3 + dp_k4)/6;
        rho = rho + (drho_k1 + 2*drho_k2 + 2*drho_k3 + drho_k4)/6;
        alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
    end

    psr     = p;
    alphasr = alpha;
    betasr  = beta;
    rhosr   = rho;
    asr = sqrt(1/(alphasr/(gamma1*(psr+pi1))+(1-alphasr)/(gamma2*(psr+pi2)))/rhosr);

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