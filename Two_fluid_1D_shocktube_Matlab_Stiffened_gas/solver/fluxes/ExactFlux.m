%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exact flux solver                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,sL,sR,Lambda] = ExactFlux(N,Wf,gamma1,gamma2,pi1,pi2,Lambda)

% eps = 1e-12;
% speeds = zeros(N+1,5);
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

    % State in beide cellen hetzelfde
    if abs(rho1-rho4)<1e-8 && abs(u1-u4)<1e-8 && abs(p1-p4)<1e-8 && abs(beta1-beta4)<1e-8 && abs(alpha1-alpha4)<1e-8
    
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
    elseif abs(p1-p4)<1e-8 && abs(u1)<1e-8 && abs(u4)<1e-8

        F(:,i) = [0;p1;0;0;0];

        sL(i-1) = 0;
        sR(i)   = 0;
        
        
    % Two-fluid problem
    else

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
        % Runga-Kutta fourth order %
        %--------------------------%
        diffp = 1;
        m = 0;

        while diffp>1e-3
            m=m+1;

            if m == 1
                uf = 1.0*uguess;
            elseif m==2
                uf = 0.8*uguess;
            end

            if p>=p1%left == 1 %shock
                [rho2,pf2,beta2,alpha2]=TFshock(rhol,ul,pl,betal,alphal,uf);
            elseif p<p1 %left == 2 % expansion
                [rho2,pf2,beta2,alpha2]=TFexpansion(rhol,ul,pl,betal,alphal,uf);
            end

            if p>=p4 %right == 1 %shock
                [rho3,pf3,beta3,alpha3]=TFshock(rhor,ur,pr,betar,alphar,uf);
            elseif p<p4 %right == 2 %expansion
                [rho3,pf3,beta3,alpha3]=TFexpansion(rhor,ur,pr,betar,alphar,uf);
            end

            dp    = pf3 - pf2;
            diffp = abs(dp);

            if m==1
                ufold = uf;
                dpold = dp;

            else
                ufnew = uf - dp*(uf - ufold+eps)/(dp - dpold+eps);
                ufold = uf;

                uf    = ufnew;
                dpold = dp;
            end

        end %while loop

        u23 = uf;
        p23 = pf2;

        %speed of sound
        a2 = sqrt(1/(alpha2/(gamma1*(p23+pi1))+(1-alpha2)/(gamma2*(p23+pi2)))/rho2);
        a3 = sqrt(1/(alpha3/(gamma1*(p23+pi1))+(1-alpha3)/(gamma2*(p23+pi2)))/rho3);

        lambda0 = u23;

        if p23>p1
            gamma = alpha1/gamma1+(1-alpha1)/gamma2;
            lambda1 = u1-a1*sqrt(1+(gamma+1)/(2*gamma)*(p23/p1-1));
            lambda2 = lambda1;
        elseif p23<=p1
            lambda1 = u1-a1;
            lambda2 = u23-a2;
        end

        if p23>p4
            gamma = alpha4/gamma1+(1-alpha4)/gamma2;
            lambda4 = u4+a4*sqrt(1+(gamma+1)/(2*gamma)*(p23/p4-1));
            lambda3 = lambda4;
        elseif p23<=p4
            lambda3 = u23+a3;
            lambda4 = u4+a4;
        end

        if lambda1>=0
            rhof = rho1;
            uf   = u1;
            pf   = p1;
            betaf = beta1;
            alphaf = alpha1;
        elseif lambda1<0 && lambda2>=0 %expansion wave
            [rho3,pf3,beta3,alpha3]=TFexpansion_sonic_R(rho3,u3,p3,beta3,alpha3,uf);
            uf   = (gamma-1)/(gamma+1)*u1+2/(gamma+1)*a1;
            af   = u23;
            pf   = p1*(af/a1)^((2*gamma)/(gamma-1));
            rhof = rho1*(af/a1)^(2/(gamma-1));
        elseif lambda2<0 && lambda0>=0
            rhof = rho2;
            uf   = u23;
            pf   = p23;
        elseif lambda0<0 && lambda3>=0
            rhof = rho3;
            uf   = u23;
            pf   = p23;
        elseif lambda3<0 && lambda4>=0
            uf   = (gamma-1)/(gamma+1)*u4-2/(gamma+1)*a4;
            af   = -u23;
            pf   = p4*(af/a4)^(2*gamma/(gamma-1));
            rhof = rho4*(af/a4)^(2/(gamma-1));
        elseif lambda4<0
            rhof = rho4;
            uf   = u23;
            pf   = p23;
        end

        F(1,i) = rhof*uf;
        F(2,i) = rhof*uf^2+pf;
        F(3,i) = uf*pf*gamma/(gamma-1)+1/2*rhof*uf^3;

        Lambda  = max(abs([lambda0 lambda1 lambda2 lambda3 lambda4 Lambda]));