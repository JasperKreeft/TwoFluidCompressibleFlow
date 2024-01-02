%% Linear Osher's flux solver

function [F,Lambda] = LinearOsherFlux(N,Wf,gamma,Lambda)

eps = 1e-12;

F  = zeros(3,N+1);

% i = j - 1/2
for i  = 2:N

    rho1   = Wf(1,2*i-2);
    u1     = Wf(2,2*i-2);
    p1     = Wf(3,2*i-2);
    a1     = sqrt(gamma*p1/rho1);
    M1     = u1/a1;

    rho4   = Wf(1,2*i-1);
    u4     = Wf(2,2*i-1);
    p4     = Wf(3,2*i-1);
    a4     = sqrt(gamma*p4/rho4);
    M4     = u4/a4;
    
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

    % cavity test
    if (u1+2/(gamma-1)*a1)<(u4-2/(gamma-1)*a4)
        disp('vacuum is created')
        pause
    end

    % State in beide cellen hetzelfde
    if abs(rho1-rho4)<1e-12 && abs(u1-u4)<1e-12 && abs(p1-p4)<1e-12

        F(:,i) = [rho1*u1;
                  rho1*u1^2+p1;
                  u1*p1*gamma/(gamma-1)+1/2*rho1*u1^3];

        Lambda  = max(abs([u1 u1-a1 u4+a4 Lambda]));


    % supersonic tests
    elseif M1>1 && M4>1 

        F(:,i) = [rho1*u1;
                  rho1*u1^2+p1;
                  u1*p1*gamma/(gamma-1)+1/2*rho1*u1^3];
 
        Lambda  = max(abs([u1 u1-a1 u4+a4 Lambda]));

    elseif M1<-1 && M4<-1

        F(:,i) = [rho4*u4;
                  rho4*u4^2+p4;
                  u4*p4*gamma/(gamma-1)+1/2*rho4*u4^3];

        Lambda  = max(abs([u1 u1-a1 u4+a4 Lambda]));

    else  

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate Linear Osher's solution                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        %Final state values
        u23 = (rho1*a1*u1+rho4*a4*u4+(p1-p4))/(rho1*a1+rho4*a4);
        p23 = (rho4*a4*p1+rho1*a1*p4+rho1*a1*rho4*a4*(u1-u4))/(rho1*a1+rho4*a4);
        rho2 = rho1+(p23-p1)/a1^2;
        rho3 = rho4+(p23-p4)/a4^2;

        %speed of sound
        a2 = sqrt(gamma*p23/rho2);
        a3 = sqrt(gamma*p23/rho3);

        %Sonic point va1ues
        [rhosl,usl,psl] = linosher_sonic_L(rho1,u1,p1,a1,gamma);
        [rhosr,usr,psr] = linosher_sonic_R(rho4,u4,p4,a4,gamma);

        F1 = [rho1*u1; rho1*u1^2+p1; (gamma/(gamma-1))*u1*p1+0.5*rho1*u1^3];
        F2 = [rho2*u23; rho2*u23^2+p23; (gamma/(gamma-1))*u23*p23+0.5*rho2*u23^3];
        F3 = [rho3*u23; rho3*u23^2+p23; (gamma/(gamma-1))*u23*p23+0.5*rho3*u23^3];
        F4 = [rho4*u4; rho4*u4^2+p4; (gamma/(gamma-1))*u4*p4+0.5*rho4*u4^3];
        Fsl = [rhosl*usl; rhosl*usl^2+psl; (gamma/(gamma-1))*usl*psl+0.5*rhosl*usl^3];
        Fsr = [rhosr*usr; rhosr*usr^2+psr; (gamma/(gamma-1))*usr*psr+0.5*rhosr*usr^3];

        lambda0 = u23;
        lambda1 = u1-a1;
        lambda2 = u23-a2;
        lambda3 = u23+a3;
        lambda4 = u4+a4;
        Lambda  = max(abs([lambda0 lambda1 lambda2 lambda3 lambda4 Lambda]));

        if lambda0 >= 0 && lambda2>=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = F1;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = F1+F4-Fsr;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = Fsl;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = Fsl-Fsr+F4;
            end

        elseif  lambda0 >= 0 && lambda2<=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = F1-Fsl+F2;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = F1-Fsl+F2-Fsr+F4;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = F2;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = F4+F2-Fsr;
            end
    
        elseif  lambda0 <= 0 && lambda3>=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = F1-Fsl+F3;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = F1-Fsl+F3-Fsr+F4;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = F3;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = F4+F3-Fsr;
            end
    
        elseif  lambda0 <= 0 && lambda3<=0
            if lambda1>=0 && lambda4>=0
                F(:,i) = F1-Fsl+Fsr;
            elseif lambda1>=0 && lambda4<=0
                F(:,i) = F1-Fsl+F4;
            elseif lambda1<=0 && lambda4>=0
                F(:,i) = Fsr;
            elseif lambda1<=0 && lambda4<=0
                F(:,i) = F4;
            end    
    
        end %flux choice

    end %if L==R
    
end %for i loop

%% Left sonic point

function [rhosl,usl,psl]=linosher_sonic_L(rho1,u1,p1,a1,gamma)

%-----------------------------------------------------------------------%
%                         Newton-Raphson solver                         %
%-----------------------------------------------------------------------%

uguess = a1;

diffu = 1;
m = 0;

while (diffu>1e-8)
    m=m+1;

    %-----------------------------------------------------------------------%
    %                               LEFT WAVE                               %
    %-----------------------------------------------------------------------%
    if m==1
        usl = uguess;
    elseif m==2
        usl = 0.8*uguess;
    end

    psl   = p1-rho1*a1*(usl-u1);
    rhosl = rho1-rho1/a1*(usl-u1);
    asl   = sqrt(gamma*psl/rhosl);

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

%% Right sonic point

function [rhosr,usr,psr]=linosher_sonic_R(rho4,u4,p4,a4,gamma)

%-----------------------------------------------------------------------%
%                         Newton-Raphson solver                         %
%-----------------------------------------------------------------------%

uguess = -a4;

diffu = 1;
m = 0;

while diffu>1e-8
    m=m+1;

    %-----------------------------------------------------------------------%
    %                               RIGHT WAVE                               %
    %-----------------------------------------------------------------------%
    if m==1
        usr = uguess;
    elseif m==2
        usr = 0.8*uguess;
    end

    psr   = p4+rho4*a4*(usr-u4);
    rhosr = rho4+rho4/a4*(usr-u4);
    asr   = sqrt(gamma*psr/rhosr);

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