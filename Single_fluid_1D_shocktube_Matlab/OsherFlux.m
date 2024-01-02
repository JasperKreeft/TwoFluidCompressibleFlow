% Osher's flux solver

function [F,Lambda] = OsherFlux(N,Wf,gamma,Lambda)

eps = 1e-12;

F = zeros(3,N+1);

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
    if abs(rho1-rho4)+abs(u1-u4)+abs(p1-p4)<1e-12

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
        % calculate Osher's solution                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Fina1 state va1ues
        z    = (gamma-1)/(2*gamma);
        p23  = ((a1+a4-(gamma-1)/2*(u4-u1))/(a1/p1^z+a4/p4^z))^(1/z);
        u23  = ((p1/p4)^z*u1/a1+u4/a4+2/(gamma-1)*((p1/p4)^z-1))/(1/a1*(p1/p4)^z+1/a4);
        rho2 = rho1*(p23/p1)^(1/gamma);
        rho3 = rho4*(p23/p4)^(1/gamma);

        % Speed of sound
        a2 = sqrt(gamma*p23/rho2);
        a3 = sqrt(gamma*p23/rho3);

        % Sonic point va1ues
        usl   = (gamma-1)/(gamma+1)*u1+2*a1/(gamma+1);
        asl   = usl;
        rhosl = rho1*(asl/a1)^(2/(gamma-1));
        psl   = p1*(rhosl/rho1)^gamma;

        usr   = (gamma-1)/(gamma+1)*u4-2*a4/(gamma+1);
        asr   = -usr;
        rhosr = rho4*(asr/a4)^(2/(gamma-1));
        psr   = p4*(rhosr/rho4)^gamma;

        F1  = [   rho1*u1;    rho1*u1^2+p1;     (gamma/(gamma-1))*u1*p1+0.5*rho1*u1^3];
        F2  = [  rho2*u23;  rho2*u23^2+p23;  (gamma/(gamma-1))*u23*p23+0.5*rho2*u23^3];
        F3  = [  rho3*u23;  rho3*u23^2+p23;  (gamma/(gamma-1))*u23*p23+0.5*rho3*u23^3];
        F4  = [   rho4*u4;    rho4*u4^2+p4;     (gamma/(gamma-1))*u4*p4+0.5*rho4*u4^3];
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

    end %if L==R or supersonic

end %for i loop