% Osher's flux solver for the x-direction

function [F,Lambda] = OsherFlux_x(NX,NY,Wfx,gamma,Lambda)

eps = 1e-12;

F  = zeros(NY,NX+1,4);

for i=1:NY
    for j=2:NX

        rho1 = Wfx(i,2*(j-1),1);
        u1   = Wfx(i,2*(j-1),2);
        v1   = Wfx(i,2*(j-1),3);
        p1   = Wfx(i,2*(j-1),4);
        a1   = sqrt(gamma*p1/rho1);
        M1   = u1/a1;

        rho4 = Wfx(i,2*j-1,1);
        u4   = Wfx(i,2*j-1,2);
        v4   = Wfx(i,2*j-1,3);
        p4   = Wfx(i,2*j-1,4);
        a4   = sqrt(gamma*p4/rho4);
        M4   = u4/a4;

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
    
        % state in beide cellen hetzelfde
        if abs(rho1-rho4)<1e-12 && abs(u1-u4)<1e-12 && abs(p1-p4)<1e-12

            F(i,j,:) = [rho1*u1;
                        rho1*u1^2+p1;
                        rho1*u1*v1;
                        u1*p1*gamma/(gamma-1)+1/2*rho1*u1*(u1^2+v1^2)];

            Lambda  = max(abs([u1 u1-a1 u4+a4 Lambda]));

        % supersonic tests
        elseif M1>1 && M4>1

            F(i,j,:) = [rho1*u1;
                        rho1*u1^2+p1;
                        rho1*u1*v1;
                        u1*p1*gamma/(gamma-1)+1/2*rho1*u1*(u1^2+v1^2)];
      
            Lambda  = max(abs([u1 u1-a1 u4+a4 Lambda]));

        elseif M1<-1 && M4<-1
    
            F(i,j,:) = [rho4*u4;
                        rho4*u4*v4;
                        rho4*u4^2+p4;
                        u4*p4*gamma/(gamma-1)+1/2*rho4*u4*(u4^2+v4^2)];

            Lambda  = max(abs([u1 u1-a1 u4+a4 Lambda]));

        else

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % calculate Osher's solution                         %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Fina1 state va1ues
            z    = (gamma-1)/(2*gamma);
            p23   = ((a1+a4-(gamma-1)/2*(u4-u1))/(a1/p1^z+a4/p4^z))^(1/z);
            u23   = ((p1/p4)^z*u1/a1+u4/a4+2/(gamma-1)*((p1/p4)^z-1))/(1/a1*(p1/p4)^z+1/a4);
            rho2 = rho1*(p23/p1)^(1/gamma);
            rho3 = rho4*(p23/p4)^(1/gamma);
            v2   = v1;
            v3   = v4;

            % Speed of sound
            a2 = sqrt(gamma*p23/rho2);
            a3 = sqrt(gamma*p23/rho3);

            % Sonic point va1ues
            usl   = (gamma-1)/(gamma+1)*u1+2*a1/(gamma+1);
            asl   = usl;
            rhosl = rho1*(asl/a1)^(2/(gamma-1));
            psl   = p1*(rhosl/rho1)^gamma;
            vsl   = v1;

            usr   = (gamma-1)/(gamma+1)*u4-2*a4/(gamma+1);
            asr   = -usr;
            rhosr = rho4*(asr/a4)^(2/(gamma-1));
            psr   = p4*(rhosr/rho4)^gamma;
            vsr   = v4;

            F1  = [   rho1*u1;    rho1*u1^2+p1;    rho1*u1*v1;       (gamma/(gamma-1))*u1*p1+0.5*rho1*u1*(u1^2+v1^2)];
            F2  = [  rho2*u23;  rho2*u23^2+p23;   rho2*u23*v2;   (gamma/(gamma-1))*u23*p23+0.5*rho2*u23*(u23^2+v2^2)];
            F3  = [  rho3*u23;  rho3*u23^2+p23;   rho3*u23*v3;   (gamma/(gamma-1))*u23*p23+0.5*rho3*u23*(u23^2+v3^2)];
            F4  = [   rho4*u4;    rho4*u4^2+p4;    rho4*u4*v4;       (gamma/(gamma-1))*u4*p4+0.5*rho4*u4*(u4^2+v4^2)];
            Fsl = [rhosl*usl; rhosl*usl^2+psl; rhosl*usl*vsl; (gamma/(gamma-1))*usl*psl+0.5*rhosl*usl*(usl^2+vsl^2)];
            Fsr = [rhosr*usr; rhosr*usr^2+psr; rhosr*usr*vsr; (gamma/(gamma-1))*usr*psr+0.5*rhosr*usr*(usr^2+vsr^2)];

            lambda0 = u23;
            lambda1 = u1-a1;
            lambda2 = u23-a2;
            lambda3 = u23+a3;
            lambda4 = u4+a4;
            Lambda  = max(abs([lambda0 lambda1 lambda2 lambda3 lambda4 Lambda]));

            if lambda0 >= 0 && lambda2>=0
                if lambda1>=0 && lambda4>=0
                    F(i,j,:) = F1;
                elseif lambda1>=0 && lambda4<=0
                    F(i,j,:) = F1+F4-Fsr;
                elseif lambda1<=0 && lambda4>=0
                    F(i,j,:) = Fsl;
                elseif lambda1<=0 && lambda4<=0
                    F(i,j,:) = Fsl-Fsr+F4;
                end

            elseif  lambda0 >= 0 && lambda2<=0
                if lambda1>=0 && lambda4>=0
                    F(i,j,:) = F1-Fsl+F2;
                elseif lambda1>=0 && lambda4<=0
                    F(i,j,:) = F1-Fsl+F2-Fsr+F4;
                elseif lambda1<=0 && lambda4>=0
                    F(i,j,:) = F2;
                elseif lambda1<=0 && lambda4<=0
                    F(i,j,:) = F4+F2-Fsr;
                end
    
            elseif  lambda0 <= 0 && lambda3>=0
                if lambda1>=0 && lambda4>=0
                    F(i,j,:) = F1-Fsl+F3;
                elseif lambda1>=0 && lambda4<=0
                    F(i,j,:) = F1-Fsl+F3-Fsr+F4;
                elseif lambda1<=0 && lambda4>=0
                    F(i,j,:) = F3;
                elseif lambda1<=0 && lambda4<=0
                    F(i,j,:) = F4+F3-Fsr;
                end

            elseif  lambda0 <= 0 && lambda3<=0
                if lambda1>=0 && lambda4>=0
                    F(i,j,:) = F1-Fsl+Fsr;
                elseif lambda1>=0 && lambda4<=0
                    F(i,j,:) = F1-Fsl+F4;
                elseif lambda1<=0 && lambda4>=0
                    F(i,j,:) = Fsr;
                elseif lambda1<=0 && lambda4<=0
                    F(i,j,:) = F4;
                end

            end %flux choice
        end %if L==R or supersonic
    end %for j loop
end %for i loop