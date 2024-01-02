% Osher's flux solver for the y-direction

function [G,Lambda] = OsherFlux_y(NX,NY,Wfy,gamma,Lambda)

eps = 1e-12;

G  = zeros(NY+1,NX,4);

for i=2:NY
    for j=1:NX

        rho1 = Wfy(2*(i-1),j,1);
        u1   = Wfy(2*(i-1),j,2);
        v1   = Wfy(2*(i-1),j,3);
        p1   = Wfy(2*(i-1),j,4);
        a1   = sqrt(gamma*p1/rho1);
        M1   = v1/a1;

        rho4 = Wfy(2*i-1,j,1);
        u4   = Wfy(2*i-1,j,2);
        v4   = Wfy(2*i-1,j,3);
        p4   = Wfy(2*i-1,j,4);
        a4   = sqrt(gamma*p4/rho4);
        M4   = v4/a4;

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
        
        % Cavity test
        if (v1+2/(gamma-1)*a1)<(v4-2/(gamma-1)*a4) 
            disp('vacuum is created')
            pause
        end

        % State in beide cellen hetzelfde
        if abs(rho1-rho4)<1e-12 && abs(v1-v4)<1e-12 && abs(p1-p4)<1e-12

            G(i,j,:) = [rho1*v1;
                        rho1*v1*u1;
                        rho1*v1^2+p1;
                        v1*p1*gamma/(gamma-1)+1/2*rho1*v1*(u1^2+v1^2)];

            Lambda  = max(abs([v1 v1-a1 v4+a4 Lambda]));

        % Supersonic tests
        elseif M1>1 && M4>1
    
            G(i,j,:) = [rho1*v1;
                        rho1*v1*u1;
                        rho1*v1^2+p1;
                        v1*p1*gamma/(gamma-1)+1/2*rho1*v1*(u1^2+v1^2)];

            Lambda  = max(abs([v1 v1-a1 v4+a4 Lambda]));

        elseif M1<-1 && M4<-1
    
            G(i,j,:) = [rho4*v4;
                        rho4*v4*u4;
                        rho4*v4^2+p4;
                        v4*p4*gamma/(gamma-1)+1/2*rho4*v4*(u4^2+v4^2)];

            Lambda  = max(abs([v1 v1-a1 v4+a4 Lambda]));

        else
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % calculate Osher's solution                         %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Fina1 state va1ues
            z    = (gamma-1)/(2*gamma);
            p23  = ((a1+a4-(gamma-1)/2*(v4-v1))/(a1/p1^z+a4/p4^z))^(1/z);
            v23  = ((p1/p4)^z*v1/a1+v4/a4+2/(gamma-1)*((p1/p4)^z-1))/(1/a1*(p1/p4)^z+1/a4);
            rho2 = rho1*(p23/p1)^(1/gamma);
            rho3 = rho4*(p23/p4)^(1/gamma);
            u2   = u1;
            u3   = u4;

            % Speed of sound
            a2 = sqrt(gamma*p23/rho2);
            a3 = sqrt(gamma*p23/rho3);

            % Sonic point va1ues
            vsl   = (gamma-1)/(gamma+1)*v1+2*a1/(gamma+1);
            asl   = vsl;
            rhosl = rho1*(asl/a1)^(2/(gamma-1));
            psl   = p1*(rhosl/rho1)^gamma;
            usl   = u1;

            vsr   = (gamma-1)/(gamma+1)*v4-2*a4/(gamma+1);
            asr   = -vsr;
            rhosr = rho4*(asr/a4)^(2/(gamma-1));
            psr   = p4*(rhosr/rho4)^gamma;
            usr   = u4;

            G1  = [   rho1*v1;    rho1*v1*u1;    rho1*v1^2+p1;       (gamma/(gamma-1))*v1*p1+0.5*rho1*v1*(u1^2+v1^2)];
            G2  = [  rho2*v23;   rho2*v23*u2;  rho2*v23^2+p23;   (gamma/(gamma-1))*v23*p23+0.5*rho2*v23*(u2^2+v23^2)];
            G3  = [  rho3*v23;   rho3*v23*u3;  rho3*v23^2+p23;   (gamma/(gamma-1))*v23*p23+0.5*rho3*v23*(u3^2+v23^2)];
            G4  = [   rho4*v4;    rho4*v4*u4;    rho4*v4^2+p4;       (gamma/(gamma-1))*v4*p4+0.5*rho4*v4*(u4^2+v4^2)];
            Gsl = [rhosl*vsl; rhosl*vsl*usl; rhosl*vsl^2+psl; (gamma/(gamma-1))*vsl*psl+0.5*rhosl*vsl*(usl^2+vsl^2)];
            Gsr = [rhosr*vsr; rhosr*vsr*usr; rhosr*vsr^2+psr; (gamma/(gamma-1))*vsr*psr+0.5*rhosr*vsr*(usr^2+vsr^2)];

            lambda0 = v23;
            lambda1 = v1-a1;
            lambda2 = v23-a2;
            lambda3 = v23+a3;
            lambda4 = v4+a4;
            Lambda  = max(abs([lambda0 lambda1 lambda2 lambda3 lambda4 Lambda]));

            if lambda0 >= 0 && lambda2>=0
                if lambda1>=0 && lambda4>=0
                    G(i,j,:) = G1;
                elseif lambda1>=0 && lambda4<=0
                    G(i,j,:) = G1+G4-Gsr;
                elseif lambda1<=0 && lambda4>=0
                    G(i,j,:) = Gsl;
                elseif lambda1<=0 && lambda4<=0
                    G(i,j,:) = Gsl-Gsr+G4;
                end

            elseif  lambda0 >= 0 && lambda2<=0
                if lambda1>=0 && lambda4>=0
                    G(i,j,:) = G1-Gsl+G2;
                elseif lambda1>=0 && lambda4<=0
                    G(i,j,:) = G1-Gsl+G2-Gsr+G4;
                elseif lambda1<=0 && lambda4>=0
                    G(i,j,:) = G2;
                elseif lambda1<=0 && lambda4<=0
                    G(i,j,:) = G4+G2-Gsr;
                end
    
            elseif  lambda0 <= 0 && lambda3>=0
                if lambda1>=0 && lambda4>=0
                    G(i,j,:) = G1-Gsl+G3;
                elseif lambda1>=0 && lambda4<=0
                    G(i,j,:) = G1-Gsl+G3-Gsr+G4;
                elseif lambda1<=0 && lambda4>=0
                    G(i,j,:) = G3;
                elseif lambda1<=0 && lambda4<=0
                    G(i,j,:) = G4+G3-Gsr;
                end

            elseif  lambda0 <= 0 && lambda3<=0
                if lambda1>=0 && lambda4>=0
                    G(i,j,:) = G1-Gsl+Gsr;
                elseif lambda1>=0 && lambda4<=0
                    G(i,j,:) = G1-Gsl+G4;
                elseif lambda1<=0 && lambda4>=0
                    G(i,j,:) = Gsr;
                elseif lambda1<=0 && lambda4<=0
                    G(i,j,:) = G4;
                end    

            end %flux choice
        end %if L==R or supersonic
    end %for j loop
end %for i loop