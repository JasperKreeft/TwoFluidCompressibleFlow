% Exact flux solver
% Based on Gasdynamics

function [F,Lambda] = ExactFlux2(N,Wf,gamma,Lambda)

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
    if abs(rho1-rho4)<1e-12 && abs(u1-u4)<1e-12 && abs(p1-p4)<1e-12  %state in beide cellen hetzelfde

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
        % calculate exact solution                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %initial guess final state pressure
        z    = (gamma-1)/(2*gamma);
        p01  = ((a1+a4-(gamma-1)/2*(u4-u1))/(a1/p1^z+a4/p4^z))^(1/z);

        if p1>p01 && p4>p01 %2 expansion waves (exact solution possible)

            p23 = p01;

        else

            p02 = (rho1*a1*p4+rho4*a4*p1-rho1*a1*rho4*a4*(u4-u1))/(rho1*a1+rho4*a4);
            % p02 = .5*(p1+p4);

            j = 0;
            pnew = 1;
            pold = 0;
            fout = 1;
            jmax = 20;

            while fout>1e-8
                j = j + 1;

                if j==1
                    p = max(eps,p02);
                elseif j==2
                    p = max(eps,p01);
                else
                    p = max(eps,pnew);
                end

                if p>=p1%left == 1 %shock
                    m1 = rho1*a1*sqrt(1+(gamma+1)/(2*gamma)*(p/p1-1));
                elseif p<p1 %left == 2 % expansion
                    m1 = rho1*a1*(gamma-1)/(2*gamma)*(1-p/p1)/(1-(p/p1)^z+eps);
                end

                if p>=p4 %right == 1 %shock
                    m4 = rho4*a4*sqrt(1+(gamma+1)/(2*gamma)*(p/p4-1));
                elseif p<p4 %right == 2 %expansion
                    m4 = rho4*a4*z*(1-p/p4)/(1-(p/p4)^z+eps);
                end

                %pressure of final state
                G = (m1*p4+m4*p1-m1*m4*(u4-u1))/(m1+m4+eps) - p;

                if j>1
                    pnew = p - G*(p-pold+eps)/(G-Gold+eps);
                    pnew = max(pnew,eps);
                    fout = abs(pnew-p)/(1/2*(pnew+p)+eps);
                end

                %save data
                Gold  = G;
                pold = p;

                if j>=jmax
                    disp('not converged')
                end
            end % while

            p23 = max(pnew,1e-10);

        end % 2 expansions

        %velocity of final state
        if p23>=p1 %left == 1 %shock
            m1 = rho1*a1*sqrt(1+(gamma+1)/(2*gamma)*(p23-p1)/p1);
        elseif p23<p1 %left == 2 %expansion
            m1 = rho1*a1*(gamma-1)/(2*gamma)*(1-p23/p1)/(1-(p23/p1)^z+eps);
        end

        if p23>=p4 %right == 1 % shock 
            m4 = rho4*a4*sqrt(1+(gamma+1)/(2*gamma)*(p23-p4)/p4);
        elseif p23<p4 %right == 2 %expansion
            m4 = rho4*a4*z*(1-p23/p4)/(1-(p23/p4)^z+eps);
        end

        u23 = (m1*u1+m4*u4-(p4-p1))/(m1+m4+eps);

        %density
        if p23>=p1 %left == 1 %shock
            rho2 = rho1*(1+(gamma+1)/(gamma-1)*p23/p1)/((gamma+1)/(gamma-1)+p23/p1);
        elseif p23<p1 %left == 2 %expansion
            rho2 = rho1*(p23/p1)^(1/gamma);
        end

        if p23>=p4 %right ==1 %shock
            rho3 = rho4*(1+(gamma+1)/(gamma-1)*p23/p4)/((gamma+1)/(gamma-1)+p23/p4);
        elseif p23<p4 %right == 2 %expansion
            rho3 = rho4*(p23/p4)^(1/gamma);
        end

        %speed of sound
        a2 = sqrt(gamma*p23/rho2);
        a3 = sqrt(gamma*p23/rho3);

        lambda0 = u23;

        if p23>p1
            lambda1 = u1-a1*sqrt(1+(gamma+1)/(2*gamma)*(p23/p1-1));
            lambda2 = lambda1;
        elseif p23<=p1
            lambda1 = u1-a1;
            lambda2 = u23-a2;
        end

        if p23>p4
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
        elseif lambda1<0 && lambda2>=0
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
            uf   = u4;
            pf   = p4;
        end

        F(1,i) = rhof*uf;
        F(2,i) = rhof*uf^2+pf;
        F(3,i) = uf*pf*gamma/(gamma-1)+1/2*rhof*uf^3;

        Lambda  = max(abs([lambda0 lambda1 lambda2 lambda3 lambda4 Lambda]));
        
    end %if L==R or supersonic

end %for i loop