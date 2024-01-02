%% flux solver
% Based on Book of Toro

function [F,Lambda] = ExactFlux(N,Wf,gamma,Lambda)

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

        z  = (gamma-1)/(2*gamma);

        % Initial guess

        % % Two-rarefaction approximation
        p    = ((a1+a4-(gamma-1)/2*(u4-u1))/(a1/p1^z+a4/p4^z))^z;
        % 
        % % Linearised solution based on primitive variables
        % p_pv = 1/2*(p1+p4)-1/8*(u4-u1)*(rho1+rho4)*(a1+a4);
        % p    = max(1e-6,p_pv);
        % 
        % % Two-shock approximation
        % p_pv = 1/2*(p1+p4)-1/8*(u4-u1)*(rho1+rho4)*(a1+a4);
        % p0   = max([1e-6 p_pv]);
        % g1   = sqrt((2/((gamma+1)*rho1))/(p0+(gamma-1)/(gamma+1)*p1));
        % g4   = sqrt((2/((gamma+1)*rho4))/(p0+(gamma-1)/(gamma+1)*p4));
        % p_ts = (g1*p1+g4*p4-(u4-u1))/(g1+g4);
        % p    = max([1e-6 p_ts]);

        % Arithmetic mean
        % p    = 1/2*(p1+p4);

        fout = 1;
        while fout>1e-8

            if p>p1      % shock left
                f1 = (p-p1)*sqrt(2/((gamma+1)*rho1)/(p+(gamma-1)/(gamma+1)*p1));
            elseif p<=p1 % expansion left
                f1 = 2*a1/(gamma-1)*((p/p1)^z-1);
            end

            if p>p4      % shock right
                f4 = (p-p4)*sqrt(2/((gamma+1)*rho4)/(p+(gamma-1)/(gamma+1)*p4));
            elseif p<=p4 % expansion right
                f4 = 2*a4/(gamma-1)*((p/p4)^z-1);
            end

            f = f1+f4+(u4-u1);

            if p>p1
                f1_p = (1-(p-p1)/(2*(p+(gamma-1)/(gamma+1)*p1)))*sqrt((2/((gamma+1)*rho1))/(p+(gamma-1)/(gamma+1)*p1));
            elseif p<=p1
                f1_p = 1/(rho1*a1)*(p/p1)^(-(gamma+1)/(2*gamma));
            end

            if p>p4
                f4_p = (1-(p-p4)/(2*(p+(gamma-1)/(gamma+1)*p4)))*sqrt((2/((gamma+1)*rho4))/(p+(gamma-1)/(gamma+1)*p4));
            elseif p<=p4
                f4_p = 1/(rho4*a4)*(p/p4)^(-(gamma+1)/(2*gamma));
            end

            f_p = f1_p+f4_p;

            p_new = p-f/f_p;

            fout = abs(p_new-p)/(1/2*(p_new+p));

            p = max(p_new,1e-10);

        end

        p23 = p;

        u23 = 1/2*(u1+u4)+1/2*(f4-f1);

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

    end %if L==R
end %for i loop