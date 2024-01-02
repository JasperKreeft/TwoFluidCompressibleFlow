% HLL flux solver

function [F,Lambda] = HLLFlux(N,Wf,gamma,Lambda)

eps = 1e-12;

F = zeros(3,N+1);

% i = j - 1/2
for i  = 2:N

    rho1   = Wf(1,2*i-2);
    u1     = Wf(2,2*i-2);
    p1     = Wf(3,2*i-2);
    a1     = sqrt(gamma*p1/rho1);

    rho4   = Wf(1,2*i-1);
    u4     = Wf(2,2*i-1);
    p4     = Wf(3,2*i-1);
    a4     = sqrt(gamma*p4/rho4);

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

    % Direct wave speed estimate
    S1 = min(u1-a1,u4-a4);
    S4 = max(u1+a1,u4+a4);

    % Direct wave speed estimate based on Roe
%     uw = (sqrt(rho1)*u1+sqrt(rho4)*u4)/(sqrt(rho1)+sqrt(rho4));
%     H1 = gamma/(gamma-1)*p1/rho1+0.5*u1^2;
%     H4 = gamma/(gamma-1)*p4/rho4+0.5*u4^2;
%     Hw = (sqrt(rho1)*H1+sqrt(rho4)*H4)/(sqrt(rho1)+sqrt(rho4));
%     aw = sqrt((gamma-1)*(Hw-0.5*uw^2));
%     S1 = uw-aw;
%     S4 = uw+aw;
    
    
    Lambda = max(abs([ S1 S4 Lambda]));
    
    Q1 = [ rho1; rho1*u1; 1/(gamma-1)*p1+0.5*rho1*u1^2];
    Q4 = [ rho4; rho4*u4; 1/(gamma-1)*p4+0.5*rho4*u4^2];
    
    F1 = [ rho1*u1; rho1*u1^2+p1; (gamma/(gamma-1))*u1*p1+0.5*rho1*u1^3];
    F4 = [ rho4*u4; rho4*u4^2+p4; (gamma/(gamma-1))*u4*p4+0.5*rho4*u4^3];

    if S1>=0
        F(:,i) = F1;
    elseif S4>=0
        F(:,i) = (S4*F1-S1*F4+S1*S4*(Q4-Q1))/(S4-S1);
    elseif S4<=0
        F(:,i) = F4;
    end
    
end %for i loop