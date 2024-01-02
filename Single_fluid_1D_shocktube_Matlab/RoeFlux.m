% Roe flux solver

function [F,Lambda] = RoeFlux(N,Wf,gamma,Lambda)

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

    rhow = sqrt(rho1*rho4);
    uw = (sqrt(rho1)*u1+sqrt(rho4)*u4)/(sqrt(rho1)+sqrt(rho4));
    H1 = gamma/(gamma-1)*p1/rho1+0.5*u1^2;
    H4 = gamma/(gamma-1)*p4/rho4+0.5*u4^2;
    Hw = (sqrt(rho1)*H1+sqrt(rho4)*H4)/(sqrt(rho1)+sqrt(rho4));
    aw = sqrt((gamma-1)*(Hw-0.5*uw^2));

    K = [ 1        1        1       ;
          uw-aw    uw       uw+aw   ;
          Hw-uw*aw 0.5*uw^2 Hw+uw*aw];

    alpha = [ ((p4-p1)-rhow*aw*(u4-u1))/(2*aw^2);
              (rho4-rho1)-(p4-p1)/aw^2          ;
              ((p4-p1)+rhow*aw*(u4-u1))/(2*aw^2)];

    % without entropy fix
%     lambda = [ uw-aw ; uw ; uw+aw ];
      
    % with entropy fix    
    % Roe-averaged states estimator
%     rhostarL = rho1+alpha(1);
%     ustar = (rho1*u1+alpha(1)*(uw-aw))/(rho1+alpha(1));
%     pstar = (gamma-1)*(p1/(rho1*(gamma-1))+0.5*u1^2+alpha(1)*(Hw-uw*aw)-0.5*rhostarL*ustar^2);
%     astarL = sqrt(gamma*pstar/rhostarL);
%     rhostarR = rho4-alpha(3);
% %     ustarR = (rho4*u4-alpha(3)*(uw+aw))/(rho4-alpha(3));
% %     pstarR = (gamma-1)*(p4/(rho4*(gamma-1))+0.5*u4^2-alpha(3)*(Hw+uw*aw)-0.5*rhostarR*ustarR^2);
%     astarR = sqrt(gamma*pstar/rhostarR);
      
    % Two-rarefaction estimator
    z = (gamma-1)/(2*gamma);
    pstar  = ((a1+a4-(gamma-1)/2*(u4-u1))/(a1/(p1^z)+a4/(p4^z)))^(1/z);
    astarL = a1*(pstar/p1)^z;
    ustar  = u1+2/(gamma-1)*(a1-astarL);
    astarR = a4*(pstar/p4)^z;
%     ustarR = u4+2/(gamma-1)*(astarR-a4);

          
    lambda1L = u1-a1;
    lambda1R = ustar-astarL;
    lambda3R = u4+a4;
    lambda3L = ustar+astarR;

    lambda = zeros(3,1);
    lambda(2) = uw;
    if lambda1L<0 && lambda1R>0
        lambda(1) = lambda1L*(lambda1R-(uw-aw))/(lambda1R-lambda1L);
    else
        lambda(1) = uw-aw;
    end
    if lambda3L<0 && lambda3R>0
        lambda(3) = lambda3L*((uw+aw)-lambda3L)/(lambda3R-lambda3L);
    else
        lambda(3) = uw+aw;
    end
         

    F1 = [ rho1*u1; rho1*u1^2+p1; (gamma/(gamma-1))*u1*p1+0.5*rho1*u1^3];
    F4 = [ rho4*u4; rho4*u4^2+p4; (gamma/(gamma-1))*u4*p4+0.5*rho4*u4^3];

    waveflux = zeros(3,1);
    for j=1:3
        waveflux = waveflux+alpha(j)*abs(lambda(j))*K(:,j);
    end
    F(:,i) = (F1+F4)/2-0.5*waveflux;
    
    Lambda = max(abs([ lambda; Lambda]));

end %for i loop