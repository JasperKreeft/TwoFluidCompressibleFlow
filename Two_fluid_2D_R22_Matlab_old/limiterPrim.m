function [W,Wfx,Wfy]=limiterPrim(dx,dy,Q,gamma1,gamma2)

%--------------------------------------------------%
%                                                  %
% Limiters, used to avoid non-monotonicity         %
%                                                  %
%                                                  %
% Primitive variant                                %
%--------------------------------------------------%

global NX NY

lim = 2;

epss = rand(1)*100*eps;

W = zeros(NY,NX,8);
Wfx  = zeros(NY,2*NX,8);
Wfy  = zeros(2*NY,NX,8);
dWLx = zeros(NY,NX,8);
dWRx = zeros(NY,NX,8);
dWAy = zeros(NY,NX,8);
dWBy = zeros(NY,NX,8);

% for i=1:NY
%     for j=1:NX
%         WRK3(i,j,1)=Q(i,j,1);
%         WRK3(i,j,2)=Q(i,j,2);
%         WRK3(i,j,3)=Q(i,j,3);
%         WRK3(i,j,4)=Q(i,j,4)/Q(i,j,3);
%         WRK3(i,j,5)=Q(i,j,5)/Q(i,j,3);
%         WRK3(i,j,6)=(gamma1-gamma2)*(Q(i,j,8)-((Q(i,j,7)*(Q(i,j,4)^(2e0)+Q(i,j,5)^(2e0)))/(2e0*Q(i,j,3)^(2e0))))+(gamma2-1e0)*(Q(i,j,6)-(Q(i,j,4)^(2e0)+Q(i,j,5)^(2e0))/(2e0*Q(i,j,3)));
%         WRK3(i,j,7)=Q(i,j,7)/Q(i,j,3);
%         WRK3(i,j,8)=((gamma1-1e0)*(Q(i,j,8)-Q(i,j,7)*(Q(i,j,4)^(2e0)+Q(i,j,5)^(2e0))/(2e0*Q(i,j,3)^(2e0))))/((gamma1-gamma2)*(Q(i,j,8)-Q(i,j,7)*(Q(i,j,4)^(2e0)+Q(i,j,5)^(2e0))/(2e0*Q(i,j,3)^(2e0)))+(gamma2-1e0)*(Q(i,j,6)-(Q(i,j,4)^(2e0)+Q(i,j,5)^(2e0))/(2e0*Q(i,j,3))));
%     end
% end

W(:,:,1)=Q(:,:,1);
W(:,:,2)=Q(:,:,2);
W(:,:,3)=Q(:,:,3);
W(:,:,4)=Q(:,:,4)./Q(:,:,3);
W(:,:,5)=Q(:,:,5)./Q(:,:,3);
W(:,:,6)=(gamma1-gamma2)*(Q(:,:,8)-((Q(:,:,7).*(Q(:,:,4).^2+Q(:,:,5).^2))./(2*Q(:,:,3).^2)))+(gamma2-1)*(Q(:,:,6)-(Q(:,:,4).^2+Q(:,:,5).^2)./(2*Q(:,:,3)));
W(:,:,7)=Q(:,:,7)./Q(:,:,3);
W(:,:,8)=((gamma1-1)*(Q(:,:,8)-Q(:,:,7).*(Q(:,:,4).^2+Q(:,:,5).^2)./(2*Q(:,:,3).^2)))./((gamma1-gamma2)*(Q(:,:,8)-Q(:,:,7).*(Q(:,:,4).^2+Q(:,:,5).^2)./(2*Q(:,:,3).^2))+(gamma2-1e0)*(Q(:,:,6)-(Q(:,:,4).^(2e0)+Q(:,:,5).^2)./(2*Q(:,:,3))));


% x-direction
% First-Order

if lim == 0
    % nothing
else
    dw =diff(W,1,2);
    dWLx(:,1,3:8)  = dw(:,1,3:8)/2;
    dWLx(:,NX,3:8) = 0;
    dWRx(:,1,3:8)  = 0;
    dWRx(:,NX,3:8) = -dw(:,NX-1,3:8)/2;
    rL = (dw(:,2:NX-1,3:8)+epss)./(dw(:,1:NX-2,3:8)+epss);
    rR = 1./rL;
    if lim==1
        % Minmod limiter
        phiL = (rL>=0).*(rL+(rL>=1).*(1-rL));
        phiR = (rR>=0).*(rR+(rR>=1).*(1-rR));
    elseif lim==2
        % Koren limiter
        phiL = (rL>=0).*( 2*rL + (1/3-4/3*rL).*(rL>=0.25) + (5/3-2/3*rL).*(rL>=2.5) );
        phiR = (rR>=0).*( 2*rR + (1/3-4/3*rR).*(rR>=0.25) + (5/3-2/3*rR).*(rR>=2.5) );
    end
    dWLx(:,2:NX-1,3:8) =  1/2*phiL.*dw(:,1:NX-2,3:8);
    dWRx(:,2:NX-1,3:8) = -1/2*phiR.*dw(:,2:NX-1,3:8);
end    
    
%y-direction
if lim==0
    % nothing
else
    dw =diff(W,1,1);
    dWAy(1,:,3:8)  = dw(1,:,3:8)/2;
    dWAy(NY,:,3:8) = 0;
    dWBy(1,:,3:8)  = 0;
    dWBy(NY,:,3:8) = -dw(NY-1,:,3:8)/2;
    rA = (dw(2:NY-1,:,3:8)+epss)./(dw(1:NY-2,:,3:8)+epss);
    rB = 1./rA;
    if lim==1
        % Minmod limiter
        phiA = (rA>=0).*(rA+(rA>=1).*(1-rA));
        phiB = (rB>=0).*(rB+(rB>=1).*(1-rB));
    elseif lim==2
        % Koren limiter
        phiA = (rA>=0).*( 2*rA + (1/3-4/3*rA).*(rA>=0.25) + (5/3-2/3*rA).*(rA>=2.5) );
        phiB = (rB>=0).*( 2*rB + (1/3-4/3*rB).*(rB>=0.25) + (5/3-2/3*rB).*(rB>=2.5) );
    end
    dWAy(2:NY-1,:,3:8) =  1/2*phiA.*dw(1:NY-2,:,3:8);
    dWBy(2:NY-1,:,3:8) = -1/2*phiB.*dw(2:NY-1,:,3:8);
end
    
    
    

% else
%     for i=2:NY-1
%         for j=2:NX-1
%             for k=1:4
%                 if j==1
%                     dWLx(i,j,k) = 1/2*(w(i,j+1,k+2)-w(i,j,k+2));
%                     dWRx(i,j,k) = 0;
%                 elseif j==NX
%                     dWLx(i,j,k) = 0;
%                     dWRx(i,j,k) = 1/2*(w(i,j-1,k+2)-w(i,j,k+2));
%                 else
%                     r = (w(i,j+1,k)-w(i,j,k)+eps)/(w(i,j,k)-w(i,j-1,k)+eps);
%                     if limiter == 1 % Minmod
%                         phiL = max(0,min(1,r));
%                         phiR = max(0,min(1,1/r));
%                     elseif limiter == 2 % Koren
%                         phiL = max(0,min(2*r,min(1/3+2/3*r,2)));
%                         phiR = max(0,min(2/r,min(1/3+2/3/r,2)));
%                     elseif limiter == 3 % Van Leer                        
%                         phiL = (r+abs(r))/(1+r);
%                         phiR = (1/r+abs(1/r))/(1+1/r);
%                     elseif limiter == 4 % Van Albada
%                         phiL = (r^2+r)/(r^2+1);
%                         phiR = ((1/r)^2+1/r)/((1/r)^2+1);
%                     elseif limiter == 5 % HQUICK                        
%                         phiL = 2*(r+abs(r))/(r+3);
%                         phiR = 2*((1/r)+abs(1/r))/((1/r)+3);
%                     end 
%                     dWLx(i,j,k) = 1/2*phiL*(w(i,j,k)-w(i,j-1,k));
%                     dWRx(i,j,k) = -1/2*r*phiR*(w(i,j,k)-w(i,j-1,k));
%                 end
% 
%                 if i==1
%                     dWAy(i,j,k) = 1/2*(w(i+1,j,k)-w(i,j,k));
%                     dWBy(i,j,k) = 0;
%                 elseif i==NY
%                     dWAy(i,j,k) = 0;
%                     dWBy(i,j,k) = 1/2*(w(i-1,j,k)-w(i,j,k));
%                 else
%                     r = (w(i+1,j,k)-w(i,j,k)+eps)/(w(i,j,k)-w(i-1,j,k)+eps);
%                     if limiter == 1 % Minmod
%                         phiA = max(0,min(1,r));
%                         phiB = max(0,min(1,1/r));
%                     elseif limiter == 2 % Koren
%                         phiA = max(0,min(2*r,min(1/3+2/3*r,2)));
%                         phiB = max(0,min(2/r,min(1/3+2/3/r,2)));
%                     elseif limiter == 3 % Van Leer
%                         phiA = (r+abs(r))/(1+r);
%                         phiB = (1/r+abs(1/r))/(1+1/r);
%                     elseif limiter == 4 % Van Albada                        
%                         phiA = (r^2+r)/(r^2+1);
%                         phiB = ((1/r)^2+1/r)/((1/r)^2+1);  
%                     elseif limiter == 5 % HQUICK
%                         phiA = 2*(r+abs(r))/(r+3);
%                         phiB = 2*((1/r)+abs(1/r))/((1/r)+3);
%                     end
%                     dWAy(i,j,k) = 1/2*phiA*(w(i,j,k)-w(i,j-1,k));
%                     dWBy(i,j,k) = -1/2*r*phiB*(w(i,j,k)-w(i,j-1,k));
%                 end
%             end
%         end
%     end
% end
    
    
    
    
    
    
    
    
    
    
    
    
% % Minmod limiter
% elseif lim == 1
%     for i=1:NY
%         for j=1:NX
%             for k=1:6
%                 if j==1
%                     dWLx(i,j,k) = 1e0/2e0*(WRK3(i,j+1,k+2)-WRK3(i,j,k+2));
%                     dWRx(i,j,k) = 0;
%                 elseif j==NX
%                     dWLx(i,j,k) = 0;
%                     dWRx(i,j,k) = 1e0/2e0*(WRK3(i,j-1,k+2)-WRK3(i,j,k+2));
%                 else
%                     rL = (WRK3(i,j+1,k+2)-WRK3(i,j,k+2)+eps)/(WRK3(i,j,k+2)-WRK3(i,j-1,k+2)+eps);
%                     rR = (WRK3(i,j-1,k+2)-WRK3(i,j,k+2)+eps)/(WRK3(i,j,k+2)-WRK3(i,j+1,k+2)+eps);
%                     if rL<0e0
%                       phiL = 0e0;
%                     elseif rL<1e0
%                       phiL = rL;
%                     else
%                       phiL = 1e0;
%                     end
%                    if rR<0e0
%                      phiR = 0e0;
%                    elseif rR<1e0
%                      phiR = rR;
%                    else
%                      phiR = 1e0;
%                    end
%                     dWLx(i,j,k) = 1e0/2e0*phiL*(WRK3(i,j,k+2)-WRK3(i,j-1,k+2));
%                     dWRx(i,j,k) = 1e0/2e0*phiR*(WRK3(i,j,k+2)-WRK3(i,j+1,k+2));
%                 end
%             end
%         end
%     end
% 
% 
% % Koren limiter
% elseif lim == 2
%     for i=1:NY
%         for j=1:NX
%             for k=1:6
%                 if j==1
%                     dWLx(i,j,k) = 1e0/2e0*(WRK3(i,j+1,k+2)-WRK3(i,j,k+2));
%                     dWRx(i,j,k) = 0;
%                 elseif j==NX
%                     dWLx(i,j,k) = 0;
%                     dWRx(i,j,k) = 1e0/2e0*(WRK3(i,j-1,k+2)-WRK3(i,j,k+2));
%                 else
%                     rL = (WRK3(i,j+1,k+2)-WRK3(i,j,k+2)+eps)/(WRK3(i,j,k+2)-WRK3(i,j-1,k+2)+eps);
%                     rR = (WRK3(i,j-1,k+2)-WRK3(i,j,k+2)+eps)/(WRK3(i,j,k+2)-WRK3(i,j+1,k+2)+eps);
%                     if rL<0e0
%                       phiL = 0e0;
%                     elseif rL<0.25e0
%                       phiL = 2e0*rL;
%                     elseif rL<2.5e0
%                       phiL = 1e0/3e0+2e0/3e0*rL;
%                     else
%                       phiL = 2e0;
%                     end
%                     if rR<0e0
%                       phiR = 0e0;
%                     elseif rR<0.25e0
%                       phiR = 2e0*rR;
%                     elseif rR<2.5e0
%                       phiR = 1e0/3e0+2e0/3e0*rR;
%                     else
%                       phiR = 2e0;
%                     end
%                     dWLx(i,j,k) = 1e0/2e0*phiL*(WRK3(i,j,k+2)-WRK3(i,j-1,k+2));
%                     dWRx(i,j,k) = 1e0/2e0*phiR*(WRK3(i,j,k+2)-WRK3(i,j+1,k+2));
%                 end
%             end
%         end
%     end
% end
% 
% 
% %y-direction
% if lim == 0
%     for i=1:NY
%         for j=1:NX
%             dWBy(i,j,3:8)=0;
%             dWAy(i,j,3:8)=0;
%         end
%     end
%     
% elseif lim == 1
%     for i=1:NY
%         for j=1:NX
%             for k=1:6
%                 if i==1
%                     dWAy(i,j,k) = 1e0/2e0*(WRK3(i+1,j,k+2)-WRK3(i,j,k+2));
%                     dWBy(i,j,k) = 0;
%                 elseif i==NY
%                     dWAy(i,j,k) = 0;
%                     dWBy(i,j,k) = 1e0/2e0*(WRK3(i-1,j,k+2)-WRK3(i,j,k+2));
%                 else
%                     rA = (WRK3(i+1,j,k+2)-WRK3(i,j,k+2)+eps)/(WRK3(i,j,k+2)-WRK3(i-1,j,k+2)+eps);
%                     rB = (WRK3(i-1,j,k+2)-WRK3(i,j,k+2)+eps)/(WRK3(i,j,k+2)-WRK3(i+1,j,k+2)+eps);
%                     if rA<0e0
%                       phiA = 0e0;
%                     elseif rA<1e0
%                       phiA = rA;
%                     else
%                       phiA = 1e0;
%                     end
%                    if rB<0e0
%                      phiB = 0e0;
%                    elseif rB<1e0
%                      phiB = rB;
%                    else
%                      phiB = 1e0;
%                    end
%                     dWAy(i,j,k) = 1e0/2e0*phiA*(WRK3(i,j,k+2)-WRK3(i-1,j,k+2));
%                     dWBy(i,j,k) = 1e0/2e0*phiB*(WRK3(i,j,k+2)-WRK3(i+1,j,k+2));
% 
%                 end
%             end
%         end
%     end
% 
% 
% elseif lim == 2
%     for i=1:NY
%         for j=1:NX
%             for k=1:6
%                 if i==1
%                     dWAy(i,j,k) = 1e0/2e0*(WRK3(i+1,j,k+2)-WRK3(i,j,k+2));
%                     dWBy(i,j,k) = 0;
%                 elseif i==NY
%                     dWAy(i,j,k) = 0;
%                     dWBy(i,j,k) = 1e0/2e0*(WRK3(i-1,j,k+2)-WRK3(i,j,k+2));
%                 else
%                     rA = (WRK3(i+1,j,k+2)-WRK3(i,j,k+2)+eps)/(WRK3(i,j,k+2)-WRK3(i-1,j,k+2)+eps);
%                     rB = (WRK3(i-1,j,k+2)-WRK3(i,j,k+2)+eps)/(WRK3(i,j,k+2)-WRK3(i+1,j,k+2)+eps);
%                     if rA<0e0
%                       phiA = 0e0;
%                     elseif rA<0.25e0
%                       phiA = 2e0*rA;
%                     elseif rA<2.5e0
%                       phiA = 1e0/3e0+2e0/3e0*rA;
%                     else
%                       phiA = 2e0;
%                     end
%                     if rB<0e0
%                       phiB = 0e0;
%                     elseif rB<0.25e0
%                       phiB = 2e0*rB;
%                     elseif rR<2.5e0
%                       phiB = 1e0/3e0+2e0/3e0*rB;
%                     else
%                       phiB = 2e0;
%                     end
%                     dWAy(i,j,k) = 1e0/2e0*phiA*(WRK3(i,j,k+2)-WRK3(i-1,j,k+2));
%                     dWBy(i,j,k) = 1e0/2e0*phiB*(WRK3(i,j,k+2)-WRK3(i+1,j,k+2));
%                 end
%             end
%         end
%     end
% end

for i=1:NY
    for j=1:NX
        Wfx(i,2*j-1,1)   = W(i,j,1)-0.5*dx;
        Wfx(i,2*j-1,2)   = W(i,j,2);
        Wfx(i,2*j,1)     = W(i,j,1)+0.5*dx;
        Wfx(i,2*j,1)     = W(i,j,2);
        Wfy(2*i-1,j,1)   = W(i,j,1);
        Wfy(2*i-1,j,2)   = W(i,j,2)-0.5*dy;
        Wfy(2*i,j,1)     = W(i,j,1);
        Wfy(2*i,j,2)     = W(i,j,2)+0.5*dy;
        for k = 1:6
            Wfx(i,2*j-1,k+2) = W(i,j,k+2) + dWRx(i,j,k);
            Wfx(i,2*j,k+2)   = W(i,j,k+2) + dWLx(i,j,k);
            Wfy(2*i-1,j,k+2) = W(i,j,k+2) + dWBy(i,j,k);
            Wfy(2*i,j,k+2)   = W(i,j,k+2) + dWAy(i,j,k);
        end
    end
end