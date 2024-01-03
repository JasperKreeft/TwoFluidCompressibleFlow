function [WRK3,Wfx,Wfy]=limiterPrim(NX,NY,dx,dy,QRK3,gamma1,gamma2,eps)

%--------------------------------------------------%
%                                                  %
% Limiters, used to avoid non-monotonicity         %
%                                                  %
%                                                  %
% Primitive variant                                %
%--------------------------------------------------%

lim = 2;

WRK3 = zeros(NY,NX,8);
Wfx  = zeros(NY,2*NX,8);
Wfy  = zeros(2*NY,NX,8);
dWLx = zeros(NY,NX,8);
dWRx = zeros(NY,NX,8);
dWAy = zeros(NY,NX,8);
dWBy = zeros(NY,NX,8);

for i=1:NY
    for j=1:NX
        WRK3(i,j,1)=QRK3(i,j,1);
        WRK3(i,j,2)=QRK3(i,j,2);
        WRK3(i,j,3)=QRK3(i,j,3);
        WRK3(i,j,4)=QRK3(i,j,4)/QRK3(i,j,3);
        WRK3(i,j,5)=QRK3(i,j,5)/QRK3(i,j,3);
        WRK3(i,j,6)=(gamma1-gamma2)*(QRK3(i,j,8)-((QRK3(i,j,7)*(QRK3(i,j,4)^(2e0)+QRK3(i,j,5)^(2e0)))/(2e0*QRK3(i,j,3)^(2e0))))+(gamma2-1e0)*(QRK3(i,j,6)-(QRK3(i,j,4)^(2e0)+QRK3(i,j,5)^(2e0))/(2e0*QRK3(i,j,3)));
        WRK3(i,j,7)=QRK3(i,j,7)/QRK3(i,j,3);
        WRK3(i,j,8)=((gamma1-1e0)*(QRK3(i,j,8)-QRK3(i,j,7)*(QRK3(i,j,4)^(2e0)+QRK3(i,j,5)^(2e0))/(2e0*QRK3(i,j,3)^(2e0))))/((gamma1-gamma2)*(QRK3(i,j,8)-QRK3(i,j,7)*(QRK3(i,j,4)^(2e0)+QRK3(i,j,5)^(2e0))/(2e0*QRK3(i,j,3)^(2e0)))+(gamma2-1e0)*(QRK3(i,j,6)-(QRK3(i,j,4)^(2e0)+QRK3(i,j,5)^(2e0))/(2e0*QRK3(i,j,3))));
    end
end




% x-direction
% First-Order

if lim == 0
    for i=1:NY
        for j=1:NX
            dWRx(i,j,3:8)=0;
            dWLx(i,j,3:8)=0;
        end
    end

% Minmod limiter
elseif lim == 1
    for i=1:NY
        for j=1:NX
            for k=1:6
                if j==1
                    dWLx(i,j,k) = 1e0/2e0*(WRK3(i,j+1,k+2)-WRK3(i,j,k+2));
                    dWRx(i,j,k) = 0;
                elseif j==NX
                    dWLx(i,j,k) = 0;
                    dWRx(i,j,k) = 1e0/2e0*(WRK3(i,j-1,k+2)-WRK3(i,j,k+2));
                else
                    rL = (WRK3(i,j+1,k+2)-WRK3(i,j,k+2)+eps)/(WRK3(i,j,k+2)-WRK3(i,j-1,k+2)+eps);
                    rR = (WRK3(i,j-1,k+2)-WRK3(i,j,k+2)+eps)/(WRK3(i,j,k+2)-WRK3(i,j+1,k+2)+eps);
                    if rL<0e0
                      phiL = 0e0;
                    elseif rL<1e0
                      phiL = rL;
                    else
                      phiL = 1e0;
                    end
                   if rR<0e0
                     phiR = 0e0;
                   elseif rR<1e0
                     phiR = rR;
                   else
                     phiR = 1e0;
                   end
                    dWLx(i,j,k) = 1e0/2e0*phiL*(WRK3(i,j,k+2)-WRK3(i,j-1,k+2));
                    dWRx(i,j,k) = 1e0/2e0*phiR*(WRK3(i,j,k+2)-WRK3(i,j+1,k+2));
                end
            end
        end
    end


% Koren limiter
elseif lim == 2
    for i=1:NY
        for j=1:NX
            for k=1:6
                if j==1
                    dWLx(i,j,k) = 1e0/2e0*(WRK3(i,j+1,k+2)-WRK3(i,j,k+2));
                    dWRx(i,j,k) = 0;
                elseif j==NX
                    dWLx(i,j,k) = 0;
                    dWRx(i,j,k) = 1e0/2e0*(WRK3(i,j-1,k+2)-WRK3(i,j,k+2));
                else
                    rL = (WRK3(i,j+1,k+2)-WRK3(i,j,k+2)+eps)/(WRK3(i,j,k+2)-WRK3(i,j-1,k+2)+eps);
                    rR = (WRK3(i,j-1,k+2)-WRK3(i,j,k+2)+eps)/(WRK3(i,j,k+2)-WRK3(i,j+1,k+2)+eps);
                    if rL<0e0
                      phiL = 0e0;
                    elseif rL<0.25e0
                      phiL = 2e0*rL;
                    elseif rL<2.5e0
                      phiL = 1e0/3e0+2e0/3e0*rL;
                    else
                      phiL = 2e0;
                    end
                    if rR<0e0
                      phiR = 0e0;
                    elseif rR<0.25e0
                      phiR = 2e0*rR;
                    elseif rR<2.5e0
                      phiR = 1e0/3e0+2e0/3e0*rR;
                    else
                      phiR = 2e0;
                    end
                    dWLx(i,j,k) = 1e0/2e0*phiL*(WRK3(i,j,k+2)-WRK3(i,j-1,k+2));
                    dWRx(i,j,k) = 1e0/2e0*phiR*(WRK3(i,j,k+2)-WRK3(i,j+1,k+2));
                end
            end
        end
    end
end


%y-direction
if lim == 0
    for i=1:NY
        for j=1:NX
            dWBy(i,j,3:8)=0;
            dWAy(i,j,3:8)=0;
        end
    end
    
elseif lim == 1
    for i=1:NY
        for j=1:NX
            for k=1:6
                if i==1
                    dWAy(i,j,k) = 1e0/2e0*(WRK3(i+1,j,k+2)-WRK3(i,j,k+2));
                    dWBy(i,j,k) = 0;
                elseif i==NY
                    dWAy(i,j,k) = 0;
                    dWBy(i,j,k) = 1e0/2e0*(WRK3(i-1,j,k+2)-WRK3(i,j,k+2));
                else
                    rA = (WRK3(i+1,j,k+2)-WRK3(i,j,k+2)+eps)/(WRK3(i,j,k+2)-WRK3(i-1,j,k+2)+eps);
                    rB = (WRK3(i-1,j,k+2)-WRK3(i,j,k+2)+eps)/(WRK3(i,j,k+2)-WRK3(i+1,j,k+2)+eps);
                    if rA<0e0
                      phiA = 0e0;
                    elseif rA<1e0
                      phiA = rA;
                    else
                      phiA = 1e0;
                    end
                   if rB<0e0
                     phiB = 0e0;
                   elseif rB<1e0
                     phiB = rB;
                   else
                     phiB = 1e0;
                   end
                    dWAy(i,j,k) = 1e0/2e0*phiA*(WRK3(i,j,k+2)-WRK3(i-1,j,k+2));
                    dWBy(i,j,k) = 1e0/2e0*phiB*(WRK3(i,j,k+2)-WRK3(i+1,j,k+2));

                end
            end
        end
    end


elseif lim == 2
    for i=1:NY
        for j=1:NX
            for k=1:6
                if i==1
                    dWAy(i,j,k) = 1e0/2e0*(WRK3(i+1,j,k+2)-WRK3(i,j,k+2));
                    dWBy(i,j,k) = 0;
                elseif i==NY
                    dWAy(i,j,k) = 0;
                    dWBy(i,j,k) = 1e0/2e0*(WRK3(i-1,j,k+2)-WRK3(i,j,k+2));
                else
                    rA = (WRK3(i+1,j,k+2)-WRK3(i,j,k+2)+eps)/(WRK3(i,j,k+2)-WRK3(i-1,j,k+2)+eps);
                    rB = (WRK3(i-1,j,k+2)-WRK3(i,j,k+2)+eps)/(WRK3(i,j,k+2)-WRK3(i+1,j,k+2)+eps);
                    if rA<0e0
                      phiA = 0e0;
                    elseif rA<0.25e0
                      phiA = 2e0*rA;
                    elseif rA<2.5e0
                      phiA = 1e0/3e0+2e0/3e0*rA;
                    else
                      phiA = 2e0;
                    end
                    if rB<0e0
                      phiB = 0e0;
                    elseif rB<0.25e0
                      phiB = 2e0*rB;
                    elseif rR<2.5e0
                      phiB = 1e0/3e0+2e0/3e0*rB;
                    else
                      phiB = 2e0;
                    end
                    dWAy(i,j,k) = 1e0/2e0*phiA*(WRK3(i,j,k+2)-WRK3(i-1,j,k+2));
                    dWBy(i,j,k) = 1e0/2e0*phiB*(WRK3(i,j,k+2)-WRK3(i+1,j,k+2));
                end
            end
        end
    end
end

for i=1:NY
    for j=1:NX
        Wfx(i,2*j-1,1)   = WRK3(i,j,1)-0.5*dx;
        Wfx(i,2*j-1,2)   = WRK3(i,j,2);
        Wfx(i,2*j,1)     = WRK3(i,j,1)+0.5*dx;
        Wfx(i,2*j,1)     = WRK3(i,j,2);
        Wfy(2*i-1,j,1)   = WRK3(i,j,1);
        Wfy(2*i-1,j,2)   = WRK3(i,j,2)-0.5*dy;
        Wfy(2*i,j,1)     = WRK3(i,j,1);
        Wfy(2*i,j,2)     = WRK3(i,j,2)+0.5*dy;
        for k = 1:6
            Wfx(i,2*j-1,k+2) = WRK3(i,j,k+2) + dWRx(i,j,k);
            Wfx(i,2*j,k+2)   = WRK3(i,j,k+2) + dWLx(i,j,k);
            Wfy(2*i-1,j,k+2) = WRK3(i,j,k+2) + dWBy(i,j,k);
            Wfy(2*i,j,k+2)   = WRK3(i,j,k+2) + dWAy(i,j,k);
        end
    end
end