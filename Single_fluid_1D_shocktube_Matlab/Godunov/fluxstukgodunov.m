% i = j - 1/2
for i  = 1:N+1

if wwl(1,i) == wwr(1,i+1) & wwl(2,i) == wwr(2,i+1) & wwl(2,i) == wwr(2,i+1)
    
    rhof = wwl(1,i);
    uf   = wwl(2,i);
    pf   = wwl(3,i);
    
    
else
    
rhol = wwl(1,i);
ul   = wwl(2,i);
pl   = wwl(3,i);

rhor = wwr(1,i+1);
ur   = wwr(2,i+1);
pr   = wwr(3,i+1);


%initial speed of sound
al   = sqrt(gamma*pl/rhol);
ar   = sqrt(gamma*pr/rhor);

ptest2 = 0.5;
ptestold2 = 0;

%SHOCK OR EXPANSION
if ul > ur
    if pl > pr
        while abs(ptest2-ptestold2)>0.0001
            ptestold2 = ptest2;            
            ptest2 = (pr + rhor*ar*sqrt(1+(gamma+1)/(2*gamma)*(ptestold2-pr)/pr)*(ul-ur))^2;
            ptest = sqrt(ptest2);
        end
        if ptest > pl
            left = 1;
            right = 1;
        elseif ptest < pl
            left = 2;
            right = 1;
        else
            left = 0;
            right = 1;
        end
    
    elseif pl < pr
        while abs(ptest2-ptestold2)>0.0001
            ptestold2 = ptest2;            
            ptest2 = (pl + rhol*al*sqrt(1+(gamma+1)/(2*gamma)*(ptestold2-pl)/pl)*(ul-ur))^2;
            ptest = sqrt(ptest2);
        end
        if ptest > pr
            left = 1;
            right = 1;
        elseif ptest < pr
            left = 1;
            right = 2;
        else
            left = 1;
            right = 0;
        end
            
    else
        left = 1;
        right = 1;
    end
    
elseif ul < ur
    if pl > pr
        while abs(real(ptest2-ptestold2))>0.0001
            ptestold2 = ptest2;
            ptest2 = (pl+rhol*al*(gamma-1)/(2*gamma)*(1-ptestold2/pl)/(1-(ptestold2/pl)^((gamma-1)/(2*gamma))+small)*(ul-ur))^2;
            ptest = sqrt(ptest2);
        end
        if ptest > pr
            left = 2;
            right = 1;
        elseif ptest < pr
            left = 2;
            right = 2;
        else
            left = 2;
            right = 0;
        end
        
    elseif pl < pr
        while abs(ptest2-ptestold2)>0.0001
            ptestold2 = ptest2;            
            ptest2 = (pr+rhor*ar*(gamma-1)/(2*gamma)*(1-ptestold2/pr)/(1-(ptestold2/pr)^((gamma-1)/(2*gamma))+small)*(ul-ur))^2;
            ptest = sqrt(ptest2);
        end
        if ptest > pl
            left = 2;
            right = 1;
        elseif ptest < pl
            left = 2;
            right = 2;
        else
            left = 2;
            right = 0;
        end
                
    else
        left = 2;
        right = 2;
    end
    
else
    if pl > pr
        left = 2;
        right = 1;
    elseif pl < pr
        left = 1;
        right = 2;
    else
        left = 0;
        right = 0;
    end
end



%initial guess final state pressure
pf01 = 10;%(((gamma-1)/2*(ul-ur)+al+ar)/(al*pl^((2*gamma)/(gamma-1))+ar*pr^((2*gamma)/(gamma-1))))^((2*gamma)/(gamma-1));
pf02 = 0.1;%(rhol*al*pr+rhor*ar*pl-rhol*al*rhor*ar*(ur-ul))/(rhol*al+rhor*ar);

k = 0;
pfnew = 0.5;
pfold = 0;

while abs(pfnew-pfold)>0.00001
k = k + 1;

if k==1
pf = pf02;
elseif k==2
pf = pf01;
else
pf = pfnew;
end

if left == 1 %shock
ml = rhol*al*sqrt(1+(gamma+1)/(2*gamma)*(pf-pl)/pl);
elseif left == 2 % expansion
ml = rhol*al*(gamma-1)/(2*gamma)*(1-pf/pl)/(1-(pf/pl)^((gamma-1)/(gamma*2))+small);
else
ml = 0;
end

if right == 1 %shock
mr = rhor*ar*sqrt(1+(gamma+1)/(2*gamma)*(pf-pr)/pr);
elseif right == 2 %expansion
mr = rhor*ar*(gamma-1)/(2*gamma)*(1-pf/pr)/(1-(pf/pr)^((gamma-1)/(gamma*2))+small);
else
mr = 0;
end

%pressure of final state
G = (ml*pr+mr*pl-ml*mr*(ur-ul))/(ml+mr+small) - pf;

if k>1
pfnew2 = (pf - (G*(pf-pfold))/(G-Gold+small))^2;
pfnew = sqrt(pfnew2);
end

%save data
Gold = G;
pfold = pf;
end
pf = real(pfnew);

%velocity of final state
if left == 1 %shock
ml = rhol*al*sqrt(1+(gamma+1)/(2*gamma)*(pf-pl)/pl);
elseif left == 2 %expansion
ml = rhol*al*(gamma-1)/(2*gamma)*(1-pf/pl)/(1-(pf/pl)^((gamma-1)/(gamma*2))+small);
else
ml = 0;
end
    
if right == 1 % shock 
mr = rhor*ar*sqrt(1+(gamma+1)/(2*gamma)*(pf-pr)/pr);
elseif right == 2 %expansion
mr = rhor*ar*(gamma-1)/(2*gamma)*(1-pf/pr)/(1-(pf/pr)^((gamma-1)/(gamma*2))+small);    
else
mr = 0;
end

if (mr+ml)>0
uf = (ml*ul+mr*ur-(pr-pl))/(ml+mr);
else
uf = ur;
mr = mr
ml = ml
pause
end

%density
if left == 1 %shock
rho2 = rhol*(1+(gamma+1)/(gamma-1)*pf/pl)/((gamma+1)/(gamma-1)+pf/pl);
elseif left == 2 %expansion
rho2 = rhol*(pf/pl)^(1/gamma);
else
rho2 = rhol;
end

if right ==1 %shock
rho3 = rhor*(1+(gamma+1)/(gamma-1)*pf/pr)/((gamma+1)/(gamma-1)+pf/pr);
elseif right == 2 %expansion
rho3 = rhor*(pf/pr)^(1/gamma);
else
rho3 = rhor;
end

%speed of sound
a2 = sqrt(gamma*pf/rho2);
a3 = sqrt(gamma*pf/rho3);


%direction characteristics
dxdtfirst1 = [];
dxdtlast1  = [];
dxdts1     = [];
dxdt0      = [];
dxdtlast4  = [];
dxdtfirst4 = [];
dxdts4     = [];

if left == 1 %shock
dxdts1 = ul - ar*sqrt(1+(gamma+1)/(2*gamma)*(pf-pl)/pl);
elseif left == 2 %expansion
dxdtfirst1 = ul - al;    
dxdtlast1  = uf - a2;
else
end

if right == 1 %shock
dxdts4 = ur + ar*sqrt(1+(gamma+1)/(2*gamma)*(pf-pr)/pr);
elseif right == 2 %expansion
dxdtfirst4 = ur + ar;
dxdtlast4  = uf + a3;
else
end

%contact discontinuity
dxdt0 = uf;


% (t,x = 0) in region:

% region 1
if (left == 0 & dxdt0 > 0) | dxdtfirst1 > 0
    pf   = pl;
    uf   = ul;
    rhof = rhol;
    
%region 2
elseif dxdtlast1 < 0 & dxdt0 > 0
    pf   = pf;
    uf   = uf;
    rhof = rho2;
    
%region 3
elseif dxdt0 < 0 & dxdtlast4 > 0
    pf   = pf;
    uf   = uf;
    rhof = rho3;
    
%region 4
elseif (right == 0 & dxdt0 < 0) | dxdtfirst4 < 0 
    pf   = pr;
    uf   = ur;
    rhof = rhor;
    
%region 1/2
elseif dxdtfirst4 < 0 & dxdtlast4 > 0
    uf   = 2/(gamma+1)*ar+(gamma-1)/(gamma+1)*ur;
    af   = -(gamma-1)/(gamma+1)*ul+2/(gamma+1)*al;
    pf   = pl*(af/al).^((2*gamma)/(gamma-1));
    rhof = rhol*(af/al).^(2/(gamma-1));
    
%region 3/4
elseif dxdtlast4 < 0 & dxdtfirst4 > 0
    uf   = 2/(gamma+1)*ar+(gamma-1)/(gamma+1)*ur;
    af   = -(gamma-1)/(gamma+1)*ur+2/(gamma+1)*ar;
    pf   = pr*(af/ar).^((2*gamma)/(gamma-1));
    rhof = rhor*(af/ar).^(2/(gamma-1));

else % no waves
    pf   = pl;
    uf   = ul;
    rhof = rhol;
end

end %end if leftstate == rightstate


f(:,i) = [rhof*uf,rhof*uf^2+pf,gamma/(gamma-1)*uf*pf+0.5*rhof*uf^3];
end