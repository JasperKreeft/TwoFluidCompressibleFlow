%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              %
% Written by: Jasper Kreeft  (2007)            %
% Updated by: Jasper Kreeft  (2015)            %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial conditions

global gamma1 gamma2
global pi1 pi2

RK_coeff{1} = 1;

RK_coeff{2} = [   1   0
                1/2 1/2 ];
            
RK_coeff{3} = [   1   0   0
                1/4 1/4   0
                1/6 1/6 2/3 ];

n = 1;

if nr<=8
    pi1 = 0; pi2 = 0;
end

if nr<=5
    test     = zeros(5,1);
    test(nr) = 1;

    gamma1         = [1.4    1.4 1.6  1.4   1.667 ]*test; % ratio of specific heat left state
    gamma2         = [1.6    1.6 1.4  1.6   1.200 ]*test; % ratio of specific heat right state

    W(1,1:(N/2))   = [1.0 1000.0 1.0 10.0   3.1748]*test; % Density left state
    W(2,1:(N/2))   = [1.0    1.0 0.0  0.0   9.4350]*test; % Velocity left state
    W(3,1:(N/2))   = [1.0    1.0 1.0 10.0 100.0   ]*test; % pressure left state
    W(4,1:(N/2))   = [1.0    1.0 1.0  1.0   1.0   ]*test; % mass fraction left state
    W(5,1:(N/2))   = [1.0    1.0 1.0  1.0   1.0   ]*test; % volume fraction left state

    W(1,(N/2+1):N) = [0.125 1.0 0.125 0.125 1.0   ]*test; % Density right state
    W(2,(N/2+1):N) = [1.0   1.0 0.0   0.0   0.0   ]*test; % Velocity right state
    W(3,(N/2+1):N) = [1.0   1.0 0.1   0.1   1.0   ]*test; % Pressure right state
    W(4,(N/2+1):N) = [0.0   0.0 0.0   0.0   0.0   ]*test; % Mass fraction right state
    W(5,(N/2+1):N) = [0.0   0.0 0.0   0.0   0.0   ]*test; % Volume fraction right state

elseif nr==6

    %Air1
    gamma1   = 1.4;
    W(1,1:(N/2)) = 1.92691;
    W(2,1:(N/2)) = 0.33361;
    W(3,1:(N/2)) = 1.56980;
    W(4,1:(N/2)) = 1.0;
    W(5,1:(N/2)) = 1.0;

    %Air2
    W(1,(N/2-7):(N/2)) = 1.4;
    W(2,(N/2-7):(N/2)) = 0.0;
    W(3,(N/2-7):(N/2)) = 1.0;
    W(4,(N/2-7):(N/2)) = 1.0;
    W(5,(N/2-7):(N/2)) = 1.0;

    %He
    gamma2   = 1.648;
    W(1,(N/2+1):N) = 0.25463;
    W(2,(N/2+1):N) = 0.0;
    W(3,(N/2+1):N) = 1.0;
    W(4,(N/2+1):N) = 0.0;
    W(5,(N/2+1):N) = 0.0;

elseif nr==7

    %Air1
    gamma1   = 1.4;
    W(1,1:(N/2)) = 1.92691;
    W(2,1:(N/2)) = 0.33361;
    W(3,1:(N/2)) = 1.56980;
    W(4,1:(N/2)) = 1.0;
    W(5,1:(N/2)) = 1.0;

    %Air2
    W(1,(N/2-7):(N/2)) = 1.4;
    W(2,(N/2-7):(N/2)) = 0.0;
    W(3,(N/2-7):(N/2)) = 1.0;
    W(4,(N/2-7):(N/2)) = 1.0;
    W(5,(N/2-7):(N/2)) = 1.0;

    %R22
    gamma2   = 1.249;
    W(1,(N/2+1):N) = 4.41540;
    W(2,(N/2+1):N) = 0.0;
    W(3,(N/2+1):N) = 1.0;
    W(4,(N/2+1):N) = 0.0;
    W(5,(N/2+1):N) = 0.0;
    
elseif nr==8
    gamma1 = 1.4;
    W(1,1:(N/2)) = 2171;
    W(2,1:(N/2)) = 0.0;
    W(3,1:(N/2)) = 20;
    W(4,1:(N/2)) = 0.3250;
    W(5,1:(N/2)) = 0.5954;
    
    gamma2 = 1.62;
    W(1,(N/2+1):N) = 2171;
    W(2,(N/2+1):N) = 0.0;
    W(3,(N/2+1):N) = 1.0;
    W(4,(N/2+1):N) = 0.3250;
    W(5,(N/2+1):N) = 0.5954;

elseif nr==9
    % Water-air problem (Murrone 5.2.1)
    gamma1 = 1.4;
    pi1    = 0.0;
    W(1,1:(N*5/10)) = 525;
    W(2,1:(N*5/10)) = 0.0;
    W(3,1:(N*5/10)) = 1e4;
    W(4,1:(N*5/10)) = 0.5*50/525;
    W(5,1:(N*5/10)) = 0.5;

    gamma2 = 4.4;
    pi2    = 6.e3;
    W(1,(N*5/10+1):N) = 525;
    W(2,(N*5/10+1):N) = 0.0;
    W(3,(N*5/10+1):N) = 1;
    W(4,(N*5/10+1):N) = 0.5*1000/525;
    W(5,(N*5/10+1):N) = 0.5;
    
elseif nr==10
    % Epoxy-Spinel
    gamma1 = 2.43;
    pi1    = 5.3e4;
    W(1,1:(N*6/10)) = 2171;
    W(2,1:(N*6/10)) = 0.0;
    W(3,1:(N*6/10)) = 2e6;
    W(4,1:(N*6/10)) = 0.3250;
    W(5,1:(N*6/10)) = 0.5954;

    gamma2 = 1.62;
    pi2    = 141.45e4;
    W(1,(N*6/10+1):N) = 2171;
    W(2,(N*6/10+1):N) = 0.0;
    W(3,(N*6/10+1):N) = 1;
    W(4,(N*6/10+1):N) = 0.3250;
    W(5,(N*6/10+1):N) = 0.5954;

end

t(n)   = 0;
Lambda = 20;

Qf = zeros(5,2*N);


if limiter==0
    CFL = 0.95;
else
    CFL = 0.48;
end


if ani==1
    han = waitbar(0,'Please wait...','Position',[500 500 288 60]);
elseif ani==3 || ani==4
%     figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6])
%     figure(winFull{:})
    figure(winHD{:})
end

if movie
    filename = ['Case' num2str(nr) '_N' num2str(N) '_T' num2str(round(T*1000)) 'ms_ERK' num2str(time) ...
                '_Flux' num2str(flux) '_Lim' num2str(limiter)];

    writerObj = VideoWriter(filename,'MPEG-4');
    writerObj.FrameRate = 15;
    writerObj.Quality = 75;

    open(writerObj);

end