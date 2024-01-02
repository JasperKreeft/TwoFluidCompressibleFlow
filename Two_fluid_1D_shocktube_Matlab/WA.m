gamma1         = 4.4; % ratio of specific heat left state
gamma2         = 1.4; % ratio of specific heat right state

W(1,1:(N/2))   = 1000; % Density left state
W(2,1:(N/2))   = ; % Velocity left state
W(3,1:(N/2))   = ; % pressure left state
W(4,1:(N/2))   = ; % mass fraction left state
W(5,1:(N/2))   = 1.0; % volume fraction left state

W(1,(N/2+1):N) = ; % Density right state
W(2,(N/2+1):N) = ; % Velocity right state
W(3,(N/2+1):N) = ; % Pressure right state
W(4,(N/2+1):N) = ; % Mass fraction right state
W(5,(N/2+1):N) = ; % Volume fraction right state