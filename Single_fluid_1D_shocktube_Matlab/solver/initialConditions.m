% Initial conditions

n = 1;

if nr~=10
test     = zeros(9,1);
test(nr) = 1;

gamma            = [ 1.4    1.4  1.4   1.4    1.4    1.4    1.4    1.4     1.4]*test;    % ratio of specific heat

W(1,1:(N/2))   = [ 1.0 1000.0  1.0  10.0    1.0    1.0    1.0    5.99924 1  ]*test;    % Density Left state
W(2,1:(N/2))   = [ 1.0    1.0  0.0   0.0   -2.0    0.0    0.0   19.5975  0  ]*test;    % Velocity Left state
W(3,1:(N/2))   = [ 1.0    1.0  1.0  10.0    0.4 1000.0    0.01 460.894   1  ]*test;    % Pressure Left state

W(1,(N/2+1):N) = [ 0.1    1.0  0.125 0.125  1.0    1.0    1.0    5.99242 1  ]*test; % Density Right state
W(2,(N/2+1):N) = [ 1.0    1.0  0.0   0.0    2.0    0.0    0.0   -6.19633 0  ]*test; % Velocity Right state
W(3,(N/2+1):N) = [ 1.0    1.0  0.1   0.1    0.4    0.01 100.0   46.0950  1  ]*test; % Pressure Right state

else
% gamma = 1.4;
% W(1,1:N,1) = 1.0;
% W(2,1:N,1) = 0.0;
% W(3,1:round(N/10),1)                 = 1000.0 ;
% W(3,(round(N/10)+1):round(9*N/10),1) =    0.01;
% W(3,(round(9*N/10)+1):N,1)           =  100.0 ;

gamma = 1.4;
W(1,1:(N*6/10)) = 2171;
W(2,1:(N*6/10)) = 0.0;
W(3,1:(N*6/10)) = 2e6 ;
W(1,(N*6/10+1):N) = 2171;
W(2,(N*6/10+1):N) = 0.0;
W(3,(N*6/10+1):N) = 1 ;
end