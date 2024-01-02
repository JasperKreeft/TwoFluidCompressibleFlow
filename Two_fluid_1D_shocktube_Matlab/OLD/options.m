%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              %
% Written by: Jasper Kreeft  (2007)            %
% Updated by: Jasper Kreeft  (2015)            %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem
disp('Choose initial conditions:')
disp('  ')
disp('     case  1: Translating interface problem')
disp('     case  2: Translating interface problem with strong density jump')
disp('     case  3: Two-fluid Sod problem')
disp('     case  4: Two-fluid high pressure high density Sod problem')
disp('     case  5: No-reflection problem')
disp('     case  6: Shock Helium-bubble interaction')
disp('     case  7: Shock R22-bubble interaction')
disp('     case  8: Mixture problem')
disp('  ')
disp('Stiffened gas problems:')
disp('     case  9: Water-air problem (Murrone 5.2.1)')
disp('     case 10: Epoxy-Spinel')
disp('  ')
nr = [];
while isempty(nr) || nr<1 || nr>10 || rem(nr,1)~=0
nr = input('enter the number of the case you want to solve: ');
end
disp('  ')

%% Grid
N   = [];
while isempty(N)
N   = input('enter the number of grid cells: ');
end
disp('  ')

%% simulation time
T = [];
while isempty(T) || T<0
T   = input('enter the simulation time (sec): ');
end
disp('  ')


%% Temporal method
time = [];
disp('Choose time discretization method:')
disp('  ')
disp('  1. Forward Euler / ERK1    (first-order accurate)')
disp('  2. Explicit Runge-Kutta 2b (second-order accurate)')
disp('  3. Explicit Runge-Kutta 3b (third-order accurate)')
disp('  ')
while isempty(time) || time<1 || time>3 || rem(time,1)~=0
time = input('enter the temperal method you want to use:  ');
end

%% Flux solver
flux = [];
disp('Choose type of flux solver:')
disp('  ')
disp('  1. exact/Godunov')
disp('  2. Osher')
disp('  3. linearized Osher')
while isempty(flux) || flux<1 || flux>3 || rem(flux,1)~=0
flux = input('enter the type of flux solver you want to use:  ');
end

%% Boundary conditions
BC_L = [];     % Left boundary
BC_R = [];     % right boundary

disp('Choose type of boundary condtion:')
disp('  ')
disp('  1. open end')
disp('  2. closed end')
disp('  ')
while isempty(BC_L) || BC_L<1 || BC_L>2 || rem(BC_L,1)~=0
BC_L = input('enter the type of boundary condition for the left boundary: ');
end
disp('  ')
while isempty(BC_R) || BC_R<1 || BC_R>2 || rem(BC_R,1)~=0
BC_R = input('enter the type of boundary condition for the right boundary: ');
end

%% Limiter
limiter = [];
disp('Choose type of limiter:')
disp('  ')
disp('   0 = no limiter')
disp('   1 = minmod limiter')
disp('   2 = Koren limiter')
disp('   3 = Van Leer limiter')
disp('   4 = Van Albada limiter')
disp('   5 = Superbee limiter')
disp('  ')
while isempty(limiter) || limiter<0 || limiter>5 ||rem(limiter,1)~=0
limiter = input('enter the number of the limiter you want to use: ');
end

%% animation
ani = [];
disp('Choose type of animation during calculation')
disp('  ')
disp('   0 = no animation')
disp('   1 = waitbar')
disp('   2 = density plot')
disp('   3 = rho,u,p,e plots')
disp('   4 = density in tube')
disp('  ')
while isempty(ani) || ani<0 || ani>4 || rem(ani,1)~=0
ani = input('enter the number of animation you want to see:  ');
end

%% Create movie
movie = [];
while isempty(movie) || ani<0 || ani>4 || rem(ani,1)~=0
movie = input('Create a movie (1=Yes, 2=No)? ');
end