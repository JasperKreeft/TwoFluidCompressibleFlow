% options

%% Problem
disp('Choose initial conditions:')
disp('  ')
disp('     case  1: translating contact discontinuity')
disp('     case  2: translating contact discontinuity with strong density jump')
disp('     case  3: Sod problem')
disp('     case  4: high pressure high density Sod problem')
disp('     case  5: double expansion / 123 problem')
disp('     case  6: left half of blast wave problem of Woodward and Colella')
disp('     case  7: right half of blast wave problem of Woodward and Colella')
disp('     case  8: collision of shocks from cases 6 & 7')
disp('     case  9: empty')
disp('     case 10: Woodward and Colella problem')
disp('  ')
nr = 0;
while nr<1 || nr>10 || rem(nr,1)~=0
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
time = 0;
disp('Choose time discretization method:')
disp('  ')
disp('  1. Forward Euler / ERK1    (first-order accurate)')
disp('  2. Explicit Runge-Kutta 3b (third-order accurate)')
disp('  ')
while time<1 || time>2 || rem(time,1)~=0
time = input('enter the temperal method you want to use:  ');
end

%% Flux solver
flux = 0;
disp('Choose type of flux solver:')
disp('  ')
disp('  1. exact/Godunov')
disp('  2. linearized Osher')
disp('  3. Osher')
disp('  4. HLL')
disp('  5. HLLC')
disp('  6. Roe')
disp('  7. PVRS')
disp('  8. Osher with general EOS (in progress)')
while flux<1 || flux>8 || rem(flux,1)~=0
flux = input('enter the type of flux solver you want to use:  ');
end

%% Boundary conditions
BC_L = 0;     % Left boundary
BC_R = 0;     % right boundary

disp('Choose type of boundary condtion:')
disp('  ')
disp('  1. open end')
disp('  2. closed end')
disp('  ')
while BC_L<1 || BC_L>2 || rem(BC_L,1)~=0
BC_L = input('enter the type of boundary condition for the left boundary: ');
end
disp('  ')
while BC_R<1 || BC_R>2 || rem(BC_R,1)~=0
BC_R = input('enter the type of boundary condition for the right boundary: ');
end

%% Variables to limit
limitedvariables = -1;
disp('Choose kind of variables to limit')
disp('  ')
disp('   1 = primitive variables')
disp('   2 = conservative variables')
disp('   3 = Riemann invariant variables')
while limitedvariables<0 || limitedvariables>3 ||rem(limitedvariables,1)~=0
limitedvariables = input('enter the number of the variables you want to limit: ');
end

%% Limiter
if time==1
    limiter = 0;
else
limiter = -1;
disp('Choose type of limiter:')
disp('  ')
disp('   0 = no limiter')
disp('   1 = minmod limiter')
disp('   2 = Koren limiter')
disp('   3 = Van Leer limiter')
disp('   4 = Van Albada limiter')
disp('   5 = Superbee limiter')
disp('  ')
while limiter<0 || limiter>5 ||rem(limiter,1)~=0
limiter = input('enter the number of the limiter you want to use: ');
end
end

%% animation
ani = 9;
disp('Choose type of animation during calculation')
disp('  ')
disp('   0 = no animation')
disp('   1 = waitbar')
disp('   2 = density plot')
disp('   3 = rho,u,p,e plots')
disp('   4 = density in tube')
disp('  ')
while ani<0 || ani>4 || rem(ani,1)~=0
ani = input('enter the number of animation you want to see:  ');
end