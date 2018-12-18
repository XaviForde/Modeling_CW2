%% Script to plot the result for question 1b
clear
clc
color1 = [0.4 0.1 0.6]; %Color for plotting
%% Set inital parameters
order = 1;  %linear order
dt = 0.1;   %timestep
Tend = 1;   %run for until 1 second
Tvec = 0:dt:Tend; %Vector of timesteps
SaveRate = 1;   %Save every time ttep to see oscilations
%% Create 10 element mesh
ne = 10;
mesh = OneDimLinearMeshGen(0,1,ne,order);
%% Set source constants
mesh.fvec = 0*mesh.nvec;    %No source
mesh.DCvec = ones(1,length(mesh.nvec)); %D = 1 everywhere
mesh.RCvec = zeros(1,length(mesh.nvec)); % No reaction
%% Time Stepping Scheme
theta = 0.5;    % Crank-Nicolson Scheme
%% Create Boundary Condition Stucture
BC(1).type = "dirichlet";   
BC(1).value = 0;
BC(2).type = "dirichlet";
BC(2).value = 1;
NBC = zeros(length(mesh.nvec), length(Tvec));  %No Neumann Condition
%Set initial condition at t=0
IC = zeros(order*mesh.ne +1, 1);

%% Solve Using Crank-Nicolson Time Stepping Scheme
SOL_CN = TransReactDiffSolver(mesh, dt, Tend, theta, BC, NBC, IC, order, SaveRate);

%% Solve using Backward Euler Timestepping
theta = 1;  %Backward Euler scheme
SOL_BE = TransReactDiffSolver(mesh, dt, Tend, theta, BC, NBC, IC, order, SaveRate);

%% Solve using the Analytical Solution
TvecAn = 0:0.01:Tend;
C80_Analytic = zeros(length(TvecAn), 1);
for i = 1:length(TvecAn)
    C80_Analytic(i) = TransientAnalyticSoln(0.8, TvecAn(i));
end
    
%% Plot how temperature varies at x = 0.8
Tvec = 0:dt:Tend;
C80_CN = SOL_CN(8*order+1, :);  %Crank-Nicolson at x = 0.80
C80_BE = SOL_BE(8*order+1, :);  %Backward Euler at x = 0.80
figure()
plot(TvecAn, C80_Analytic, '-')     %Plot Analytic Solution
hold on
plot(Tvec(2:SaveRate:end),C80_CN(2:end), '-*r')     %Plot Crank-Nicolson
hold on
plot(Tvec(2:SaveRate:end),C80_BE(2:end), '-x', 'color', color1)     %Plot Backward Euler

%% Format the Figure
title({'Solution of $\frac{\delta c}{\delta t} = \frac{ \delta^2 c}{\delta x^2}$ at $x = 0.8m$, Time Step: 0.1s'}, 'interpreter' ,'latex', 'FontSize', 14)
legend({'Analytical', 'Crank-Nicolson', 'Backward Euler'},'Location', 'southeast', 'interpreter', 'latex');
xlabel('Time \ in \ Seconds','interpreter','latex', 'FontSize', 12);
ylabel('$c(0.8,t)$', 'interpreter','latex', 'FontSize', 12);
grid on
%% Saving 0.1s Time Step Figure
%print -depsc \Users\xav_m\OneDrive\Documents\XAVI\University\Final_Year\Systems_Mod\Modeling_CW2\Report\Figures\epsQ1bdt05


%% %%%%  Repeat with smaller timestep 
dt = .02;   %New timestep
Tvec = 0:dt:Tend; %Vector of timesteps
NBC = zeros(length(mesh.nvec), length(Tvec));  %No Neumann Condition
%% Solve Using Crank-Nicolson Time Stepping Scheme
theta = 0.5;    %Crank-Nicolson Scheme
SOL_CN = TransReactDiffSolver(mesh, dt, Tend, theta, BC, NBC, IC, order, SaveRate);

%% Solve using Backward Euler Timestepping
theta = 1;  %Backward Euler scheme
SOL_BE = TransReactDiffSolver(mesh, dt, Tend, theta, BC, NBC, IC, order, SaveRate);

%% Plot how temperature varies at x = 0.8
Tvec = 0:dt:Tend;   %Time Vector for FEM solutions
C80_CN = SOL_CN(8*order+1, :);  %Crank-Nicolson at x = 0.80
C80_BE = SOL_BE(8*order+1, :);  %Backward Euler at x = 0.80
figure()
plot(TvecAn, C80_Analytic, '-')     %Plot Analytic Solution
hold on
plot(Tvec(2:SaveRate:end),C80_CN(2:end), '-r')     %Plot Crank-Nicolson
hold on
plot(Tvec(2:SaveRate:end),C80_BE(2:end), '-', 'color', color1)     %Plot Backward Euler

%% Format the Figure
title({'Solution of $\frac{\delta c}{\delta t} = \frac{ \delta^2 c}{\delta x^2}$ at $x = 0.8m$, Time Step: 0.02s'}, 'interpreter' ,'latex', 'FontSize', 14)
legend({'Analytical', 'Crank-Nicolson', 'Backward Euler'},'Location', 'southeast', 'interpreter', 'latex');
xlabel('Time \ in \ Seconds','interpreter','latex', 'FontSize', 12);
ylabel('$c(0.8,t)$', 'interpreter','latex', 'FontSize', 12);
ylim([0 0.9])
grid on
%% Saving 0.02s Time Step Figure
%print -depsc \Users\xav_m\OneDrive\Documents\XAVI\University\Final_Year\Systems_Mod\Modeling_CW2\Report\Figures\epsQ1bdt025
