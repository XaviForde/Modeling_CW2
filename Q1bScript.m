%% Script to plot the result for question 1b

clear
clc
color1 = [0.4 0.1 0.6];
%% Set inital parameters
order = 1;  %linear order
dt = .1;   %timestep
Tend = 1;   %run for until 1 second
%% Create 10 element mesh
ne = 10;
mesh = OneDimLinearMeshGen(0,1,ne,order);
%% Crank-Nicolson Scheme
theta = 0.5;
%% Create Boundary Condition Stucture
BC(1).type = "dirichlet";   
BC(1).value = 0;
BC(2).type = "dirichlet";
BC(2).value = 1;

%Set initial condition at t=0
IC = zeros(order*mesh.ne +1, 1);

%Set source constants
f_constant = 0;
f_linear = 0;

%% Solve Using Crank-Nicolson Time Stepping Scheme
SOL_CN = TransReactDiffSolver(mesh, dt, Tend, theta, BC, IC, f_constant, f_linear, order);

%% Solve using Backward Euler Timestepping
theta = 1;  %Backward Euler scheme
SOL_BE = TransReactDiffSolver(mesh, dt, Tend, theta, BC, IC, f_constant, f_linear, order);

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
plot(Tvec(2:end),C80_CN(2:end), '-*r')     %Plot Crank-Nicolson
hold on
plot(Tvec(2:end),C80_BE(2:end), '-x', 'color', color1)     %Plot Backward Euler

%% Format the Figure
title({'Solution of $\frac{\delta c}{\delta t} = \frac{ \delta^2 c}{\delta x^2}$ at $x = 0.8m$, Time Step: 0.1s'}, 'interpreter' ,'latex', 'FontSize', 14)
lgd = legend({'Analytical', 'Crank-Nicolson', 'Backward Euler', 't = 1'},'Location', 'southeast', 'interpreter', 'latex');
%lgd.Title.String = '-- \ -- FEM \ , --- Analytical';
xlabel('Time \ in \ Seconds','interpreter','latex', 'FontSize', 12);
ylabel('$c(0.8,t)$', 'interpreter','latex', 'FontSize', 12);
grid on
%% Saving Linear Basis Time Stepping Scheme Figure
print -depsc \Users\xav_m\OneDrive\Documents\XAVI\University\Final_Year\Systems_Mod\Modeling_CW2\Report\Figures\epsQ1bdt05


%% %%%%  Repeat plot with reduced timestep %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = .02;   %New timestep
%% Solve Using Crank-Nicolson Time Stepping Scheme
theta = 0.5;    %Crank-Nicolson Euler scheme
SOL_CN = TransReactDiffSolver(mesh, dt, Tend, theta, BC, IC, f_constant, f_linear, order);

%% Solve using Backward Euler Timestepping
theta = 1;  %Backward Euler scheme
SOL_BE = TransReactDiffSolver(mesh, dt, Tend, theta, BC, IC, f_constant, f_linear, order);

%% Plot how temperature varies at x = 0.8
Tvec = 0:dt:Tend;   %Time Vector for FEM solutions
C80_CN = SOL_CN(8*order+1, :);  %Crank-Nicolson at x = 0.80
C80_BE = SOL_BE(8*order+1, :);  %Backward Euler at x = 0.80
figure(4)
plot(TvecAn, C80_Analytic, '-')     %Plot Analytic Solution
hold on
plot(Tvec(2:end),C80_CN(2:end), '-r')     %Plot Crank-Nicolson
hold on
plot(Tvec(2:end),C80_BE(2:end), '-', 'color', color1)     %Plot Backward Euler

%% Format the Figure
title({'Solution of $\frac{\delta c}{\delta t} = \frac{ \delta^2 c}{\delta x^2}$ at $x = 0.8m$, Time Step: 0.02s'}, 'interpreter' ,'latex', 'FontSize', 14)
lgd = legend({'Analytical', 'Crank-Nicolson', 'Backward Euler'},'Location', 'southeast', 'interpreter', 'latex');
%lgd.Title.String = '-- \ -- FEM \ , --- Analytical';
xlabel('Time \ in \ Seconds','interpreter','latex', 'FontSize', 12);
ylabel('$c(0.8,t)$', 'interpreter','latex', 'FontSize', 12);
grid on
%% Saving Linear Basis Time Stepping Scheme Figure
print -depsc \Users\xav_m\OneDrive\Documents\XAVI\University\Final_Year\Systems_Mod\Modeling_CW2\Report\Figures\epsQ1bdt025
