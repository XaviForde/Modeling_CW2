clear
clc

%% This script generates the figure for question 1b - T at x = 0.8m

%Set Save Rate to limit Size of solution
SaveRate = 1; %Every time step is for good time resolution

%Solve the FEM solution:
%% Set inital parameters
order = 2;  %quadratic order
dt = 0.01;   %timestep
Tend = 1;   %run for until 1 second
Tvec = 0:dt:Tend; %Vector of timesteps
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



%% Solve the Transcient Diffusion Reation equation using FEM method
SOL = TransReactDiffSolver(mesh, dt, Tend, theta, BC, NBC, IC, order, SaveRate);


%% Solve using the Analytical Solution:
TvecAn = 0:0.01:Tend;
C80_Analytic = zeros(length(TvecAn), 1);
for i = 1:length(TvecAn)
    C80_Analytic(i) = TransientAnalyticSoln(0.8, TvecAn(i));
end

%% Plot how temperature varies at x = 0.8
Tvec = 0:dt*SaveRate:Tend;  %Set time vector for analytical plot
C80 = SOL(8*order+1, :);  %Crank-Nicolson at x = 0.80
fig = figure();
plot(TvecAn, C80_Analytic, '-')     %Plot Analytic Solution
hold on
plot(Tvec(2:end),C80(2:end), '-r')     %Plot Crank-Nicolson

%% Format the Figure
title({'Solution of $\frac{\delta c}{\delta t} = \frac{ \delta^2 c}{\delta x^2}$ at $x = 0.8m$, Time Step: 0.01s'}, 'interpreter' ,'latex', 'FontSize', 14)
lgd = legend({'Analytical', 'Crank-Nicolson'},'Location', 'southeast', 'interpreter', 'latex');
lgd.Title.String = '-- \ -- FEM \ , --- Analytical';
xlabel('Time \ in \ Seconds','interpreter','latex', 'FontSize', 12);
ylabel('$c(0.8,t)$', 'interpreter','latex', 'FontSize', 12);
ylim([0 0.9])
grid on
%% Saving Linear Basis Time Stepping Scheme Figure
%print -depsc \Users\xav_m\OneDrive\Documents\XAVI\University\Final_Year\Systems_Mod\Modeling_CW2\Report\Figures\epsQ1b

%Creat and save a zoomed plot
figzoom = gcf;
ylim auto
xlim([0.05 0.052])
%% Saving Zoomed Figure
%print -depsc \Users\xav_m\OneDrive\Documents\XAVI\University\Final_Year\Systems_Mod\Modeling_CW2\Report\Figures\epsQ1bZoom



