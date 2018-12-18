%% Script to evaluate the L2 norm for the Transcient solver for quadratic
% and linear cases
clear
clc
%Set Save Rate to limit Size of solution
SaveRate = 1; %Every time step is for good time resolution

%Quad Order
order = 2;  %order
dt = .0001;  % timestep
Tend = 1;   %run for until 1 second
Tvec = 0:dt:Tend; %Vector of timesteps
%Crank Scheme
theta = 1;
%% Create Boundary Condition Stucture
BC(1).type = "dirichlet";
BC(1).value = 0;
BC(2).type = "dirichlet";
BC(2).value = 1;

%% Caulculate L2 norm for different mesh resolutions
ne = [5, 8, 10, 12];    %The mesh sizes to used
for i = 1:length(ne)
    mesh = OneDimLinearMeshGen(0,1,ne(i),order);    %Create the mesh
    %% Set source constants
    mesh.fvec = 0*mesh.nvec;
    mesh.DCvec = ones(1,length(mesh.nvec));
    mesh.RCvec = 0*mesh.nvec;
    %Set Neumann Condition
    NBC = zeros(length(mesh.nvec), length(Tvec));  %No Neumann Condition
    %Set initial condition at t=0
    IC = zeros(order*mesh.ne +1, 1);
    %Calculate Solution
    SOL = TransReactDiffSolver(mesh, dt, Tend, theta, BC, IC, order, SaveRate);
    %Get solution at x = 0.9m
    C_FEM = SOL(:,(0.9/(dt*SaveRate))+1);
    %Calculate L2 norm
    L2(i) = L2norm(mesh, C_FEM, order, 0.9);
    %Calculate Logs
    L2log(i) = log(L2(i));
    nelog(i) = log(ne(i));
end
%% Calculate line of best fit
fitvars = polyfit(nelog, L2log, 1);
%Get the gradient of the line
grad = fitvars(1)
%Plot log log graph
figure()
loglog(L2,ne)

%% Format the Figure
title({'Quadratic L2 Norm at $x = 0.8m, \ t = 0.9s$','Time Step: 0.0001s'}, 'interpreter' ,'latex', 'FontSize', 14)
xlabel('Log of Numeber of Elements','interpreter','latex', 'FontSize', 12);
ylabel('Log of L2 Norm', 'interpreter','latex', 'FontSize', 12);
grid on
%% Saving Linear Basis Time Stepping Scheme Figure
%print -depsc \Users\xav_m\OneDrive\Documents\XAVI\University\Final_Year\Systems_Mod\Modeling_CW2\Report\Figures\epsQ1L2O2

%% Repeating for Order 1
order = 1;  %linear order
%% Caulculate L2 norm for different mesh resolutions
ne = [5, 8, 10, 12];    %The mesh sizes to used
for i = 1:length(ne)
    mesh = OneDimLinearMeshGen(0,1,ne(i),order);    %Create the mesh
    %% Set source constants
    mesh.fvec = 0*mesh.nvec;
    mesh.DCvec = ones(1,length(mesh.nvec));
    mesh.RCvec = zeros(1,length(mesh.nvec));
    %Set Neumann Condition
    NBC = zeros(length(mesh.nvec), length(Tvec));  %No Neumann Condition
    IC = zeros(order*mesh.ne +1, 1);    %Set initial condition at t=0
    %Calculate Solution
    SOL = TransReactDiffSolver(mesh, dt, Tend, theta, BC, IC, order, SaveRate);
    %Get solution at x = 0.9m
    C_FEM = SOL(:,(0.9/(dt*SaveRate))+1);
    %Calculate L2 norm
    L2(i) = L2norm(mesh, C_FEM, order, 0.9);
    %Calculate logs
    L2log(i) = log(L2(i));
    nelog(i) = log(ne(i));
end
%% Calculate line of best fit
fitvars2 = polyfit(nelog, L2log, 1);
%Get the gradient of the line
grad2 = fitvars2(1)
%Plot log log graph
figure()
loglog(L2,ne)

%% Format the Figure
title({'Linear L2 Norm at $x = 0.8m \ t = 0.9s$', 'Time Step: 0.0001s'}, 'interpreter' ,'latex', 'FontSize', 14)
xlabel('Log of Numeber of Elements','interpreter','latex', 'FontSize', 12);
ylabel('Log of L2 Norm', 'interpreter','latex', 'FontSize', 12);
grid on
%% Saving Linear Basis Time Stepping Scheme Figure
%print -depsc \Users\xav_m\OneDrive\Documents\XAVI\University\Final_Year\Systems_Mod\Modeling_CW2\Report\Figures\epsQ1L2O1
