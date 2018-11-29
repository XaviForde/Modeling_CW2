%% Script to evaluate the L2 norm for the Transcient solver
clear
clc

%Quad Order
order = 2;  %order
dt = .0001;  % timestep
Tend = 1;   %run for until 1 second
%Crank Scheme
theta = 1;
%% Create Boundary Condition Stucture
BC(1).type = "dirichlet";
BC(1).value = 0;
BC(2).type = "dirichlet";
BC(2).value = 1;

%Set source constants
f_constant = 0;
f_linear = 0;

%% Caulculate L2 norm for different mesh resolutions
ne = [5, 8, 10, 12];    %The mesh sizes to used
for i = 1:length(ne)
    mesh = OneDimLinearMeshGen(0,1,ne(i),order);    %Create the mesh
    %Set initial condition at t=0
    IC = zeros(order*mesh.ne +1, 1);
    SOL = TransReactDiffSolver(mesh, dt, Tend, theta, BC, IC, f_constant, f_linear, order);
    C_FEM = SOL(:,(0.9/(dt*5))+1);
    L2(i) = L2norm(mesh, C_FEM, order, 0.9);
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
    IC = zeros(order*mesh.ne +1, 1);    %Set initial condition at t=0
    SOL = TransReactDiffSolver(mesh, dt, Tend, theta, BC, IC, f_constant, f_linear, order);
    C_FEM = SOL(:,(0.9/(dt*5))+1);
    L2(i) = L2norm(mesh, C_FEM, order, 0.9);
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
