%% This script generates the result for Part 2 Question 2a
clear
clc
%Create the mesh
order = 2;  %Linear order
xmin = 0;
xmax = 0.01;
Ne = 100;    % Number of elements
mesh = OneDimLinearMeshGen(xmin, xmax, Ne, order); %Create mesh

%Get Material and source coefficients
bloodflow = true; % no blood flow
mesh = setMatCoeffVectors(mesh, bloodflow, order);

%Set time stepping properties
theta = 0.5;     %Backward Euler Time Stepping Scheme
dt = 0.01;     %Time step
Tend = 50;
Tvec = 0:dt*5:Tend;
%Set Initial Condition
IC = 310.15 * ones((mesh.ne*order + 1), 1);
%% Create Boundary Condition Stucture
BC(1).type = "dirichlet";
BC(1).value = 329.3;
BC(2).type = "dirichlet";
BC(2).value = 310.15;

%Get bloodflow solution
%SOL = TransReactDiffSolver(mesh, dt, Tend, theta, BC, IC, order);
%Write result to csv file to save having to re-run
%dlmwrite('Q3profiles2.csv', SOL, 'delimiter', ',', 'precision', 17);

%Read in blood flow solution
SOL = csvread('Q3profiles2.csv');

%Read in NO blood flow solution
SOL_NoBlood = csvread('Q1aResult.csv');
figure()
%Plot skin boundaries
%Epidermis Boundary
plot([15/9000 15/9000], [310 410], '--k', 'HandleVisibility', 'off')
hold on
text(0.0002, 337, 'Epidermis', 'interpreter', 'latex', 'FontSize', 10)
%Dermis Boundary
plot([.005 .005], [310 410], '--k', 'HandleVisibility', 'off')
hold on
text(.003, 337, 'Dermis', 'interpreter', 'latex', 'FontSize', 10)
%Sub-cutaneous Region Text
text(.0065, 337, 'Sub-curtaneous', 'interpreter', 'latex', 'FontSize', 10)

%Plot Temperature Profiles Including Blood Flow
plot(mesh.nvec, SOL(:,(0.5/dt)/5), '-', 'color', [0.,0.55,0.55])  %Plot profile at t = 0.5s
hold on
plot(mesh.nvec, SOL(:,(1/dt)/5), '-b')    %Plot profile at t = 1s
hold on
plot(mesh.nvec, SOL(:,(2/dt)/5), '-m')    %Plot profile at t = 2s
hold on
plot(mesh.nvec, SOL(:,end), '-r')         %Plot profile at t =50s

%% Format the figure
title('Temperature Variation Through the Skin', 'interpreter', 'latex')
lgd = legend({'t = 0.5s', 't = 1s' , 't = 2s',  't = 50s'}, 'interpreter', 'latex', 'location', 'east');
xlabel('$x$ \ in \ meters','interpreter','latex', 'FontSize', 12);
ylabel('Temperature in Kelvin', 'interpreter','latex', 'FontSize', 12);
ylim([310 340]);
grid on
%% Saving The Figure
print -depsc \Users\xav_m\OneDrive\Documents\XAVI\University\Final_Year\Systems_Mod\Modeling_CW2\Report\Figures\epsQ3Profiles2
