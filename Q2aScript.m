%% This script generates the result for Part 2 Question 2a
clear
clc

%Set Save Rate to limit Size of solution
SaveRate = 5; %Every fith time step is sufficient resolution


%% Set up the problem
%Create mesh
order = 2;  %quadratic order
xmin = 0;
xmax = 0.01;
Ne = 100;    % Number of elements
mesh = OneDimLinearMeshGen(xmin, xmax, Ne, order); %Create mesh

%Add Material and source coefficients to mesh
bloodflow = false; % no blood flow
mesh = setMatCoeffVectors(mesh, bloodflow, order); 

%Set time stepping properties
theta = 1;     %Backward Euler Time Stepping Scheme
dt = 0.005;     %Time step
Tend = 50;  %Run over 50s domain
Tvec = 0:dt:Tend; %Vector of timesteps
%Set Initial Condition
IC = 310.15 * ones((mesh.ne*order + 1), 1);
%% Create Boundary Condition Stucture
BC(1).type = "dirichlet";
BC(1).value = 393.15;
BC(2).type = "dirichlet";
BC(2).value = 310.15;
NBC = zeros(length(mesh.nvec), length(Tvec));  %No Neumann Condition

%SOL = TransReactDiffSolver(mesh, dt, Tend, theta, BC, NBC, IC, order, SaveRate);
%Write result to csv file to save having to re-run:
%dlmwrite('Q1aResult.csv', SOL, 'delimiter', ',', 'precision', 17);
%Read in previous result so don't have to run whilst formatting
SOL = csvread('Q1aResult.csv');

%% Plot Temperature Profiles
figure()
%Plot skin boundary and regions
%Epidermis Boundary and text
plot([15/9000 15/9000], [310 410], '--k', 'HandleVisibility', 'off')
hold on
text(0.0002, 405, 'Epidermis', 'interpreter', 'latex', 'FontSize', 10)
%Dermis Boundary and text
plot([.005 .005], [310 410], '--k', 'HandleVisibility', 'off')
hold on
text(.003, 405, 'Dermis', 'interpreter', 'latex', 'FontSize', 10)
%Sub-cutaneous Region Text
text(.0065, 405, 'Sub-curtaneous', 'interpreter', 'latex', 'FontSize', 10)

%Plot Temperature Profiles
plot(mesh.nvec, SOL(:,(0.5/dt)/5))  %Plot profile at t = 0.5s
hold on
plot(mesh.nvec, SOL(:,(1/dt)/5))    %Plot profile at t = 1s
hold on
plot(mesh.nvec, SOL(:,(2/dt)/5))    %Plot profile at t = 2s
hold on
plot(mesh.nvec, SOL(:,(5/dt)/5))    %Plot profile at t = 5s
hold on
plot(mesh.nvec, SOL(:,end))         %Plot profile at t =50s

%% Format the figure
title('Temperature Variation Through the Skin', 'interpreter', 'latex')
legend({'t = 0.5s', 't = 1s' , 't = 2s', 't = 5s', 't = 50s'}, 'interpreter', 'latex', 'location', 'east')
xlabel('$x$ \ in \ meters','interpreter','latex', 'FontSize', 12);
ylabel('Temperature in Kelvin', 'interpreter','latex', 'FontSize', 12);
ylim([310 410]);
grid on
%% Saving The Figure
%print -depsc \Users\xav_m\OneDrive\Documents\XAVI\University\Final_Year\Systems_Mod\Modeling_CW2\Report\Figures\epsQ21TempXT

%% Plotting variation of T at single depth point with time
figure()
plot(0:.005*SaveRate:Tend, SOL(41,:))
hold on
plot(0:.005*SaveRate:Tend, SOL(101,:))
%% Format the figure
title('Temperature Variation With Time', 'interpreter', 'latex')
legend({'x = 0.002m', 'x = 0.005m'}, 'interpreter', 'latex', 'location', 'best')
xlabel('Time  in  seconds','interpreter','latex', 'FontSize', 12);
ylabel('Temperature in Kelvin', 'interpreter','latex', 'FontSize', 12);
ylim([310 410]);
grid on
%print -depsc \Users\xav_m\OneDrive\Documents\XAVI\University\Final_Year\Systems_Mod\Modeling_CW2\Report\Figures\epsTimePlot

%% Create surf plot
figure()
TimeVec = 0:dt*SaveRate:Tend;
%Space out time plots to improve clarity of figure
Tidx = [1:4:41 , 47:6:101, 113:12:221 ,261:40:length(TimeVec)];
%Get X, Y and Z 
Xsurf = mesh.nvec(1:5:end);
Ysurf = TimeVec(Tidx);
Zsurf = SOL(1:5:end, Tidx)';
%Plot
surf(Xsurf , Ysurf,  Zsurf)
colormap jet
%% Format the figure
title('Temperature Variation Through the Skin Over Time', 'interpreter', 'latex')
xlabel('$x$ \ in \ meters','interpreter','latex', 'FontSize', 12);
ylabel('Time  in  seconds','interpreter','latex', 'FontSize', 12);
zlabel('Temperature in Kelvin', 'interpreter','latex', 'FontSize', 12);
%Set Limits
zlim([310 390])
ylim([0 50])
xlim([0,0.01])
% grid on
%print -depsc \Users\xav_m\OneDrive\Documents\XAVI\University\Final_Year\Systems_Mod\Modeling_CW2\Report\Figures\epsSurf1


