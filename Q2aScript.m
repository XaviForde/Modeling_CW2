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
bloodflow = false; % no blood flow
mesh = setMatCoeffVectors(mesh, bloodflow, order);

%Set time stepping properties
theta = 1;     %Backward Euler Time Stepping Scheme
dt = 0.01;     %Time step
Tend = 50;

%Set Initial Condition
IC = 310.15 * ones((mesh.ne*order + 1), 1);
%% Create Boundary Condition Stucture
BC(1).type = "dirichlet";   
BC(1).value = 393.15;
BC(2).type = "dirichlet";
BC(2).value = 310.15;


%SOL = TransReactDiffSolver(mesh, dt, Tend, theta, BC, IC, mesh.fvec, order);
%Write result to csv file to save having to re-run
%dlmwrite('Q1aResult.csv', SOL, 'delimiter', ',', 'precision', 17);
SOL = csvread('Q1aResult.csv');
figure()
%Plot skin boundaries
%Epidermis Boundary
plot([15/9000 15/9000], [310 410], '--k', 'HandleVisibility', 'off')
hold on
text(0.0002, 405, 'Epidermis', 'interpreter', 'latex', 'FontSize', 10)
%Dermis Boundary
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


%% Calculate Gamma
%locate epidermis boundary
i = 1;  %start ID counter
while mesh.DCvec(i) == mesh.DCvec(i+1)
    i = i+1;
end
EpiBndID = i; %Epidermis Boundary index

%locate timestep where burn temperatue is reached
Tburn = 317.15;
i = 1;

while SOL(EpiBndID,i) < Tburn
    i = i+1;
end

BurnStID = i;

%Get Vector of temperatures after burning temp is reached
BurnVec = SOL(EpiBndID, BurnStID:end);

%Calculate gamma for each temperature in the burn vector

for i = 1:length(BurnVec)
    dGammadt(i) = 2*10^98 *exp(-12017/(BurnVec(i)-273.15));
end

gamma = trapz(dGammadt) * dt - Tburn*dt*(length(BurnVec)-1);



