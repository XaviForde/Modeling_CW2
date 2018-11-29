%% Investigating the L2 nrom with repect to time


%% Set inital parameters
order = 1;  %linear order
dt = .01;   %timestep
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

%% Solve the Transcient Diffusion Reation equation using FEM method
SOL = TransReactDiffSolver(mesh, dt, Tend, theta, BC, IC, f_constant, f_linear, order);

sol5 = SOL(:,2);        %Result at T = 0.05
sol10 = SOL(:,3);       %Result at T = 0.10
sol20 = SOL(:,5);       %Result at T = 0.20
sol30 = SOL(:,7);       %Result at T = 0.30
sol40 = SOL(:,9);       %Result at T = 0.40
sol50 = SOL(:,11);       %Result at T = 0.50
sol60 = SOL(:,13);       %Result at T = 0.60
sol70 = SOL(:,15);       %Result at T = 0.70
sol80 = SOL(:,17);       %Result at T = 0.80
sol90 = SOL(:,19);       %Result at T = 0.90
sol100 = SOL(:,end);    %Result at T = 1

Tvec = 2:1:19;
L2_t = zeros(length(Tvec),1);
for i = 1:length(Tvec)
    C_FEM = SOL(:,Tvec(i));
    L2_t(i) = L2norm(mesh,C_FEM,1,Tvec(i));
end

plot(Tvec,L2_t)
%% Format the Figure
title({'Linear L2 Norm at $x = 0.8m$, Time Step: 0.01s'}, 'interpreter' ,'latex', 'FontSize', 12)
xlabel('Log of Numeber of Elements','interpreter','latex', 'FontSize', 12);
ylabel('Log of L2 Norm', 'interpreter','latex', 'FontSize', 12);
grid on
%% Saving Linear Basis Time Stepping Scheme Figure
print -depsc \Users\xav_m\OneDrive\Documents\XAVI\University\Final_Year\Systems_Mod\Modeling_CW2\Report\Figures\epsQ1L2Time

