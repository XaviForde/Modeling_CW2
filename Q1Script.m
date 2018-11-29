clear
clc

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
sol100 = SOL(:,end);    %Result at T = 1

tempx = mesh.nvec;

figure()
%% Plot FEM reults for linear basis functions:
plot(mesh.nvec, sol5(1:order:end), '--r','HandleVisibility','off')
hold on
plot(mesh.nvec, sol10(1:order:end), '--', 'color', [1 0.5 1],'HandleVisibility','off')
hold on
plot(mesh.nvec, sol20(1:order:end), '--', 'color', [0 0.5 0.5],'HandleVisibility','off')
hold on
plot(mesh.nvec, sol100(1:order:end), '--b','HandleVisibility','off')
hold on


%% Plot analytical solution:
%High resolution for analytical solution
Xvec = 0:0.001:Tend;
S5 = zeros(length(Xvec),1);
S10 = zeros(length(Xvec),1);
S20 = zeros(length(Xvec),1);
S100 = zeros(length(Xvec),1);

for i = 1:length(Xvec)
    S5(i) = TransientAnalyticSoln(Xvec(i), 0.05);
    S10(i) = TransientAnalyticSoln(Xvec(i), 0.1);
    S20(i) = TransientAnalyticSoln(Xvec(i), 0.2);
    S100(i) = TransientAnalyticSoln(Xvec(i), 0.9);
end
%Plotting the Analytical Results:
plot(Xvec, S5, '-r', 'HandleVisibility','on')
hold on
plot(Xvec, S10, '-', 'color',  [1 0.5 1], 'HandleVisibility','on')
hold on
plot(Xvec, S20, '-', 'color', [0 0.5 0.5], 'HandleVisibility','on')
hold on
plot(Xvec, S100, '-b', 'HandleVisibility','on')

%% Format the figure
title('Solution of $\frac{\delta c}{\delta t} = \frac{ \delta^2 c}{\delta x^2}$ Using Linear Basis Funcitons', 'interpreter' ,'latex', 'FontSize', 14)
lgd = legend({'t = 0.05', 't = 0.10', 't = 0.20', 't = 1'},'Location', 'northwest', 'interpreter', 'latex');
lgd.Title.String = '-- \ -- FEM \ , --- Analytical';
xlabel('$x$ \ in \ Metres','interpreter','latex', 'FontSize', 12);
ylabel('$c(x)$', 'interpreter','latex', 'FontSize', 12);
grid on
%% Saving Linear Basis Function Figure
%print -depsc \Users\xav_m\OneDrive\Documents\XAVI\University\Final_Year\Systems_Mod\Modeling_CW2\Report\Figures\epsQ1a

%% Now repeat for order 2 (quadratic)
order  = 2;
mesh = OneDimLinearMeshGen(0,1,ne,order);
%Set initial condition at t=0
IC = zeros(order*mesh.ne +1, 1);
%% Solve the Transcient Diffusion Reation equation using FEM method
SOLQ = TransReactDiffSolver(mesh, dt, Tend, theta, BC, IC, f_constant, f_linear, order);

sol5Q = SOLQ(:,2);        %Result at T = 0.05
sol10Q = SOLQ(:,3);       %Result at T = 0.10
sol20Q = SOLQ(:,5);       %Result at T = 0.20
sol100Q = SOLQ(:,end);    %Result at T = 1

%% Plot FEM reults for linear basis functions:
figure()
plot(mesh.nvec, sol5Q, '--r','HandleVisibility','off')
hold on
plot(mesh.nvec, sol10Q, '--', 'color', [1 0.5 1],'HandleVisibility','off')
hold on
plot(mesh.nvec, sol20Q, '--', 'color', [0 0.5 0.5],'HandleVisibility','off')
hold on
plot(mesh.nvec, sol100Q, '--b','HandleVisibility','off')
hold on


%Plotting the Analytical Results:
plot(Xvec, S5, '-r', 'HandleVisibility','on')
hold on
plot(Xvec, S10, '-', 'color',  [1 0.5 1], 'HandleVisibility','on')
hold on
plot(Xvec, S20, '-', 'color', [0 0.5 0.5], 'HandleVisibility','on')
hold on
plot(Xvec, S100, '-b', 'HandleVisibility','on')
%% Format the quadratic figure
title('Solution of $\frac{\delta c}{\delta t} = \frac{ \delta^2 c}{\delta x^2}$ Using Quadratic Basis Funcitons', 'interpreter' ,'latex', 'FontSize', 14)
lgd = legend({'t = 0.05', 't = 0.10', 't = 0.20', 't = 1'},'Location', 'northwest', 'interpreter', 'latex');
lgd.Title.String = '-- \ -- FEM \ , --- Analytical';
xlabel('$x$ \ in \ Metres','interpreter','latex', 'FontSize', 12);
ylabel('$c(x)$', 'interpreter','latex', 'FontSize', 12);
grid on
%% Saving Linear Basis Function Figure
%print -depsc \Users\xav_m\OneDrive\Documents\XAVI\University\Final_Year\Systems_Mod\Modeling_CW2\Report\Figures\epsQ1aQuad


