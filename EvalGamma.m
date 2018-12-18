function [Gamma, EpiBndVec] = EvalGamma(T_x0, bloodflow)

%Set Save Rate to limit Size of solution
SaveRate = 5; %Every time step is for good time resolution

%Create the mesh
order = 1;  %Linear order
xmin = 0;
xmax = 0.01;
Ne = 100;    % Number of elements
mesh = OneDimLinearMeshGen(xmin, xmax, Ne, order); %Create mesh

%Get Material and source coefficients

mesh = setMatCoeffVectors(mesh, bloodflow, order);

%Set time stepping properties
theta = 0.5;     %Backward Euler Time Stepping Scheme
dt = 0.01;     %Time step
Tend = 50;
Tvec = 0:dt:Tend; %Vector of timesteps

%Set Initial Condition
IC = 310.15 * ones((mesh.ne*order + 1), 1);
%% Create Boundary Condition Stucture
BC(1).type = "dirichlet";
BC(1).value = T_x0;
BC(2).type = "dirichlet";
BC(2).value = 310.15;
NBC = zeros(length(mesh.nvec), length(Tvec));  %No Neumann Condition

SOL = TransReactDiffSolver(mesh, dt, Tend, theta, BC, NBC, IC, order, SaveRate);

%locate epidermis boundary
i = 1;  %start ID counter
while mesh.DCvec(i) == mesh.DCvec(i+1)
    i = i+1;
end
EpiBndID = i; %Epidermis Boundary index
EpiBndVec = SOL(EpiBndID, :);   %Epidermis Boundry Temp Profile for plots
%locate timestep where burn temperatue is reached
Tburn = 317.15;
i = 1;
%Creat logical NoBurn to deal with cases where Tburn is not exceeded
NoBurn = false; %Assume Tburn exceeded

while SOL(EpiBndID,i) < Tburn
    i = i+1;
    if i+1 > length(mesh.DCvec) %Tburn is never exceeded
        NoBurn = true;  %will use this to set gamma to zero
        i = length(mesh.DCvec)-1;   %Set to minimise integration time
        break
    end
end

BurnStID = i;

%Get Vector of temperatures after burning temp is reached
BurnVec = SOL(EpiBndID, BurnStID:end);

%Calculate gamma for each temperature in the burn vector
GammaBurnVec = zeros(1, length(BurnVec)); %Initialise
for i = 1:length(BurnVec)
    GammaBurnVec(i) = 2*10^98 *exp(-12017/(BurnVec(i)-273.15));
end

%Perform Integral including limits
Gamma = (trapz(GammaBurnVec) * dt);

%Set gamma to zero if Tburn was never exceeded
if NoBurn == true
    Gamma = 0;
end
    

