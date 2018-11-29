function [Gamma, EpiBndVec] = EvalGamma(T_x0, bloodflow)


%Create the mesh
order = 1;  %Linear order
xmin = 0;
xmax = 0.01;
Ne = 100;    % Number of elements
mesh = OneDimLinearMeshGen(xmin, xmax, Ne, order); %Create mesh

%Get Material and source coefficients

mesh = setMatCoeffVectors(mesh, bloodflow, order);

%Set time stepping properties
theta = 1;     %Backward Euler Time Stepping Scheme
dt = 0.01;     %Time step
Tend = 50;

%Set Initial Condition
IC = 310.15 * ones((mesh.ne*order + 1), 1);
%% Create Boundary Condition Stucture
BC(1).type = "dirichlet";
BC(1).value = T_x0;
BC(2).type = "dirichlet";
BC(2).value = 310.15;


SOL = TransReactDiffSolver(mesh, dt, Tend, theta, BC, IC, mesh.fvec, order);

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
% GammaBurnVec = zeros(1, length(BurnVec)); %Initialise
for i = 1:length(BurnVec)
    GammaBurnVec(i) = 2*10^98 *exp(-12017/(BurnVec(i)-273.15));
end

%Calculate gamma at tburn
%GammaTburn = zeros(1, length(BurnVec)); %Initialise
for i = 1:length(BurnVec)
    GammaTburn(i) = 2*10^98 *exp(-12017/(Tburn-273.15));
end

%Perform Integral including limits
Gamma = (trapz(GammaBurnVec) * dt) - (trapz(GammaTburn) * dt);

%Set gamma to zero if Tburn was never exceeded
if NoBurn == true
    Gamma = 0;
end
    
