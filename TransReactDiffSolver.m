function [Ct_matrix] = TransReactDiffSolver(mesh, dt, Tend, theta, BC, IC, F, order)
%% This function solves the transient diffusion-reation equation using the
%Finite Element Method
% Inputs:
%   mesh -  1D mesh containing mesh information (structure)
%   dt -    Timestep
%   Tend -  Function will calculate solution until this time is reached
%   theta - Method of timestepping (0 for forward Euler, 0.5 for
%           Crank-Niclson, 1 for backward Euler)
%   BC -    Boundary Conditions (structure)
%   IC -    Initial Condition (vector)
%   order - The order of Basis Functions (1 for Linear, 2 for Quadratic)

%% Set Number of Guass points
Ngp = order+1;

%% Initalise time
Tvec = 0:dt:Tend;

%% Initial condition, t = 0
Cnext = IC;
Ct_matrix(:,1) = Cnext;

%% Iterate through remaining timesteps
for i = 2:length(Tvec)
    
    %Set most previous time step solution to current solution
    Ccurrent = Cnext;
    % Calculate Stiffness and Mass Global Matrices at this timestep
    K = GlobalStiffnessMatrixGen(mesh, Ngp, order);  %Stiffness
    M = GlobalMassMatrixGen(mesh, Ngp, order);  %Mass
    GM = M + theta*dt*K;
    tempM = (M - (1-theta)*dt*K);
    GVec = tempM*Ccurrent;
    GVec = GVec + dt*(theta*F + (1-theta)*F);
    
    %Apply Boundary Conditions
    [GM, GVec] = ApplyBCs(BC, GM, GVec, mesh.ne, order);
    %Calculate next C solution vector
    Cnext = GM\GVec;
    
    
    %% Save Every 5th timestep in matlab matrix
    save_every = 5;
    if mod(i-1,save_every) == 0
        Ct_matrix(:,((i-1)/save_every)+1) = Cnext;
    end
    
end
end
