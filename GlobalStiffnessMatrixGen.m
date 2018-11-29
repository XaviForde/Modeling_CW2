function [GlobalMatrix] = GlobalStiffnessMatrixGen(mesh, Ngp, order)
%%Returns the Global Stiffness Matrix by iterating through elements and adding
% the local element matrices for diffusion and reaction into the global
% matrix in the correct position
% Inputs:
%   mesh -  1D mesh containing mesh parameters (structure)
%   Ngp - Number of Gauss points for the Gauss scheme
%   order - The order of Basis Functions (1 for Linear, 2 for Quadratic)

%% Initiate Global Stiffness Matrix
ne = mesh.ne;       %Number of Elements
GlobalMatrix = zeros((order*ne+1),(order*ne+1));  %Global Stiffness Matrix
%% Loop over Elements and Assemble Local Element Matrices into Global Matrix
for eID = 1:ne
    
    %Calculate Local Element Matrix for Diffusion Term
    DiffusionLocal = DiffLEM(mesh, Ngp, eID, order);

    %Calculate Local Element Matrix for Linear Reaction Term
    ReactionLocal = ReactLEM(mesh, Ngp, eID, order);
    
    %Local Stiffness Element Matrix is the Diffusion subtract the Linear Reation
    StiffnessLEM = DiffusionLocal - ReactionLocal;
    
    %Find start and end index for the element
    StIdx = eID + (eID-1)*(order-1);
    EndIdx = StIdx + order;
    
    %Add Local Element Matrix into the correct location within the Global Matrix
    GlobalMatrix(StIdx:EndIdx, StIdx:EndIdx) = GlobalMatrix(StIdx:EndIdx, StIdx:EndIdx)...
                                            + StiffnessLEM;
end