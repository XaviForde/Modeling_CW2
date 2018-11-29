function [GlobalMass] = GlobalMassMatrixGen(mesh, Ngp, order)
%%Returns the Global mass matrix by iterating through elements and adding
% the local element matrices for diffusion and reaction into the global
% matrix in the correct position
% Inputs:
%   mesh -  1D mesh containing mesh parameters (structure)
%   Ngp - Number of Gauss points for the Gauss scheme
%   order - The order of Basis Functions (1 for Linear, 2 for Quadratic)


%% Initiate Global Matrix
ne = mesh.ne;       %Number of Elements
GlobalMass = zeros((order*ne +1),(order*ne +1));  
%% Loop over Elements and Assemble Local Element Matrices into Global Matrix
for eID = 1:ne
        
    %Local Element Matrix is the Diffusion subtract the Linear Reation
    MassLocal = MassLEM(mesh, Ngp, eID, order);
    
    %Find start and end index for the element
    StIdx = eID + (eID-1)*(order-1);
    EndIdx = StIdx + order;
    
    %Add Local Element Matrix into the correct location within the Global Matrix
    GlobalMass(StIdx:EndIdx, StIdx:EndIdx) = GlobalMass(StIdx:EndIdx, StIdx:EndIdx)...
                                            + MassLocal;
end