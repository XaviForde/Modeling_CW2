function [F] = SourceVectorGen(mesh, Ngp, order)
%Generates source vector
% Inputs:
%   mesh -  1D mesh containing mesh parameters (structure)
%   Ngp - Number of Gauss points for the Gauss scheme
%   order - The order of Basis Functions (1 for Linear, 2 for Quadratic)

ne = mesh.ne;   %Number of elements

%% Initalise Gloabal Source Vector F
F = zeros(order*ne+1,1);     %Initialise Source Vector

%% Add in the Constant and Linear Element Vectors and Add to F
for eID = 1:ne
    
    %Caluculate Source Local Element Vectors
    LocalSourceVector = SourceLEVec(mesh, eID, Ngp, order);
    
    %Find start and end index for the element
    StIdx = eID + (eID-1)*(order-1);
    EndIdx = StIdx + order;
    
    %Add Local Source Vector into Global Source Vector at correct location
    F(StIdx:EndIdx) = F(StIdx:EndIdx) + LocalSourceVector;
end


