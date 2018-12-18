function [LocalSourceVec] = SourceLEVec(mesh, eID, Ngp, order)
% Returns the Local Element Vector (LEVec) for the source term
%
% Inputs:
%   mesh - 1D mesh containing mesh parameters (structure)
%   Ngp - Number of Gauss points for the Gauss scheme 
%   eID - The elements unique index
%   order - The order of Basis Functions (1 for Linear, 2 for Quadratic)

gq = CreateGQScheme(Ngp);
LocalSourceVec = zeros(order+1,1); %Initialize local element matrix

%% Create local element, evaluate using Guassian Quadrature
for row = 1:(order+1)
    for i = 1:Ngp
        %Get the value of the current Gauss point
        xipt = gq.xipts(i);
        %Get source coefficient
        SourceCoeff = EvalField(mesh, mesh.fvec, eID,xipt,order);
        %Get the psi value common for the current row
        psi_row = EvalBasis(row-1,xipt,order);
        %Add each integration into loval vector
        LocalSourceVec(row,1) = LocalSourceVec(row,1) + SourceCoeff*psi_row;
    end
end

%% Multiply by J and f
J =mesh.elem(eID).J;        %Get Jacobi for the element
LocalSourceVec = J*LocalSourceVec;

end
