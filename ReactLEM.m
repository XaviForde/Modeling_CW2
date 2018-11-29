function [ReactLEM] = ReactLEM(mesh, Ngp, eID, order)
%% Returns the Local Element Reaction Matrix
% Inputs:
%   mesh - 1D mesh containing mesh parameters (structure)
%   Ngp - Number of Gauss points for the Gauss scheme 
%   eID - The elements unique index
%   order - The order of Basis Functions (1 for Linear, 2 for Quadratic)

%% Initiate Gauss Scheme and Rection LEM (Local Element Matrix)
gq = CreateGQScheme(Ngp);   %Gauss Scheme
ReactLEM = zeros(order+1);  %Initialize local element matrix

%% Evaluate local element matrix cells using Guassian Quadrature
for row = 1:(order+1)
    for col = 1:(order+1)
        for i = 1:Ngp
            %Get the value of the current Gauss point
            xipt = gq.xipts(i);
            %Get matierial coefficient
            MatCoeff = EvalField(mesh,mesh.RCvec, eID,xipt,order);
            %Get the psi value common for the current row
            psi_row = EvalBasis(row-1,xipt,order);
            %Get the psi value common for the current column
            psi_col = EvalBasis(col-1,xipt,order);
            %Add integral for each Gauss point to the current cell
            ReactLEM(row,col) = ReactLEM(row,col) + MatCoeff*psi_col*psi_row;
        end
    end
end

%% Multiple by the Jacobian
J = mesh.elem(eID).J;   %Jacobian 
ReactLEM = J*ReactLEM;  %Reaction LEM solution
end
