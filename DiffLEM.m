function [DiffLEM] = DiffLEM(msh, Ngp, eID, order)
%% Returns the Local Element Diffusion Matrix
% Inputs:
%   mesh - 1D mesh containing mesh parameters (structure)
%   Ngp - Number of Gauss points for the Gauss scheme 
%   eID - The elements unique index
%   order - The order of Basis Functions (1 for Linear, 2 for Quadratic)

gq = CreateGQScheme(Ngp);
DiffLEM = zeros(order+1); %Initialize local element matrix

%% Evaluate local element matrix cells using Guassian Quadrature
for row = 1:(order+1)
    for col = 1:(order+1)
        for i = 1:Ngp
            %Get the value of the current Gauss point
            xipt = gq.xipts(i);
            %Evaluate Material Coefficient
            MatCoeff = EvalField(msh,msh.DCvec, eID,xipt,order);
            %Evaluate Gradient common for the current row
            row_grad = EvalBasisGrad(row-1,xipt,order);
            %Evaluate Gradient common for the current column
            col_grad = EvalBasisGrad(col-1,xipt,order);
            %Add integral for each Gauss point to the current cell
            DiffLEM(row,col) = DiffLEM(row,col) + MatCoeff*col_grad*row_grad;
        end
    end
end

%%Divide by the Jacobian
J = msh.elem(eID).J;
DiffLEM = DiffLEM/J;
end