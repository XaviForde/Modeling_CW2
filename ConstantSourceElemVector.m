function [ConstSourceVec] = ConstantSourceElemVector(mesh, eID, f_constant, Ngp, order)
% Returns the Local Constant Source Vector for an element
%
% Inputs:
% mesh - Mesh which contains local elements within it's structure and
%           related variables
% eID - index for the element within the mesh
% f_const - Constant source term 

J =mesh.elem(eID).J;        %Get Jacobi for the element

%ElemVector = [f_constant*J; f_constant*J];  %each term is fJ

%%%%Could take FieldVal out of the loop and multiply at end??
gq = CreateGQScheme(Ngp);
N = gq.npts;
ConstSourceVec = zeros(order+1,1); %Initialize local element matrix

%% Create local element, evaluate using Guassian Quadrature
for row = 1:(order+1)
        for i = 1:N
            xipt = gq.xipts(i);
            ConstSourceVec(row,1) = ConstSourceVec(row,1) + EvalBasis(row-1,xipt,order);
        end
end

%% Multiply by J and f

ConstSourceVec = J*f_constant*ConstSourceVec;

end
