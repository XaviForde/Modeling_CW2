function [ElemVector] = LinearSourceElemVector(mesh, eID, f_linear, Ngp, order)
% Returns the Local Linear Source Vector for an element
%
% Inputs:
% mesh - Mesh which contains local elements within it's structure and
%           related variables
% eID - index for the element within the mesh
% f_linear - Linear source term multiplier

%% Extract element variables from mesh
J =mesh.elem(1).J;      %Get Jacobi for the element
x0 = mesh.elem(eID).x(1);   %Get he x0 value for the element
x1 = mesh.elem(eID).x(2);   %Get he x1 value for the element
xVector = [x0;x1];          %Create a vector of the x values

%% Local Source Vector with J =1 and f_linear = 1
StandardMatrix = [(2/3), (1/3);...
                  (1/3), (2/3)];
              
%% Multiply by J, f_linear and the x vector to get Local Linear Source 
%Vector for the specific element
ElemVector = J * f_linear *  StandardMatrix * xVector;
end