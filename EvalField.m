function [mat_value] = EvalField(mesh, MatVec, eID, xipt, order)
%% Returns material property for an element
%Could be diffusion or reaction coefficient that is returned.
% Inputs:
%   mesh - 1D mesh containing mesh parameters (structure)
%   MatVec - Vector containing material property values at each node
%   eID - The elements unique index
%   xipt - Gauss point
%   order - The order of Basis Function (1 for Linear, 2 for Quadratic)


switch order
    case 1  %Linear Basis Function
        psi0 = EvalBasis(0, xipt, order);   %basis function value at x0
        psi1 = EvalBasis(1, xipt, order);   %Basis function value at x1
        psi = [psi0, psi1];
        coeffs = [MatVec(eID); MatVec(eID+1)];  %Material coeffs at nodes
        mat_value = psi*coeffs; %Coefficient for the element
        
    case 2  %Quadratic Basis Function
        psi = [EvalBasis(0, xipt, order) EvalBasis(1, xipt, order) EvalBasis(2, xipt, order)];
        coeffs = [MatVec(2*eID-1); MatVec(2*eID); MatVec(2*eID+1)]; %Material coeffs at nodes
        mat_value = psi*coeffs; %Coefficient for the element
end