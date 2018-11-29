function [ psi ] = EvalBasis(lnid,xipt,order)
%Returns the value of the Laplace Basis functions at a given value of x
% Inputs:
%   lnid - Local node index 
%   xipt - Gauss point
%   order - The order of Basis Function (1 for Linear, 2 for Quadratic)

switch order
    case 1  %Linear Basis Function
        sign = (-1)^(lnid+1);
        psi = (1 + sign*xipt)/2;
        
    case 2  %Quadratic Basis Function
        switch lnid
            case 0
                psi = xipt*(xipt-1) / 2;    %psi0
            case 1
                psi = 1 - xipt^2;           %psi1
            case 2
                psi = xipt*(1+xipt) / 2;    %psi2
        end
        
end

