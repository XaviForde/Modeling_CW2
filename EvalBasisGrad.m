function [ dpsidxi ] = EvalBasisGrad(lnid,xipt, order)
%EvalBasisGrad Returns gradient of basis functions
% Inputs:
%   lnid - Local node index 
%   xipt - Gauss point
%   order - The order of Basis Function (1 for Linear, 2 for Quadratic)

%Use the node id to generate the sign of the basis gradient - ie.
%either + or -. when lnid=0, sign is -ve, when lnid=1, sign is +ve.
switch order
    case 1  %Linear basis function
        sign = (-1)^(lnid+1);
        dpsidxi = 0.5 * sign;
       
    case 2  %Quadrtatic basis function
        switch lnid
            case 0
                dpsidxi = xipt - 0.5;
            case 1
                dpsidxi = -2*xipt;
            case 2
                dpsidxi = xipt + 0.5;                 
        end
end


end

