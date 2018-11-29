function [MassLEM] = MassLEM(mesh, Ngp, eID, order)
%% Returns the Local Element Mass Matrix
% Inputs:
%   mesh - 1D mesh containing mesh parameters (structure)
%   Ngp - Number of Gauss points for the Gauss scheme
%   eID - The elements unique index
%   order - The order of Basis Functions (1 for Linear, 2 for Quadratic)

%% Initialise Guass Scheme and Local Element Matrix
gq = CreateGQScheme(Ngp); %Create the Gauss scheme
MassLEM = zeros(order+1); %Initialize local element matrix

%% Evaluate local element matrix cells using Guassian Quadrature
for row = 1:(order+1)   %Iterate over rows
    
    for col = 1:(order+1)   %Iterate over the columns in the row
        
        %Use Guassian Quadrature to intergrate the eqn for the cell
        for i = 1:Ngp
            %Get the value of the current Gauss point
            xipt = gq.xipts(i);
            %Get the psi value common for the row
            psi_row = EvalBasis(row-1,xipt,order);
            %Get the psi value common for the column
            psi_col = EvalBasis(col-1,xipt,order);
            %Add integral for each Gauss point to the current cell
            MassLEM(row,col) = MassLEM(row,col) + psi_col*psi_row;
            
        end
    end
end

%% Multiply by Jacobian
J = mesh.elem(eID).J;    %Jaocbian
MassLEM = J*MassLEM;    %Local element matrix solution

end