function [GlobalMatrix, F] = ApplyBCs(BC, GlobalMatrix, F, ne, order)
% Applies boundary conditions to F and Global Matrix as appropriate
% for the specified boundary condition
%
%Inputs:
%BC - Structure which contains the following:
%       BC(1) - Holds data for minimum x boundary
%       BC(2) - Holds data for maximum x boundary
%       BC().type - Type of boundary condition: "neumann", "dirichlet" or
%                   "none" (must be a lower case string)
%       BC().value - Value of the boundary condition (float or int)
%
%GlobalMatrix - formation of local element matrices (NxN matrix)
%F - Source Vector (size Nx1 vector)

%% If user inputted character arrays change these to strings
BC(1).type = string(BC(1).type);
BC(2).type = string(BC(2).type);

%% Solve xmin Boundary Condition 
if BC(1).type == "none" %No action needed if no boundary condition
    % pass 
    
elseif BC(1).type == "dirichlet"  %solve for a dirichlet BC
    
    %set all first row elements to 0 except first element  
    GlobalMatrix(1,:) = [1, zeros(1,ne*order)];   
    %set first element of the source vector to the xmin BC value
    F(1) = BC(1).value;
    
elseif  BC(1).type == "nuemann"    %solve for a nuemann BC
    
    %subtract the xmin BC value from the first element of the source vector 
    F(1) = F(1) - BC(1).value;
end

%% Solve xmax Boundary Condition
if BC(2).type == "none" %No action needed if no boundary condition
    % pass
    
elseif BC(2).type == "dirichlet"  %solve for dirichlet BC
    
    %set all end row elements to 0 except last element 
    GlobalMatrix(end,:) = [zeros(1,ne*order), 1];
    %set last element of the source vector to the xmax BC value
    F(end) = BC(2).value;
    
elseif BC(2).type == "nuemann"    %solve nuemann BC
    
    %subtract the xmin BC value from the last element of the source vector 
    F(end) = F(end) + BC(2).value;
end   
end




