function [L2norm] = L2norm(mesh, C_FEM, order, t)
%% Calculates the L2 norm error for a vector
% Inputs:
%   mesh - A mesh containing various information (structure)
%   C_FEM - The numerical solution (vector)
%   order - The order of Basis Functions (1 for Linear, 2 for Quadratic)
%   t - The time at which the L2 norm is being calculated for (seconds)


%% Set up Gauss Scheme
Ngp = order+2; %Number of gauss ponts
gq = CreateGQScheme(Ngp);   %Creat Gauss Scheme
Nel = mesh.ne; %Number of elements
%% Initialise error to 0
error = 0;

%% Integrate using Gaussian Quadrature
switch order
    case 1  %% Solve for Linear Basis Functions
        for eID = 1:Nel
            
            %Get Jacobian for the element
            J = mesh.elem(eID).J;
            %x values element nodes
            x0 = mesh.elem(eID).x(1);
            x1 = mesh.elem(eID).x(2);
            %c values at element nodes
            c0 = C_FEM(eID); 
            c1 = C_FEM(eID+1); 
            
            %Loop through gauss points for the element
            for i2 = 1:Ngp
                %Get Gaussian weight
                wi = gq.gsw(i2);
                %Evaluate Basis functions at the Gauss Point
                psi0 = EvalBasis(0,gq.xipts(i2),order);
                psi1 = EvalBasis(1,gq.xipts(i2),order);
                %Evaluate x and c
                x = x0*psi0 + x1*psi1;
                C = c0*psi0 + c1*psi1;
                %Get exact solution
                Ce = TransientAnalyticSoln(x,t);
                %Add error at this Gauss point to the cumilating error
                error = error + J*wi*(Ce - C)^2;
            end
        end
        
    case 2  %% Solve for Quadratic Basis Functions
        for eID = 1:Nel
            
            %Get Jacobian for the element
            J = mesh.elem(eID).J;
            %x values element nodes
            x0 = mesh.elem(eID).x(1);
            x1 = mesh.elem(eID).x(2);
            x2 = mesh.elem(eID).x(3);
            %c values at element nodes
            c0 = C_FEM(mesh.elem(eID).n(1)); 
            c1 = C_FEM(mesh.elem(eID).n(2));
            c2 = C_FEM(mesh.elem(eID).n(3));
            
            %% Loop through gauss points for the element
            for i2 = 1:Ngp
                %Get Gaussian weight
                wi = gq.gsw(i2);
                %Evaluate Basis functions at the Gauss Point
                psi0 = EvalBasis(0,gq.xipts(i2),order);
                psi1 = EvalBasis(1,gq.xipts(i2),order);
                psi2 = EvalBasis(2,gq.xipts(i2),order);
                %Evaluate x and c
                x = x0*psi0 + x1*psi1 + x2*psi2;
                C = c0*psi0 + c1*psi1 + c2*psi2;
                %Get exact solution
                Ce = TransientAnalyticSoln(x,t);
                
                %Add error at this Gauss point to the cumilating error
                error = error + J*wi*(Ce - C)^2;
            end
        end
end

%% Square root the error to get the L2 norm
L2norm = sqrt(error);

end