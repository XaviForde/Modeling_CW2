function [DiffLEM] = DiffLEM(msh, Ngp, eID, order)
%% Returns local diffusion element matrix
%%%%Could take FieldVal out of the loop and multiply at end??
gq = CreateGQScheme(Ngp);
Ngp = gq.npts;
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