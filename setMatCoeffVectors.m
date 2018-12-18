function [mesh] = setMatCoeffVectors(mesh, bloodflow, order)
%% This function sets the Material Coefficients for Diffusion, Reaction
%   and source. The skin properties for each layer are found using
%   getSkinProperties.m and then coefficients are set at each node
% Inputs:
%   mesh -  1D mesh containing mesh information (structure)
%   bloodflow - logical to indicate if there is bloodflow
%               (takes the value 'true' or 'false')

%% Initialise the vectors for coefficent

mesh.DCvec = zeros(order*mesh.ne,1);      %Diffusion Coefficient vector
mesh.RCvec = zeros(order*mesh.ne,1);      %Reaction Coefficient vector
mesh.fvec = zeros(order*mesh.ne,1);       %Source Coefficient vector


%% Calculate Coefficients at each layer
%Get structure containing porperties of the skin and each layer
skin = getSkinProperties(bloodflow);

%Sub Curtaneous Coefficients
skin.SubCut.D = skin.SubCut.k / (skin.SubCut.rho * skin.SubCut.c);

skin.SubCut.lambda = -(skin.SubCut.G * skin.SubCut.rhoB * skin.SubCut.cB)/ ...
    (skin.SubCut.rho * skin.SubCut.c);

skin.SubCut.f = (-1) * skin.SubCut.lambda * skin.SubCut.TB;

%Dermis Coefficients
skin.dermis.D = skin.dermis.k / (skin.dermis.rho * skin.dermis.c);

skin.dermis.lambda = -(skin.dermis.G * skin.dermis.rhoB * skin.dermis.cB)/ ...
    (skin.dermis.rho * skin.dermis.c);

skin.dermis.f = (-1) * skin.dermis.lambda * skin.dermis.TB;

%Epiermis Coefficients
skin.epidermis.D = skin.epidermis.k / (skin.epidermis.rho * skin.epidermis.c);

skin.epidermis.lambda = -(skin.epidermis.G * skin.epidermis.rhoB * skin.epidermis.cB)/ ...
    (skin.epidermis.rho * skin.epidermis.c);

skin.epidermis.f = (-1) * skin.epidermis.lambda * skin.epidermis.TB;

%% Now cycle through nodes and set coefficients as appropiate
for i = 1:length(mesh.nvec)
    
    if mesh.nvec(i) >= skin.epidermis.start && mesh.nvec(i) <= skin.epidermis.end
        mesh.DCvec(i) = skin.epidermis.D;
        mesh.RCvec(i) = skin.epidermis.lambda;
        mesh.fvec(i) = skin.epidermis.f;
        
    elseif mesh.nvec(i) > skin.dermis.start && mesh.nvec(i) <= skin.dermis.end
        
        mesh.DCvec(i) = skin.dermis.D;
        mesh.RCvec(i) = skin.dermis.lambda;
        mesh.fvec(i) = skin.dermis.f;
        
    elseif mesh.nvec(i) > skin.SubCut.start && mesh.nvec(i) <= skin.SubCut.end
        mesh.DCvec(i) = skin.SubCut.D;
        mesh.RCvec(i) = skin.SubCut.lambda;
        mesh.fvec(i) = skin.SubCut.f ;
    else
        error('Node x position is not within the skin layers given')
    end
end
end

