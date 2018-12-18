function [skin] = getSkinProperties(bloodflow)
%%%% Setting properties through the skin for the specific case defined
%%%% in part 2 of Coursework 2. 

%% Set Sub-Cutaneous Properties
skin.SubCut.start = 0.005;  %Start point (m)
skin.SubCut.end = 0.01;     %End point (m)
skin.SubCut.k = 20;     %Thermal Conductivity (W / m K)
skin.SubCut.G = 0.0375; %Blood Flow Rate(m^3/s)
skin.SubCut.c = 3300;   %Specific Heat Capacity of tissue (J/kgK)
skin.SubCut.rho = 1200;     %Density of the tissue (kg/m^3)
skin.SubCut.rhoB = 1060;    %Density of blood (kg/m^3)
skin.SubCut.cB = 3770;  %Specific Heat Capacity of blood (J/kgK)
skin.SubCut.TB = 310.15;    %Temperature of blood (K)


%% Set Dermis Properties
skin.dermis.start = 15/9000;  %Start point (m)
skin.dermis.end = 0.005;      %End point (m)
skin.dermis.k = 40;     %Thermal Conductivity (W / m K)
skin.dermis.G = 0.0375; %Blood Flow Rate(m^3/s)
skin.dermis.c = 3300;   %Specific Heat Capacity of tissue (J/kgK)
skin.dermis.rho = 1200;     %Density of the tissue (kg/m^3)
skin.dermis.rhoB = 1060;    %Density of blood (kg/m^3)
skin.dermis.cB = 3770;  %Specific Heat Capacity of blood (J/kgK)
skin.dermis.TB = 310.15;    %Temperature of blood (K)

%% Set epidermis Properties
skin.epidermis.start = 0;       %Start point (m)
skin.epidermis.end = 15/9000;   %End point (m)
skin.epidermis.k = 25;      %Thermal Conductivity (W / m K)
skin.epidermis.G = 0;       %Blood Flow Rate(m^3/s)
skin.epidermis.c = 3300;   %Specific Heat Capacity of tissue (J/kgK)
skin.epidermis.rho = 1200;     %Density of the tissue (kg/m^3)
skin.epidermis.rhoB = 0;    %Density of blood (kg/m^3)
skin.epidermis.cB = 0;  %Specific Heat Capacity of blood (J/kgK)
skin.epidermis.TB = 0;    %Temperature of blood (K)

%% Set all bloodflow to zero for no bloodflow condition
if bloodflow == false
    skin.SubCut.G = 0;
    skin.dermis.G = 0;
end
       