function [GammaGrad] = EvalGammaGrad(T)

%Calculates the gradient of Gamma at the point T

GammaGrad = (2*10^98 * exp(-12017))*(- (400 * exp(20/(20*T - 5464))) / (20*T - 5463)^2);