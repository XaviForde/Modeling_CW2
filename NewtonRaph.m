bloodflow = true;

Tnext = 394;
Tcurrent = 0;   % Known result of Gamma >> 1
i = 0
while abs(Tcurrent - Tnext) > 0.5
    
    i = i+1
    
    Tcurrent = Tnext;
    
    [Gamma , ~] = EvalGamma(Tcurrent, bloodflow);
    Gamma = Gamma - 1
    GammaGrad = EvalGammaGrad(Tcurrent)*50;
    i = 1;
    %Calculate the gradient
    Tnext = Tcurrent - Gamma/GammaGrad
    
end
