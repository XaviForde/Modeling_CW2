bloodflow = true;

T_low = 316;   % Know result of Gamma < 1

%[Gamma , ~] = EvalGamma(T_low, bloodflow);
T_high = 394;   % Known result of Gamma >> 1
[Gamma , ~] = EvalGamma(T_high, bloodflow);

i = 1;
%Average the two inital guesses to get a third guess
T_av = (T_low+T_high)/2

%% While the difference between the new guess and previous guess is
%   greater than 0.5K continue iterating toward solution og Gamma = 1
while T_high-T_low > 0.5  
    
    [Gamma, ~]  = EvalGamma(T_av, bloodflow)
    
    if Gamma > 1   %
        T_high = T_av;
    else
        T_low = T_av;
    end
    i = i+1
    T_av = (T_low + T_high)/2
    
end
% Find values of Gamma for the solution Boundary Condition for completeness
[Gamma, ~]  = EvalGamma(T_av, bloodflow)

%Save to csv file for post processing data
%dlmwrite('Q22TempVec.csv', TempVec, 'delimiter', ',', 'precision', 12);



