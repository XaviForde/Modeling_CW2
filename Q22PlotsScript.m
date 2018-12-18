%% Script to plot the temperature profiles found at iterations
%   of the FindBC.m algorithm to show how the solution was 
%   converged upon

TempMatrix = csvread('Q22TempVec.csv');

TimeVec = linspace(0, 50, 1001);
%TimeVec = TimeVec(1:200);

%% Plot the results

figure()
plot([TimeVec(1) TimeVec(end)], [317.15, 317.15], '--k') %Plot the Tburn line
text(20, 319, '$T_{BURN}$', 'interpreter', 'latex', 'FontSize', 12)
hold on
plot(TimeVec, TempMatrix(3,:), '-' ) %Plot the low temp initial guess
hold on 
plot(TimeVec, TempMatrix(5,:), '-') %Plot the high temp initial guess
hold on
plot(TimeVec, TempMatrix(7,:), '-r') %Plot the first iteration
hold on
% plot(TimeVec, TempMatrix(8,:)) %Plot the third iteration
% hold on
plot(TimeVec, TempMatrix(end,:), '-b', 'LineWidth', 0.8) %Plot the final result

%% Format the figure
title('Temperature Profile As Algorithm Iterates', 'interpreter' ,'latex', 'FontSize', 14)
lgd = legend({'$T_{BURN}$', 'Iteration 1', 'Iteration 3','Iteration 5', 'Final Result'},'Location', 'northeast', 'interpreter', 'latex');
xlabel('Time \ in \ seconds','interpreter','latex', 'FontSize', 12);
ylabel('Temperature in Kelvin', 'interpreter','latex', 'FontSize', 12);
grid on
%Set axis limits
ylim([310 365]);
xlim([0 50]);
%% Saving The Figure
%print -depsc \Users\xav_m\OneDrive\Documents\XAVI\University\Final_Year\Systems_Mod\Modeling_CW2\Report\Figures\epsQ22BCTempProfs
