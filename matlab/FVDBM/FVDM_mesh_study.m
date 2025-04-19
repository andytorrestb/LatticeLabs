clc;
clear all;

i = 0;
results_file = 'tau_resolution_study_' + string(i) + '.txt';

while isfile(results_file)
    i=i+1;
    results_file = 'tau_resolution_study_' + string(i) + '.txt';
end

diary(results_file)

% Tau = [1, 0.5, 0.25];
Tau = [1];
dx = [11, 21, 41, 81, 161, 321];

for i=1:length(dx)
    for j=1:length(Tau)
        fprintf("Running Study for dx = %d amd Tau = %f\n", dx(i), Tau(j))
        LBM_Heat_Conduction_AT(dx(i), Tau(j))
        
        for k=1:3
            fprintf('=================================================\n')
        end
    end
end