clc;
clear all;

i = 0;
results_file = 'resolution_study_' + string(i) + '.txt';

while isfile(results_file)
    i=i+1;
    results_file = 'resolution_study_' + string(i) + '.txt';
end

diary(results_file)


dx = [100, 200, 300, 400];



for i=1:length(dx)

        fprintf("Running Study for dx = %d\n", dx(i))
        FVDBM(dx(i));
        
        for k=1:3
            fprintf('=================================================\n')
        end

end