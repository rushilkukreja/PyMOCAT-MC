% Test Full Default scenario with fixed angl_vec function
addpath('../../');
addpath('../../supporting_functions');
addpath('../../supporting_functions/new_analytic_propagator');
addpath('../../supporting_data/TLEhistoric');

fprintf('Testing Full Default scenario with fixed angl_vec function\n');

try
    tic;
    cfgMC = setup_MCconfig(1, '2020.mat');
    % Full default - no modifications needed
    
    initial_pop = size(cfgMC.mat_sats, 1);
    fprintf('Initial population: %d\n', initial_pop);
    fprintf('Time steps: %d\n', cfgMC.n_time);
    
    [nS, nD, nN, nB, ~] = main_mc(cfgMC, 1);
    
    elapsed = toc;
    total = nS + nD + nN + nB;
    
    fprintf('SUCCESS: Full Default completed!\n');
    fprintf('Results: S=%d, D=%d, N=%d, B=%d\n', nS, nD, nN, nB);
    fprintf('Total: %d\n', total);
    fprintf('Time: %.3f seconds\n', elapsed);
    
    % Save result
    result = struct('test_name', 'Full Default', 'seed', 1, 'success', true, ...
                   'initial_pop', initial_pop, 'final_pop', total, ...
                   'nS', nS, 'nD', nD, 'nN', nN, 'nB', nB, ...
                   'elapsed_time', elapsed, 'n_time_steps', cfgMC.n_time);
    
    json_str = jsonencode(result);
    fid = fopen('matlab_full_default_result.json', 'w');
    fprintf(fid, '%s', json_str);
    fclose(fid);
    
catch ME
    fprintf('ERROR: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
end