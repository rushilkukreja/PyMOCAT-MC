%% Run all MATLAB benchmarks with error handling
% Add paths
addpath('../../');
addpath('../../supporting_functions');
addpath('../../supporting_functions/new_analytic_propagator');
addpath('../../supporting_data/TLEhistoric');

fprintf('MATLAB Complete Benchmark Suite\n');
fprintf('===============================\n');

% Test configurations
tests = {
    'Basic Propagation', 42, 5, 1, 0, 'no_launch';
    'Collision Test', 123, 10, 0, 50, 'no_launch';
    'Atmospheric Drag', 333, 10, 1, 0, 'no_launch';
    'Full Default', 1, 74, 0, 0, '';
};

results = {};

for i = 1:size(tests, 1)
    test_name = tests{i, 1};
    seed = tests{i, 2};
    n_time = tests{i, 3};
    skip_collisions = tests{i, 4};
    cube_res = tests{i, 5};
    launch_model = tests{i, 6};
    
    fprintf('\n=== %s ===\n', test_name);
    
    try
        tic;
        
        % Setup config
        cfgMC = setup_MCconfig(seed, '2020.mat');
        cfgMC.n_time = n_time;
        cfgMC.skipCollisions = skip_collisions;
        
        if cube_res > 0
            cfgMC.CUBE_RES = cube_res;
        end
        
        if ~isempty(launch_model)
            cfgMC.launch_model = launch_model;
        end
        
        if strcmp(test_name, 'Atmospheric Drag')
            cfgMC.altitude_limit_low = 200;
            cfgMC.altitude_limit_up = 600;
            cfgMC.dt_days = 2;
        elseif strcmp(test_name, 'Collision Test')
            cfgMC.dt_days = 5;
        end
        
        initial_pop = size(cfgMC.mat_sats, 1);
        fprintf('Running with %d objects, %d time steps...\n', initial_pop, n_time);
        
        % Run simulation
        [nS, nD, nN, nB, ~] = main_mc(cfgMC, seed);
        
        elapsed = toc;
        total = nS + nD + nN + nB;
        ratio = nS / total;
        
        % Store result
        result = struct( ...
            'test_name', test_name, ...
            'seed', seed, ...
            'success', true, ...
            'initial_pop', initial_pop, ...
            'final_pop', total, ...
            'nS', nS, ...
            'nD', nD, ...
            'nN', nN, ...
            'nB', nB, ...
            'satellite_ratio', ratio, ...
            'elapsed_time', elapsed, ...
            'n_time_steps', n_time ...
        );
        
        results{end+1} = result;
        
        fprintf('SUCCESS: Time=%.3fs, Objects=%d\n', elapsed, total);
        
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        
        % Store failed result
        result = struct( ...
            'test_name', test_name, ...
            'seed', seed, ...
            'success', false, ...
            'error', ME.message, ...
            'elapsed_time', 0 ...
        );
        
        results{end+1} = result;
    end
end

% Save all results
fprintf('\n=== SUMMARY ===\n');
json_str = jsonencode(results);
fid = fopen('matlab_all_benchmarks.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);

for i = 1:length(results)
    r = results{i};
    if r.success
        fprintf('%s: %.3f seconds\n', r.test_name, r.elapsed_time);
    else
        fprintf('%s: FAILED\n', r.test_name);
    end
end

fprintf('\nResults saved to matlab_all_benchmarks.json\n');