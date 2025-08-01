%% Simple MATLAB Benchmark - Run one test to get timing data
% Just run Basic Propagation test to get MATLAB timing

% Add necessary paths
addpath('../../');
addpath('../../supporting_functions');
addpath('../../supporting_functions/new_analytic_propagator');
addpath('../../supporting_data/TLEhistoric');

fprintf('MATLAB Basic Propagation Benchmark\n');
fprintf('===================================\n');

try
    tic;
    
    % Setup config for Basic Propagation (matching Python exactly)
    cfgMC = setup_MCconfig(42, '2020.mat');
    cfgMC.n_time = 5;  % Very short
    cfgMC.skipCollisions = 1;
    cfgMC.P_frag = 0;
    cfgMC.launch_model = 'no_launch';
    
    initial_pop = size(cfgMC.mat_sats, 1);
    fprintf('Initial population: %d\n', initial_pop);
    fprintf('Time steps: %d\n', cfgMC.n_time);
    
    % Run simulation
    [nS, nD, nN, nB, mat_sats] = main_mc(cfgMC, 42);
    
    elapsed = toc;
    total = nS + nD + nN + nB;
    ratio = nS / total;
    
    fprintf('\nResults:\n');
    fprintf('S=%d, D=%d, N=%d, B=%d\n', nS, nD, nN, nB);
    fprintf('Total: %d\n', total);
    fprintf('Satellite ratio: %.6f\n', ratio);
    fprintf('Elapsed time: %.4f seconds\n', elapsed);
    
    % Save just this result
    result = struct( ...
        'test_name', 'Basic Propagation', ...
        'seed', 42, ...
        'success', true, ...
        'initial_pop', initial_pop, ...
        'final_pop', total, ...
        'nS', nS, ...
        'nD', nD, ...
        'nN', nN, ...
        'nB', nB, ...
        'satellite_ratio', ratio, ...
        'elapsed_time', elapsed, ...
        'n_time_steps', cfgMC.n_time ...
    );
    
    % Export to JSON for Python
    json_str = jsonencode(result);
    fid = fopen('matlab_basic_benchmark.json', 'w');
    fprintf(fid, '%s', json_str);
    fclose(fid);
    fprintf('\nResult saved to matlab_basic_benchmark.json\n');
    
catch ME
    fprintf('ERROR: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
end