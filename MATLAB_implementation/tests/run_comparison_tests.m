%% Run MATLAB Comparison Tests
% Tests matching the Python comparison tests

function run_comparison_tests()
    % Change to Quick_Start directory where all functions work
    cd('Examples/Quick_Start');
    
    fprintf('MATLAB MOCAT-MC Comparison Test Suite\n');
    fprintf('%s\n', repmat('=', 1, 50));
    
    % Define tests matching Python exactly
    tests = {
        'Basic Propagation', 42, @basic_propagation_config;
        'Collision Test', 123, @collision_test_config;
        'Atmospheric Drag', 333, @atmospheric_drag_config;
        'Full Default', 1, @full_default_config;
    };
    
    results = [];
    
    for i = 1:size(tests, 1)
        test_name = tests{i, 1};
        seed = tests{i, 2};
        config_func = tests{i, 3};
        
        fprintf('\n%s\n', repmat('=', 1, 50));
        fprintf('Test: %s\n', test_name);
        fprintf('Seed: %d\n', seed);
        fprintf('%s\n', repmat('=', 1, 50));
        
        try
            tic;
            
            % Setup config
            cfgMC = setup_MCconfig(seed, '2020.mat');
            
            % Apply test-specific config
            cfgMC = config_func(cfgMC);
            
            initial_pop = size(cfgMC.mat_sats, 1);
            fprintf('Initial population: %d\n', initial_pop);
            fprintf('Time steps: %d\n', cfgMC.n_time);
            
            % Run simulation
            [nS, nD, nN, nB, mat_sats] = main_mc(cfgMC, seed);
            
            elapsed = toc;
            total = nS + nD + nN + nB;
            ratio = nS / total;
            
            % Store result
            result.test_name = test_name;
            result.seed = seed;
            result.success = true;
            result.initial_pop = initial_pop;
            result.final_pop = total;
            result.nS = nS;
            result.nD = nD;
            result.nN = nN;
            result.nB = nB;
            result.satellite_ratio = ratio;
            result.elapsed_time = elapsed;
            result.n_time_steps = cfgMC.n_time;
            
            fprintf('Results: S=%d, D=%d, N=%d, B=%d\n', nS, nD, nN, nB);
            fprintf('Total: %d, Ratio: %.6f\n', total, ratio);
            fprintf('Time: %.2fs\n', elapsed);
            
        catch ME
            fprintf('ERROR: %s\n', ME.message);
            result.test_name = test_name;
            result.seed = seed;
            result.success = false;
            result.error = ME.message;
        end
        
        results = [results, result];
    end
    
    % Print summary
    fprintf('\n%s\n', repmat('=', 1, 60));
    fprintf('MATLAB TEST RESULTS SUMMARY\n');
    fprintf('%s\n', repmat('=', 1, 60));
    
    fprintf('\n%-20s %-10s %-8s %-8s %-6s %-6s %-6s %-6s %-6s\n', ...
        'Test', 'Status', 'Initial', 'Final', 'S', 'D', 'N', 'B', 'Time');
    fprintf('%s\n', repmat('-', 1, 80));
    
    for i = 1:length(results)
        r = results(i);
        if r.success
            fprintf('%-20s %-10s %-8d %-8d %-6d %-6d %-6d %-6d %-6.1fs\n', ...
                r.test_name, 'PASS', r.initial_pop, r.final_pop, ...
                r.nS, r.nD, r.nN, r.nB, r.elapsed_time);
        else
            fprintf('%-20s %-10s Error: %s\n', r.test_name, 'FAIL', r.error);
        end
    end
    
    % Save results
    save('matlab_comparison_results.mat', 'results');
    fprintf('\nResults saved to matlab_comparison_results.mat\n');
end

function cfgMC = basic_propagation_config(cfgMC)
    % Short propagation test
    cfgMC.n_time = 5;  % Very short
    cfgMC.skipCollisions = 1;
    cfgMC.P_frag = 0;
    cfgMC.launch_model = 'no_launch';
end

function cfgMC = collision_test_config(cfgMC)
    % Collision detection test
    cfgMC.n_time = 10;
    cfgMC.dt_days = 5;
    cfgMC.skipCollisions = 0;
    cfgMC.CUBE_RES = 50;
    cfgMC.launch_model = 'no_launch';
end

function cfgMC = atmospheric_drag_config(cfgMC)
    % Atmospheric drag test
    cfgMC.altitude_limit_low = 200;
    cfgMC.altitude_limit_up = 600;
    cfgMC.dt_days = 2;
    cfgMC.n_time = 10;
    cfgMC.skipCollisions = 1;
end

function cfgMC = full_default_config(cfgMC)
    % Full default test - no modifications
end