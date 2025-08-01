%% Quick Test Scenarios for MOCAT-MC MATLAB
% Reduced test suite for faster execution

function test_quick_scenarios()
    % Add paths
    addpath(genpath('../supporting_data/')); 
    addpath(genpath('../supporting_functions'));
    
    fprintf('MOCAT-MC MATLAB Quick Test Suite\n');
    fprintf('%s\n', repmat('=', 1, 50));
    
    % Quick tests with reduced time steps
    tests = {
        'Basic Propagation',  42;
        'Collision Test',     123;
        'Launch Test',        200;
        'Fragmentation Test', 999;
        'Full Integration',   1;
    };
    
    results = [];
    
    for i = 1:size(tests, 1)
        test_name = tests{i, 1};
        seed = tests{i, 2};
        
        fprintf('\n%s\n', repmat('=', 1, 50));
        fprintf('Test: %s\n', test_name);
        fprintf('Seed: %d\n', seed);
        fprintf('%s\n', repmat('=', 1, 50));
        
        try
            tic;
            
            % Setup base configuration
            cfgMC = setup_MCconfig(seed, '2020.mat');
            
            % Apply test-specific settings (reduced for speed)
            switch test_name
                case 'Basic Propagation'
                    cfgMC.n_time = 5;  % Very short
                    cfgMC.skipCollisions = 1;
                    cfgMC.P_frag = 0;
                    cfgMC.launch_model = 'no_launch';
                    
                case 'Collision Test'
                    cfgMC.n_time = 10;
                    cfgMC.dt_days = 5;
                    cfgMC.skipCollisions = 0;
                    cfgMC.CUBE_RES = 50;
                    cfgMC.launch_model = 'no_launch';
                    
                case 'Launch Test'
                    cfgMC.n_time = 20;
                    cfgMC.launch_model = 'matsat';
                    cfgMC.launchRepeatYrs = [2018, 2020];
                    cfgMC.skipCollisions = 1;
                    
                case 'Fragmentation Test'
                    cfgMC.n_time = 10;
                    cfgMC.P_frag = 1e-4;  % Higher for testing
                    cfgMC.max_frag = 50;
                    cfgMC.skipCollisions = 1;
                    
                case 'Full Integration'
                    cfgMC.n_time = 15;
                    % All defaults
            end
            
            % Get initial stats
            initial_pop = size(cfgMC.mat_sats, 1);
            
            % Run simulation
            [nS, nD, nN, nB, mat_sats] = main_mc(cfgMC, seed);
            
            % Calculate results
            elapsed_time = toc;
            total_objects = nS + nD + nN + nB;
            
            % Store result
            result.test = test_name;
            result.success = true;
            result.initial = initial_pop;
            result.final = total_objects;
            result.nS = nS;
            result.nD = nD;
            result.nN = nN;
            result.nB = nB;
            result.time = elapsed_time;
            
            fprintf('Initial Population: %d\n', initial_pop);
            fprintf('Final Population: %d\n', total_objects);
            fprintf('Change: %+d\n', total_objects - initial_pop);
            fprintf('S=%d, D=%d, N=%d, B=%d\n', nS, nD, nN, nB);
            fprintf('Time: %.2fs\n', elapsed_time);
            
        catch ME
            fprintf('ERROR: %s\n', ME.message);
            result.test = test_name;
            result.success = false;
            result.error = ME.message;
        end
        
        results = [results, result];
    end
    
    % Print summary
    fprintf('\n%s\n', repmat('=', 1, 50));
    fprintf('SUMMARY\n');
    fprintf('%s\n', repmat('=', 1, 50));
    
    fprintf('\n%-20s %-10s %-8s %-8s %-8s\n', 'Test', 'Status', 'Initial', 'Final', 'Time');
    fprintf('%s\n', repmat('-', 1, 60));
    
    for i = 1:length(results)
        r = results(i);
        if r.success
            fprintf('%-20s %-10s %-8d %-8d %-8.2fs\n', ...
                r.test, 'PASS', r.initial, r.final, r.time);
        else
            fprintf('%-20s %-10s Error: %s\n', r.test, 'FAIL', r.error);
        end
    end
    
    % Save results
    save('quick_test_results.mat', 'results');
    fprintf('\nResults saved to quick_test_results.mat\n');
end