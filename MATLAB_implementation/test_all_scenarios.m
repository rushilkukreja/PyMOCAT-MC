%% MOCAT-MC MATLAB Test Suite
% Comprehensive test scenarios matching Python implementation
% Run all tests and save results for comparison

function test_all_scenarios()
    % Add paths
    addpath(genpath('../supporting_data/')); 
    addpath(genpath('../supporting_functions'));
    
    fprintf('\n%s\n', repmat('=', 1, 80));
    fprintf('MOCAT-MC MATLAB IMPLEMENTATION TEST SUITE\n');
    fprintf('Started at: %s\n', datestr(now));
    fprintf('%s\n', repmat('=', 1, 80));
    
    % Initialize results structure
    results = {};
    
    % Define all test scenarios matching Python
    test_scenarios = {
        % Test Name,              Seed, Config Function
        'Basic Propagation',      42,   @basic_propagation_config;
        'Collision Detection',    123,  @collision_detection_config;
        'Fragmentation Test',     999,  @fragmentation_config;
        'Orbit Control',          777,  @orbit_control_config;
        'No Launch Baseline',     100,  @no_launch_config;
        'MatSat Launch',          200,  @matsat_launch_config;
        'MatSat Smoothed',        300,  @matsat_smooth_config;
        'Atmospheric Drag',       333,  @atmospheric_drag_config;
        'Extreme Altitudes',      1100, @extreme_altitude_config;
        'High Activity Period',   500,  @high_activity_config;
        'Mixed Launch/Collision', 600,  @mixed_scenario_config;
        'Full Integration',       1,    @full_integration_config;
    };
    
    % Run each test
    for i = 1:size(test_scenarios, 1)
        test_name = test_scenarios{i, 1};
        seed = test_scenarios{i, 2};
        config_func = test_scenarios{i, 3};
        
        result = run_test(test_name, seed, config_func);
        results{end+1} = result;
    end
    
    % Print summary
    print_summary(results);
    
    % Save results for Python comparison
    save_results(results);
end

%% Individual test runner
function result = run_test(test_name, seed, config_func)
    fprintf('\n%s\n', repmat('=', 1, 60));
    fprintf('Running Test: %s\n', test_name);
    fprintf('Seed: %d\n', seed);
    fprintf('%s\n', repmat('=', 1, 60));
    
    try
        % Start timer
        tic;
        
        % Setup configuration
        cfg_mc = setup_MCconfig(seed, '2020.mat');
        
        % Apply test-specific modifications
        cfg_mc = config_func(cfg_mc);
        
        % Get initial stats
        initial_pop = size(cfg_mc.mat_sats, 1);
        n_time = cfg_mc.n_time;
        dt_days = cfg_mc.dt_days;
        
        % Run simulation
        [nS, nD, nN, nB, mat_sats] = main_mc(cfg_mc, seed);
        
        % Calculate results
        elapsed_time = toc;
        total_objects = nS + nD + nN + nB;
        
        % Store results
        result = struct();
        result.test_name = test_name;
        result.seed = seed;
        result.success = true;
        result.initial_pop = initial_pop;
        result.final_pop = total_objects;
        result.nS = nS;
        result.nD = nD;
        result.nN = nN;
        result.nB = nB;
        result.n_time_steps = n_time;
        result.dt_days = dt_days;
        result.elapsed_time = elapsed_time;
        result.satellite_ratio = nS / total_objects;
        result.population_change = total_objects - initial_pop;
        result.error = '';
        
        % Sample orbital elements for comparison (first 10 objects)
        n_sample = min(10, size(mat_sats, 1));
        result.sample_oe = mat_sats(1:n_sample, 1:6);
        result.sample_rv = mat_sats(1:n_sample, 17:22);
        
        fprintf('Initial Population: %d\n', initial_pop);
        fprintf('Final Population: %d (change: %+d)\n', total_objects, total_objects - initial_pop);
        fprintf('Final Counts - S: %d, D: %d, N: %d, B: %d\n', nS, nD, nN, nB);
        fprintf('Satellite Ratio: %.4f\n', result.satellite_ratio);
        fprintf('Execution Time: %.2f seconds\n', elapsed_time);
        
    catch ME
        fprintf('ERROR in test %s: %s\n', test_name, ME.message);
        result = struct();
        result.test_name = test_name;
        result.seed = seed;
        result.success = false;
        result.error = ME.message;
    end
end

%% Test configuration functions (matching Python exactly)

function cfg = basic_propagation_config(cfg)
    % Test 1: Basic propagation without collisions
    cfg.n_time = 10;
    cfg.skipCollisions = 1;
    cfg.launch_model = 'no_launch';
    cfg.P_frag = 0;
end

function cfg = collision_detection_config(cfg)
    % Test 2: Collision detection
    cfg.dt_days = 1;
    cfg.DeltaT = cfg.dt_days * cfg.DAY2MIN;
    cfg.tsince = 0:cfg.DeltaT:(36*cfg.DeltaT);  % ~6 months
    cfg.n_time = length(cfg.tsince);
    cfg.CUBE_RES = 100;
    cfg.alph = 0.01;
    cfg.alph_a = 0;
    cfg.skipCollisions = 0;
end

function cfg = fragmentation_config(cfg)
    % Test 3: Fragmentation/explosion test
    cfg.P_frag = 1e-5;
    cfg.P_frag_cutoff = 20;
    cfg.max_frag = 100;
    cfg.n_time = 36;  % ~6 months
end

function cfg = orbit_control_config(cfg)
    % Test 4: Orbit control test
    cfg.orbtol = 10;
    cfg.step_control = 5;
    cfg.PMD = 0.90;
    cfg.missionlifetime = 5;
end

function cfg = no_launch_config(cfg)
    % Test 5: No launch baseline
    cfg.launch_model = 'no_launch';
    cfg.n_time = 73;  % 1 year
end

function cfg = matsat_launch_config(cfg)
    % Test 6: MatSat launch model
    cfg.launch_model = 'matsat';
    cfg.launchRepeatYrs = [2018, 2022];
    cfg.launchRepeatSmooth = 0;
    cfg.n_time = 146;  % 2 years
end

function cfg = matsat_smooth_config(cfg)
    % Test 7: MatSat with smoothing
    cfg.launch_model = 'matsat';
    cfg.launchRepeatYrs = [2015, 2020];
    cfg.launchRepeatSmooth = 1;
    cfg.n_time = 219;  % 3 years
end

function cfg = atmospheric_drag_config(cfg)
    % Test 8: Atmospheric drag focus
    cfg.altitude_limit_low = 200;
    cfg.altitude_limit_up = 600;
    cfg.dt_days = 0.5;
    cfg.DeltaT = cfg.dt_days * cfg.DAY2MIN;
    cfg.tsince = 0:cfg.DeltaT:(182*cfg.DeltaT);  % 3 months
    cfg.n_time = length(cfg.tsince);
end

function cfg = extreme_altitude_config(cfg)
    % Test 10: Extreme altitude test
    cfg.altitude_limit_low = 150;
    cfg.altitude_limit_up = 50000;
    cfg.collision_alt_limit = 50000;
    cfg.n_time = 73;  % 1 year
end

function cfg = high_activity_config(cfg)
    % Test 9: High activity period
    cfg.launch_model = 'matsat';
    cfg.launchRepeatYrs = [2019, 2022];
    cfg.launchRepeatSmooth = 0;
    cfg.missionlifetime = 5;
    cfg.PMD = 0.90;
    cfg.alph = 0.02;
    cfg.n_time = 365;  % 5 years
end

function cfg = mixed_scenario_config(cfg)
    % Test 11: Mixed launch and collision
    cfg.launch_model = 'matsat';
    cfg.launchRepeatYrs = [2018, 2020];
    cfg.skipCollisions = 0;
    cfg.alph = 0.05;
    cfg.CUBE_RES = 100;
    cfg.n_time = 146;  % 2 years
end

function cfg = full_integration_config(cfg)
    % Test 12: Full integration with defaults
    % Use all defaults - 1 year simulation
end

%% Summary printer
function print_summary(results)
    fprintf('\n%s\n', repmat('=', 1, 80));
    fprintf('TEST RESULTS SUMMARY\n');
    fprintf('%s\n', repmat('=', 1, 80));
    
    successful = sum(cellfun(@(x) x.success, results));
    total = length(results);
    
    fprintf('\nTotal Tests: %d\n', total);
    fprintf('Successful: %d\n', successful);
    fprintf('Failed: %d\n', total - successful);
    
    fprintf('\nDetailed Results:\n');
    fprintf('%-25s %-10s %-8s %-8s %-6s %-6s %-6s %-6s %-8s\n', ...
        'Test Name', 'Status', 'Initial', 'Final', 'S', 'D', 'N', 'B', 'Time(s)');
    fprintf('%s\n', repmat('-', 1, 100));
    
    for i = 1:length(results)
        r = results{i};
        if r.success
            fprintf('%-25s %-10s %-8d %-8d %-6d %-6d %-6d %-6d %-8.2f\n', ...
                r.test_name, 'PASS', r.initial_pop, r.final_pop, ...
                r.nS, r.nD, r.nN, r.nB, r.elapsed_time);
        else
            fprintf('%-25s %-10s Error: %s\n', r.test_name, 'FAIL', r.error);
        end
    end
end

%% Save results for comparison
function save_results(results)
    % Save as MAT file for easy Python loading
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('matlab_test_results_%s.mat', timestamp);
    
    % Convert cell array to structure array for better compatibility
    results_struct = [];
    for i = 1:length(results)
        results_struct = [results_struct, results{i}];
    end
    
    save(filename, 'results_struct', '-v7');
    fprintf('\nResults saved to: %s\n', filename);
    
    % Also save as CSV for easy viewing
    csv_filename = sprintf('matlab_test_summary_%s.csv', timestamp);
    fid = fopen(csv_filename, 'w');
    fprintf(fid, 'Test_Name,Status,Seed,Initial_Pop,Final_Pop,nS,nD,nN,nB,Time_s\n');
    
    for i = 1:length(results)
        r = results{i};
        if r.success
            fprintf(fid, '%s,PASS,%d,%d,%d,%d,%d,%d,%d,%.2f\n', ...
                r.test_name, r.seed, r.initial_pop, r.final_pop, ...
                r.nS, r.nD, r.nN, r.nB, r.elapsed_time);
        else
            fprintf(fid, '%s,FAIL,%d,,,,,,,\n', r.test_name, r.seed);
        end
    end
    fclose(fid);
    fprintf('Summary saved to: %s\n', csv_filename);
end