%% Direct comparison tests in Quick_Start directory
function test_comparison_direct()
    fprintf('MATLAB Direct Comparison Tests\n');
    fprintf('%s\n', repmat('=', 1, 50));
    
    % Test 1: Basic Propagation (seed 42, 5 steps)
    fprintf('\nTest 1: Basic Propagation (Seed 42)\n');
    cfgMC = setup_MCconfig(42, '2020.mat');
    cfgMC.n_time = 5;
    cfgMC.skipCollisions = 1;
    cfgMC.P_frag = 0;
    cfgMC.launch_model = 'no_launch';
    
    tic;
    [nS, nD, nN, nB, ~] = main_mc(cfgMC, 42);
    t1 = toc;
    total1 = nS + nD + nN + nB;
    fprintf('Results: S=%d, D=%d, N=%d, B=%d (Total: %d) Time: %.2fs\n', nS, nD, nN, nB, total1, t1);
    
    % Test 2: Collision Test (seed 123, 10 steps)  
    fprintf('\nTest 2: Collision Test (Seed 123)\n');
    cfgMC = setup_MCconfig(123, '2020.mat');
    cfgMC.n_time = 10;
    cfgMC.dt_days = 5;
    cfgMC.skipCollisions = 0;
    cfgMC.CUBE_RES = 50;
    cfgMC.launch_model = 'no_launch';
    
    tic;
    [nS, nD, nN, nB, ~] = main_mc(cfgMC, 123);
    t2 = toc;
    total2 = nS + nD + nN + nB;
    fprintf('Results: S=%d, D=%d, N=%d, B=%d (Total: %d) Time: %.2fs\n', nS, nD, nN, nB, total2, t2);
    
    % Test 3: Atmospheric Drag (seed 333, 10 steps)
    fprintf('\nTest 3: Atmospheric Drag (Seed 333)\n');
    cfgMC = setup_MCconfig(333, '2020.mat');
    cfgMC.altitude_limit_low = 200;
    cfgMC.altitude_limit_up = 600;
    cfgMC.dt_days = 2;
    cfgMC.n_time = 10;
    cfgMC.skipCollisions = 1;
    
    tic;
    [nS, nD, nN, nB, ~] = main_mc(cfgMC, 333);
    t3 = toc;
    total3 = nS + nD + nN + nB;
    fprintf('Results: S=%d, D=%d, N=%d, B=%d (Total: %d) Time: %.2fs\n', nS, nD, nN, nB, total3, t3);
    
    % Test 4: Full Default (seed 1, full year)  
    fprintf('\nTest 4: Full Default (Seed 1)\n');
    cfgMC = setup_MCconfig(1, '2020.mat');
    % No modifications - use defaults
    
    tic;
    [nS, nD, nN, nB, ~] = main_mc(cfgMC, 1);
    t4 = toc;
    total4 = nS + nD + nN + nB;
    fprintf('Results: S=%d, D=%d, N=%d, B=%d (Total: %d) Time: %.2fs\n', nS, nD, nN, nB, total4, t4);
    
    % Summary
    fprintf('\n%s\n', repmat('=', 1, 60));
    fprintf('MATLAB COMPARISON RESULTS SUMMARY\n');
    fprintf('%s\n', repmat('=', 1, 60));
    fprintf('Test 1 (Basic):      S=%4d, D=%4d, N=%4d, B=%4d, Total=%5d\n', ...
        results1.nS, results1.nD, results1.nN, results1.nB, total1);
    fprintf('Test 2 (Collision):  Results will be shown above\n');
    fprintf('Test 3 (Atm Drag):   Results will be shown above\n');  
    fprintf('Test 4 (Full):       Results will be shown above\n');
    
    % Compare with Python results
    fprintf('\nPython Results (for comparison):\n');
    fprintf('Test 1: S=1421, D=1843, N=9421, B=1024, Total=13709\n');
    fprintf('Test 2: S=1421, D=1825, N=9407, B=1019, Total=13672\n');
    fprintf('Test 3: S=1421, D=1826, N=9407, B=1019, Total=13673\n');
    fprintf('Test 4: S=1382, D=1805, N=9363, B=1008, Total=13558\n');
end