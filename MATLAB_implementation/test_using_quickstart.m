%% Test using existing Quick Start example
function test_using_quickstart()
    try
        % Change to Quick_Start directory where all files exist
        cd('Examples/Quick_Start');
        
        % Test basic functionality - single time step
        seed = 42;
        ICfile = '2020.mat';
        
        fprintf('Setting up config...\n');
        cfgMC = setup_MCconfig(seed, ICfile);
        
        % Make it very short
        cfgMC.n_time = 2;  % Just 2 time steps
        cfgMC.skipCollisions = 1;
        cfgMC.P_frag = 0;
        
        fprintf('Initial population: %d\n', size(cfgMC.mat_sats, 1));
        fprintf('Running simulation...\n');
        
        [nS, nD, nN, nB, mat_sats] = main_mc(cfgMC, seed);
        
        fprintf('SUCCESS! Final counts: S=%d, D=%d, N=%d, B=%d\n', nS, nD, nN, nB);
        
    catch ME
        fprintf('ERROR: %s\n', ME.message);
        fprintf('Stack:\n');
        for i = 1:length(ME.stack)
            fprintf('  In %s at line %d\n', ME.stack(i).name, ME.stack(i).line);
        end
    end
end