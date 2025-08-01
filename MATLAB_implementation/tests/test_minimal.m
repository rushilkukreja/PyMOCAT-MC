%% Minimal MATLAB Test
% Verify MOCAT-MC MATLAB is working

function test_minimal()
    try
        % Add paths
        addpath(genpath('../supporting_data/')); 
        addpath(genpath('../supporting_functions'));
        fprintf('Successfully added paths\n');
        
        % Try to setup config
        seed = 1;
        cfgMC = setup_MCconfig(seed, '2020.mat');
        fprintf('Successfully created config with %d initial objects\n', size(cfgMC.mat_sats, 1));
        
        % Set minimal simulation
        cfgMC.n_time = 2;  % Just 2 time steps
        cfgMC.skipCollisions = 1;
        cfgMC.launch_model = 'no_launch';
        cfgMC.P_frag = 0;
        
        fprintf('Running minimal simulation...\n');
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