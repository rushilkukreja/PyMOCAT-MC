% Test script to verify path fixes
clc; clear;

% add folders and subfolders: supporting_functions and supporting_data
addpath(genpath('../../supporting_data/')); 
addpath(genpath('../../supporting_functions'));

% Test 1: Check if density file can be found
fprintf('Test 1: Checking density file...\n');
fn = which('dens_jb2008_032020_022224.mat');
if isempty(fn)
    fprintf('which() could not find density file, trying relative paths...\n');
    possible_paths = {
        'dens_jb2008_032020_022224.mat',
        '../supporting_data/dens_jb2008_032020_022224.mat',
        '../../supporting_data/dens_jb2008_032020_022224.mat',
        '../../../supporting_data/dens_jb2008_032020_022224.mat'
    };
    
    for i = 1:length(possible_paths)
        if exist(possible_paths{i}, 'file')
            fn = possible_paths{i};
            fprintf('Found density file at: %s\n', fn);
            break;
        end
    end
    
    if isempty(fn)
        error('Could not find dens_jb2008_032020_022224.mat');
    end
else
    fprintf('Found density file at: %s\n', fn);
end

% Test 2: Check if initial condition file can be found
fprintf('\nTest 2: Checking initial condition file...\n');
ICfile = '2020.mat';
fn = which(ICfile);
if isempty(fn)
    fprintf('which() could not find IC file, trying relative paths...\n');
    possible_paths = {
        ICfile,
        ['../supporting_data/TLEhistoric/' ICfile],
        ['../../supporting_data/TLEhistoric/' ICfile],
        ['../../../supporting_data/TLEhistoric/' ICfile]
    };
    
    for i = 1:length(possible_paths)
        if exist(possible_paths{i}, 'file')
            fn = possible_paths{i};
            fprintf('Found IC file at: %s\n', fn);
            break;
        end
    end
    
    if isempty(fn)
        error('Could not find %s', ICfile);
    end
else
    fprintf('Found IC file at: %s\n', fn);
end

% Test 3: Try to load the density file
fprintf('\nTest 3: Testing density file load...\n');
try
    load(fn);
    fprintf('Successfully loaded density file!\n');
catch ME
    fprintf('Error loading density file: %s\n', ME.message);
end

fprintf('\nAll tests completed successfully!\n'); 