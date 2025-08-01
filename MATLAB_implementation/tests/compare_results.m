%% Compare MATLAB and Python MOCAT-MC Results
% Load and compare test results from both implementations

function compare_results(matlab_file, python_file)
    fprintf('\n%s\n', repmat('=', 1, 80));
    fprintf('MOCAT-MC MATLAB vs Python Comparison\n');
    fprintf('%s\n', repmat('=', 1, 80));
    
    % Load MATLAB results
    if nargin < 1
        matlab_file = 'matlab_test_results.mat';
    end
    if nargin < 2
        python_file = '../python_implementation/python_test_results.mat';
    end
    
    % Load results
    fprintf('\nLoading results...\n');
    matlab_data = load(matlab_file);
    python_data = load(python_file);
    
    matlab_results = matlab_data.results_struct;
    python_results = python_data.results;
    
    % Compare each test
    fprintf('\nTest-by-Test Comparison:\n');
    fprintf('%s\n', repmat('-', 1, 80));
    
    for i = 1:length(matlab_results)
        m_res = matlab_results(i);
        
        % Find matching Python test
        p_res = [];
        for j = 1:length(python_results)
            if strcmp(python_results(j).test_name, m_res.test_name) && ...
               python_results(j).seed == m_res.seed
                p_res = python_results(j);
                break;
            end
        end
        
        if isempty(p_res)
            fprintf('\nTest: %s (Seed: %d) - No Python match found\n', ...
                m_res.test_name, m_res.seed);
            continue;
        end
        
        fprintf('\nTest: %s (Seed: %d)\n', m_res.test_name, m_res.seed);
        
        if m_res.success && p_res.success
            % Compare object counts
            fprintf('  Object Counts:\n');
            fprintf('    MATLAB: S=%d, D=%d, N=%d, B=%d (Total: %d)\n', ...
                m_res.nS, m_res.nD, m_res.nN, m_res.nB, m_res.final_pop);
            fprintf('    Python: S=%d, D=%d, N=%d, B=%d (Total: %d)\n', ...
                p_res.nS, p_res.nD, p_res.nN, p_res.nB, p_res.final_pop);
            
            % Calculate differences
            diff_S = abs(m_res.nS - p_res.nS);
            diff_D = abs(m_res.nD - p_res.nD);
            diff_N = abs(m_res.nN - p_res.nN);
            diff_B = abs(m_res.nB - p_res.nB);
            diff_total = abs(m_res.final_pop - p_res.final_pop);
            
            fprintf('    Differences: S=%d, D=%d, N=%d, B=%d (Total: %d)\n', ...
                diff_S, diff_D, diff_N, diff_B, diff_total);
            
            % Percentage differences
            if m_res.final_pop > 0
                pct_diff = 100 * diff_total / m_res.final_pop;
                fprintf('    Total difference: %.2f%%\n', pct_diff);
            end
            
            % Compare satellite ratio
            fprintf('  Satellite Ratio:\n');
            fprintf('    MATLAB: %.4f\n', m_res.satellite_ratio);
            fprintf('    Python: %.4f\n', p_res.satellite_ratio);
            fprintf('    Difference: %.4f\n', abs(m_res.satellite_ratio - p_res.satellite_ratio));
            
            % Compare execution time
            fprintf('  Execution Time:\n');
            fprintf('    MATLAB: %.2f seconds\n', m_res.elapsed_time);
            fprintf('    Python: %.2f seconds\n', p_res.elapsed_time);
            
            % Compare orbital elements (if available)
            if isfield(m_res, 'sample_oe') && isfield(p_res, 'sample_oe')
                compare_orbital_elements(m_res.sample_oe, p_res.sample_oe);
            end
            
        else
            fprintf('  Test failed in one or both implementations\n');
            if ~m_res.success
                fprintf('    MATLAB Error: %s\n', m_res.error);
            end
            if ~p_res.success
                fprintf('    Python Error: %s\n', p_res.error);
            end
        end
    end
    
    % Overall summary
    print_overall_summary(matlab_results, python_results);
end

function compare_orbital_elements(matlab_oe, python_oe)
    fprintf('  Orbital Elements Comparison (first %d objects):\n', size(matlab_oe, 1));
    
    % Calculate mean absolute differences
    if size(matlab_oe) == size(python_oe)
        diff_oe = abs(matlab_oe - python_oe);
        mean_diff = mean(diff_oe, 1);
        max_diff = max(diff_oe, [], 1);
        
        labels = {'a', 'e', 'i', 'Omega', 'omega', 'M'};
        fprintf('    Mean absolute differences:\n');
        for i = 1:6
            fprintf('      %s: mean=%.6e, max=%.6e\n', labels{i}, mean_diff(i), max_diff(i));
        end
    else
        fprintf('    Warning: Different number of objects, cannot compare\n');
    end
end

function print_overall_summary(matlab_results, python_results)
    fprintf('\n%s\n', repmat('=', 1, 80));
    fprintf('OVERALL SUMMARY\n');
    fprintf('%s\n', repmat('=', 1, 80));
    
    % Count successful tests
    m_success = sum([matlab_results.success]);
    p_success = sum([python_results.success]);
    
    fprintf('\nTest Success Rate:\n');
    fprintf('  MATLAB: %d/%d tests passed\n', m_success, length(matlab_results));
    fprintf('  Python: %d/%d tests passed\n', p_success, length(python_results));
    
    % Calculate overall differences for successful tests
    total_diff = 0;
    count = 0;
    
    for i = 1:length(matlab_results)
        if matlab_results(i).success
            % Find matching Python test
            for j = 1:length(python_results)
                if strcmp(python_results(j).test_name, matlab_results(i).test_name) && ...
                   python_results(j).success
                    diff = abs(matlab_results(i).final_pop - python_results(j).final_pop);
                    total_diff = total_diff + diff;
                    count = count + 1;
                    break;
                end
            end
        end
    end
    
    if count > 0
        avg_diff = total_diff / count;
        fprintf('\nAverage population difference: %.2f objects\n', avg_diff);
    end
    
    fprintf('\nConclusion:\n');
    if avg_diff < 10
        fprintf('  Excellent agreement between MATLAB and Python implementations\n');
    elseif avg_diff < 50
        fprintf('  Good agreement with minor differences (likely due to RNG or numerical precision)\n');
    else
        fprintf('  Significant differences detected - further investigation recommended\n');
    end
end