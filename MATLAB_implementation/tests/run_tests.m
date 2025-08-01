%% Run MOCAT-MC MATLAB Tests
% Main script to run all tests

fprintf('\nMOCAT-MC MATLAB Test Runner\n');
fprintf('Select test to run:\n');
fprintf('1. Minimal test (verify setup)\n');
fprintf('2. Quick test scenarios\n');
fprintf('3. Full test suite (all scenarios)\n');
fprintf('4. Exit\n');

choice = input('Enter choice (1-4): ');

switch choice
    case 1
        fprintf('\nRunning minimal test...\n');
        test_minimal();
        
    case 2
        fprintf('\nRunning quick test scenarios...\n');
        test_quick_scenarios();
        
    case 3
        fprintf('\nRunning full test suite (this may take a while)...\n');
        test_all_scenarios();
        
    case 4
        fprintf('Exiting...\n');
        
    otherwise
        fprintf('Invalid choice\n');
end