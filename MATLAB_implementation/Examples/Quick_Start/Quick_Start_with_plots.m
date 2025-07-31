% Quick Start with Comprehensive Plots
clc;clear;

% add folders and subfolders: supporting_functions and supporting_data
addpath(genpath('../../supporting_data/')); 
addpath(genpath('../../supporting_functions'));

% initial condition file
ICfile = '2020.mat'

% MOCAT MC configuration
seed = 1;% random number generator seed

disp('MC configuration starting...');
cfgMC = setup_MCconfig(seed,ICfile);
fprintf('Seed %i\n', seed);

% MOCAT MC evolution
fprintf('Initial Population:  %i sats\n', size(cfgMC.mat_sats,1));
fprintf('Launches per year: %i\n', size(cfgMC.repeatLaunches,1));
disp('Starting main_mc...');
[nS,nD,nN,nB,mat_sats]=main_mc(cfgMC,seed);

% MOCAT MC postprocess: ratio of satellite (SR) among all space objects
ratio = nS/(nS+nD+nN+nB);
fprintf('Quick Start under no launch scenario done!\n')
fprintf('Satellite ratio in all space objects after evolution: %f\n', ratio)

%% Create Comprehensive Plots
disp('Creating plots...');

% Define indices for easier access
idx_a = 1; idx_ecco = 2; idx_inclo = 3; idx_nodeo = 4; idx_argpo = 5; 
idx_mo = 6; idx_bstar = 7; idx_mass = 8; idx_radius = 9; idx_controlled = 11;
idx_objectclass = 23; idx_r = [17 18 19]; idx_v = [20 21 22];

% Figure 1: Population Summary
figure(1);
clf;
set(gcf, 'Color', 'white');
subplot(2,2,1);
bar([nS, nD, nN, nB]);
set(gca, 'XTickLabel', {'Satellites', 'Derelicts', 'Debris', 'Rocket Bodies'});
title('Final Population Distribution');
ylabel('Number of Objects');
grid on;

subplot(2,2,2);
pie([nS, nD, nN, nB], {'Satellites', 'Derelicts', 'Debris', 'Rocket Bodies'});
title('Population Distribution (Pie Chart)');

subplot(2,2,3);
altitudes = (mat_sats(:,idx_a) - 1) * 6378.137; % Convert to km
histogram(altitudes, 20);
xlabel('Altitude (km)');
ylabel('Number of Objects');
title('Altitude Distribution');
grid on;

subplot(2,2,4);
inclinations = rad2deg(mat_sats(:,idx_inclo));
histogram(inclinations, 20);
xlabel('Inclination (degrees)');
ylabel('Number of Objects');
title('Inclination Distribution');
grid on;

% Save Figure 1
print('matlab_quick_start_figure_1_population_summary', '-dpng', '-r150');
fprintf('Saved: matlab_quick_start_figure_1_population_summary.png\n');

% Figure 2: Orbital Elements Analysis
figure(2);
clf;
set(gcf, 'Color', 'white');
subplot(1,3,1);
plot(mat_sats(:,idx_ecco), altitudes, '.');
xlabel('Eccentricity');
ylabel('Altitude (km)');
title('Eccentricity vs Altitude');
grid on;

subplot(1,3,2);
plot(inclinations, altitudes, '.');
xlabel('Inclination (degrees)');
ylabel('Altitude (km)');
title('Inclination vs Altitude');
grid on;

subplot(1,3,3);
histogram(mat_sats(:,idx_objectclass), 1:11);
xlabel('Object Class');
ylabel('Number of Objects');
title('Object Class Distribution');
grid on;

% Save Figure 2
print('matlab_quick_start_figure_2_orbital_elements', '-dpng', '-r150');
fprintf('Saved: matlab_quick_start_figure_2_orbital_elements.png\n');

% Figure 3: 3D Position Plot
figure(3);
clf;
set(gcf, 'Color', 'white');
r = mat_sats(:,idx_r);
scatter3(r(:,1), r(:,2), r(:,3), 20, altitudes, 'filled');
colorbar;
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('3D Object Positions (colored by altitude)');
grid on;
axis equal;

% Save Figure 3
print('matlab_quick_start_figure_3_3d_positions', '-dpng', '-r150');
fprintf('Saved: matlab_quick_start_figure_3_3d_positions.png\n');

% Figure 4: Altitude vs Eccentricity with Object Types
figure(4);
clf;
set(gcf, 'Color', 'white');
hold on;
colors = {'b', 'r', 'g', 'm'};
labels = {'Satellites', 'Derelicts', 'Debris', 'Rocket Bodies'};

% Find indices for each object type
sat_idx = mat_sats(:,idx_controlled) == 1;
derelict_idx = (mat_sats(:,idx_controlled) == 0) & (mat_sats(:,idx_objectclass) == 1);
debris_idx = mat_sats(:,idx_objectclass) == 3;
rb_idx = mat_sats(:,idx_objectclass) == 2;

indices = {sat_idx, derelict_idx, debris_idx, rb_idx};

for i = 1:4
    if any(indices{i})
        plot(mat_sats(indices{i}, idx_ecco), altitudes(indices{i}), '.', ...
             'Color', colors{i}, 'DisplayName', labels{i});
    end
end

xlabel('Eccentricity');
ylabel('Altitude (km)');
title('Altitude vs Eccentricity by Object Type');
legend('show');
grid on;

% Save Figure 4
print('matlab_quick_start_figure_4_altitude_eccentricity', '-dpng', '-r150');
fprintf('Saved: matlab_quick_start_figure_4_altitude_eccentricity.png\n');

% Figure 5: Summary Statistics
figure(5);
clf;
set(gcf, 'Color', 'white');
subplot(2,2,1);
stats_data = [mean(altitudes), std(altitudes), min(altitudes), max(altitudes)];
bar(stats_data);
set(gca, 'XTickLabel', {'Mean', 'Std', 'Min', 'Max'});
title('Altitude Statistics (km)');
grid on;

subplot(2,2,2);
stats_ecc = [mean(mat_sats(:,idx_ecco)), std(mat_sats(:,idx_ecco)), ...
             min(mat_sats(:,idx_ecco)), max(mat_sats(:,idx_ecco))];
bar(stats_ecc);
set(gca, 'XTickLabel', {'Mean', 'Std', 'Min', 'Max'});
title('Eccentricity Statistics');
grid on;

subplot(2,2,3);
stats_mass = [mean(mat_sats(:,idx_mass)), std(mat_sats(:,idx_mass)), ...
              min(mat_sats(:,idx_mass)), max(mat_sats(:,idx_mass))];
bar(stats_mass);
set(gca, 'XTickLabel', {'Mean', 'Std', 'Min', 'Max'});
title('Mass Statistics (kg)');
grid on;

subplot(2,2,4);
stats_radius = [mean(mat_sats(:,idx_radius)), std(mat_sats(:,idx_radius)), ...
                min(mat_sats(:,idx_radius)), max(mat_sats(:,idx_radius))];
bar(stats_radius);
set(gca, 'XTickLabel', {'Mean', 'Std', 'Min', 'Max'});
title('Radius Statistics (m)');
grid on;

% Save Figure 5
print('matlab_quick_start_figure_5_summary_statistics', '-dpng', '-r150');
fprintf('Saved: matlab_quick_start_figure_5_summary_statistics.png\n');

% Print summary statistics
fprintf('\n=== SUMMARY STATISTICS ===\n');
fprintf('Total Objects: %d\n', nS+nD+nN+nB);
fprintf('Satellites: %d (%.1f%%)\n', nS, 100*nS/(nS+nD+nN+nB));
fprintf('Derelicts: %d (%.1f%%)\n', nD, 100*nD/(nS+nD+nN+nB));
fprintf('Debris: %d (%.1f%%)\n', nN, 100*nN/(nS+nD+nN+nB));
fprintf('Rocket Bodies: %d (%.1f%%)\n', nB, 100*nB/(nS+nD+nN+nB));
fprintf('Mean Altitude: %.1f km\n', mean(altitudes));
fprintf('Altitude Range: %.1f - %.1f km\n', min(altitudes), max(altitudes));
fprintf('Mean Eccentricity: %.3f\n', mean(mat_sats(:,idx_ecco)));
fprintf('Mean Mass: %.1f kg\n', mean(mat_sats(:,idx_mass)));
fprintf('Mean Radius: %.2f m\n', mean(mat_sats(:,idx_radius)));

disp('All plots created and saved successfully!');