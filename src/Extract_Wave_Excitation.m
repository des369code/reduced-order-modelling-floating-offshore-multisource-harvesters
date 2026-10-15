% =========================================================================
% EXTRACT WAVE EXCITATION FROM AQWA DATA (REPLACES K*eta APPROXIMATION)
% =========================================================================
% Reconstructs the total environmental excitation force from the AQWA
% 6-DOF response data by running the equation of motion in reverse:
%
%   F_total(t) = M * ddX_AQWA(t) + C * dX_AQWA(t) + K * X_AQWA(t)
%
% This F_total includes wave + AQWA steady wind. The steady wind is
% small (~20 N at 6 m/s) and negligible compared to the DMST aero
% forces (~1000 N) that replace it in the coupled simulation.
%
% Self-consistency: when fed back through the ODE with the same M,C,K
% and AQWA initial conditions, this reproduces the AQWA response.
% =========================================================================
clear; clc; close all;

%% 1. LOAD SYSTEM MATRICES AND AQWA RESPONSE DATA
if ~isfile('../outputs/sysid/MasterSystemData.mat')
    error('MasterSystemData.mat not found. Run Final_Matrix_Coefficients_MCK.m first.');
end
load('../outputs/sysid/MasterSystemData.mat'); % M_mat, C_mat, K_mat (6x6 diagonal)

PosTable = readtable('../outputs/clean/Clean_6DOF_Position.csv');
VelTable = readtable('../outputs/clean/Clean_6DOF_Velocity.csv');
AccTable = readtable('../outputs/clean/Clean_6DOF_Acceleration.csv');

t_data = PosTable.Time;  % 5:0.1:55 (501 points)

% Extract 6-DOF arrays [time x DOF]
X_aqwa   = table2array(PosTable(:, 2:7));  % [Surge Sway Heave Roll Pitch Yaw]
dX_aqwa  = table2array(VelTable(:, 2:7));
ddX_aqwa = table2array(AccTable(:, 2:7));

n_steps = length(t_data);
fprintf('Loaded AQWA response data: %d time steps (%.1f s to %.1f s)\n', n_steps, t_data(1), t_data(end));

%% 2. RECONSTRUCT TOTAL EXCITATION FORCE FROM AQWA RESPONSE
% F_total = M * ddX + C * dX + K * X
% This is the combined wave + wind + mooring-nonlinearity force that
% was acting on the platform in the AQWA simulation.
F_total = zeros(n_steps, 6);

for i = 1:n_steps
    F_total(i, :) = (M_mat * ddX_aqwa(i, :)' + ...
                     C_mat * dX_aqwa(i, :)'  + ...
                     K_mat * X_aqwa(i, :)')';
end

fprintf('Reconstructed total AQWA excitation force.\n');

%% 3. ESTIMATE AND SUBTRACT AQWA STEADY WIND DRAG
% The AQWA simulation included a steady 6 m/s wind at 135 deg. This drag
% is embedded in F_total. We estimate and subtract it so that the DMST
% aero forces in Coupled_Simulation.m replace (not add to) this wind load.
%
% Estimated AQWA wind drag: F_drag = 0.5 * rho_air * Cd * A_frontal * V^2
%   rho_air = 1.225 kg/m^3, Cd ≈ 1.0 (bluff body), A_frontal ≈ 1.5 m^2
%   F_drag = 0.5 * 1.225 * 1.0 * 1.5 * 6^2 ≈ 33 N (at 135 deg)
%
% Decompose into surge/sway at 135 deg: cos(135) = -0.707, sin(135) = 0.707
rho_air    = 1.225;
Cd_aqwa    = 1.0;
A_frontal  = 1.5;
V_wind_aqwa = 6;
theta_wind_aqwa = deg2rad(135);
F_drag_total = 0.5 * rho_air * Cd_aqwa * A_frontal * V_wind_aqwa^2;

F_wind_aqwa = zeros(n_steps, 6);
F_wind_aqwa(:, 1) = F_drag_total * cos(theta_wind_aqwa);  % Surge
F_wind_aqwa(:, 2) = F_drag_total * sin(theta_wind_aqwa);  % Sway
% Heave, Roll, Pitch, Yaw wind drag assumed negligible at model scale

%% 4. SAVE EXTRACTED WAVE EXCITATION (with AQWA wind removed)
F_wave = F_total - F_wind_aqwa;
fprintf('Estimated AQWA steady wind drag: %.1f N (total), removed from F_wave.\n', F_drag_total);
fprintf('  Surge contribution: %.1f N, Sway contribution: %.1f N\n', ...
    F_drag_total*cos(theta_wind_aqwa), F_drag_total*sin(theta_wind_aqwa));

%% 4. VISUALIZATION
figure('Name', 'Extracted Wave Excitation Forces', 'Color', 'w', 'Position', [50 50 1200 850]);
force_labels = {'F_{Surge} [N]', 'F_{Sway} [N]', 'F_{Heave} [N]', ...
                'M_{Roll} [Nm]', 'M_{Pitch} [Nm]', 'M_{Yaw} [Nm]'};
colors = {[0 0.45 0.74], [0.85 0.33 0.1], [0.93 0.69 0.13], ...
          [0.49 0.18 0.56], [0.47 0.67 0.19], [0.3 0.75 0.93]};

for i = 1:6
    subplot(3, 2, i);
    plot(t_data, F_wave(:, i), 'LineWidth', 1.5, 'Color', colors{i});
    title(['\bf ' force_labels{i}]);
    xlabel('Time [s]'); ylabel('Force / Moment');
    grid on; xlim([t_data(1) t_data(end)]);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
    th = get(gca, 'Title'); set(th, 'Color', 'k');
end
h_sg = sgtitle('\bf Extracted AQWA Excitation Forces (Wave + Steady Wind)', 'FontSize', 14);
set(h_sg, 'Color', 'k');

%% 5. COMPARE WITH OLD K*eta APPROXIMATION (if available)
if isfile('../outputs/forces/Clean_Wave_Force.csv')
    OldWaveTable = readtable('../outputs/forces/Clean_Wave_Force.csv');
    F_wave_old = table2array(OldWaveTable(:, 2:7));

    figure('Name', 'Wave Force: Extracted vs K*eta Approximation', ...
           'Color', 'w', 'Position', [100 100 1200 850]);

    for i = 1:6
        subplot(3, 2, i);
        plot(t_data, F_wave(:, i), 'b-', 'LineWidth', 1.5); hold on;
        plot(t_data, F_wave_old(:, i), 'r--', 'LineWidth', 1.5);
        title(['\bf ' force_labels{i}]);
        xlabel('Time [s]'); ylabel('Force / Moment');
        legend('Extracted (AQWA)', 'K\times\eta approx', 'Location', 'best');
        grid on; xlim([t_data(1) t_data(end)]);
        set(gca, 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
        th = get(gca, 'Title'); set(th, 'Color', 'k');
    end
    h_sg2 = sgtitle('\bf Comparison: Extracted Wave Excitation vs K\times\eta Approximation', ...
            'FontSize', 14);
    set(h_sg2, 'Color', 'k');

    % Print RMS difference
    fprintf('\n--- Comparison with old K*eta approximation ---\n');
    for i = 1:6
        rms_diff = sqrt(mean((F_wave(:, i) - F_wave_old(:, i)).^2));
        pct_range = 100 * rms_diff / (max(F_wave(:, i)) - min(F_wave(:, i)));
        fprintf('  %s: RMS diff = %.2f (%.1f%% of range)\n', ...
                force_labels{i}, rms_diff, pct_range);
    end
end

%% 6. EXPORT
headers = {'Time', 'F_Surge', 'F_Sway', 'F_Heave', 'F_Roll', 'F_Pitch', 'F_Yaw'};
T_wave = array2table([t_data, F_wave], 'VariableNames', headers);
writetable(T_wave, '../outputs/forces/Clean_Wave_Force_Extracted.csv');
fprintf('\nExported Clean_Wave_Force_Extracted.csv for coupling simulation.\n');
