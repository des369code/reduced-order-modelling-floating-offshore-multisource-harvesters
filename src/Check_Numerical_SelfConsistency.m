% =========================================================================
% CHECK NUMERICAL SELF-CONSISTENCY: ODE WITH F_aero=0 VS AQWA
% =========================================================================
% CAVEAT: This is a NUMERICAL self-consistency check, NOT a physical
% validation. F_wave is computed FROM AQWA data using M, C, K. Feeding
% the same F_wave, M, C, K back through the ODE MUST reproduce the AQWA
% response — this is a mathematical identity for ANY choice of M, C, K
% (not just physically correct ones). A PASS here proves:
%   (a) ODE solver accuracy is adequate
%   (b) No coding errors in the equation of motion
% It does NOT prove M, C, K represent the actual platform physics.
% Any errors in M, C, K are absorbed into F_wave and become invisible.
%
% For physical validation, independently compare:
%   - K_moor against quasi-static pull-test or known geometry
%   - F_wave spectrum against theoretical JONSWAP wave forces
%   - Damping ratios against free-decay analysis
% =========================================================================
clear; clc; close all;

%% 1. LOAD SYSTEM MATRICES AND FORCING
if ~isfile('../outputs/sysid/MasterSystemData.mat')
    error('MasterSystemData.mat not found. Run Final_Matrix_Coefficients_MCK.m first.');
end
load('../outputs/sysid/MasterSystemData.mat'); % M_mat, C_mat, K_mat

if ~isfile('../outputs/forces/Clean_Wave_Force_Extracted.csv')
    error('Clean_Wave_Force_Extracted.csv not found. Run Extract_Wave_Excitation.m first.');
end
WaveTable = readtable('../outputs/forces/Clean_Wave_Force_Extracted.csv');
t_ref = WaveTable.Time;
F_wave_ref = table2array(WaveTable(:, 2:7));

%% 2. LOAD AQWA REFERENCE DATA (GROUND TRUTH)
PosAQWA = readtable('../outputs/clean/Clean_6DOF_Position.csv');
t_aqwa = PosAQWA.Time;
X_aqwa = table2array(PosAQWA(:, 2:7));

%% 3. RUN ODE45 WITH WAVE-ONLY FORCING (F_aero = 0)
% Start from AQWA initial conditions at t=5s (not zero) to avoid
% an artificial startup transient that would inflate error metrics.
% The AQWA simulation had already reached quasi-steady state by t=5s.
VelAQWA = readtable('../outputs/clean/Clean_6DOF_Velocity.csv');
dX_aqwa_ic = table2array(VelAQWA(1, 2:7));  % Velocities at first timestep
X_aqwa_ic  = X_aqwa(1, :);                   % Positions at first timestep

Y0 = [X_aqwa_ic'; dX_aqwa_ic'];  % 12x1 state vector
tspan = [t_ref(1), t_ref(end)];

fprintf('Solving ODE with wave-only excitation (F_aero = 0)...\n');
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);

[t_out, Y_out] = ode45(@(t, Y) EoM_WaveOnly(t, Y, M_mat, C_mat, K_mat, ...
                        t_ref, F_wave_ref), tspan, Y0, options);

% Extract positions
X_sim = Y_out(:, 1:6);  % [Surge Sway Heave Roll Pitch Yaw]

% Interpolate simulation output onto AQWA time grid for comparison
X_sim_interp = interp1(t_out, X_sim, t_aqwa, 'linear');

fprintf('ODE solved: %d output timesteps.\n', length(t_out));

%% 4. COMPUTE ERROR METRICS FOR EACH DOF
dof_names = {'Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'};
units     = {'m', 'm', 'm', 'rad', 'rad', 'rad'};

fprintf('\n==============================================================\n');
fprintf(' SYSTEM IDENTIFICATION VALIDATION RESULTS\n');
fprintf('==============================================================\n');
fprintf('%-8s %12s %12s %12s %8s\n', 'DOF', 'RMS Error', 'Max Abs Err', 'Signal Range', 'RMS%%');
fprintf('--------------------------------------------------------------\n');

rms_error  = zeros(6, 1);
max_error  = zeros(6, 1);
sig_range  = zeros(6, 1);
rms_pct    = zeros(6, 1);

for i = 1:6
    err = X_sim_interp(:, i) - X_aqwa(:, i);
    rms_error(i) = sqrt(mean(err.^2));
    max_error(i) = max(abs(err));
    sig_range(i) = max(X_aqwa(:, i)) - min(X_aqwa(:, i));
    if sig_range(i) > 1e-6
        rms_pct(i) = 100 * rms_error(i) / sig_range(i);
    else
        rms_pct(i) = 0;
    end

    fprintf('%-8s %12.6f %12.6f %12.6f %7.2f%%', ...
            dof_names{i}, rms_error(i), max_error(i), sig_range(i), rms_pct(i));

    if rms_pct(i) < 5
        fprintf('  PASS\n');
    elseif rms_pct(i) < 10
        fprintf('  WARN\n');
    else
        fprintf('  FAIL\n');
    end
end
fprintf('--------------------------------------------------------------\n');
fprintf('PASS: RMS < 5%% of signal range   WARN: 5-10%%   FAIL: > 10%%\n');
fprintf('==============================================================\n');

%% 5. PLOT OVERLAID COMPARISONS
figure('Name', 'System ID Validation: AQWA vs Reconstructed', ...
       'Color', 'w', 'Position', [50 50 1400 900]);

for i = 1:6
    subplot(3, 2, i);
    plot(t_aqwa, X_aqwa(:, i), 'b-', 'LineWidth', 1.5); hold on;
    plot(t_aqwa, X_sim_interp(:, i), 'r--', 'LineWidth', 1.5);
    title(sprintf('\\bf %s (RMS = %.2f%%)', dof_names{i}, rms_pct(i)));
    xlabel('Time [s]');
    ylabel(sprintf('%s [%s]', dof_names{i}, units{i}));
    legend('AQWA (source)', 'Reconstructed (ODE)', 'Location', 'best');
    grid on;
    xlim([t_aqwa(1) t_aqwa(end)]);
end
sgtitle('\bf System Identification Validation: AQWA Source vs ODE Reconstruction (F_{aero}=0)', ...
        'FontSize', 14);

%% 6. PLOT ERROR TIME HISTORIES
figure('Name', 'Reconstruction Error Over Time', ...
       'Color', 'w', 'Position', [100 100 1400 900]);

for i = 1:6
    subplot(3, 2, i);
    err = X_sim_interp(:, i) - X_aqwa(:, i);
    plot(t_aqwa, err, 'k-', 'LineWidth', 1);
    title(sprintf('\\bf %s Error', dof_names{i}));
    xlabel('Time [s]');
    ylabel(sprintf('Error [%s]', units{i}));
    grid on;
    xlim([t_aqwa(1) t_aqwa(end)]);
    % Add RMS reference line
    hold on;
    yline(rms_error(i), 'r--', 'LineWidth', 1);
    yline(-rms_error(i), 'r--', 'LineWidth', 1);
end
sgtitle('\bf Reconstruction Error: AQWA - ODE (F_{aero}=0)', 'FontSize', 14);

% =========================================================================
% LOCAL FUNCTION: ODE RIGHT-HAND SIDE (WAVE ONLY, NO AERO)
% =========================================================================
function dYdt = EoM_WaveOnly(t, Y, M, C, K, t_ref, Fw_ref)
    X  = Y(1:6);
    dX = Y(7:12);

    % Interpolate extracted wave excitation at current time
    Fw = interp1(t_ref, Fw_ref, t, 'linear', 'extrap')';

    % M * ddX + C * dX + K * X = F_wave
    ddX = M \ (Fw - C * dX - K * X);

    dYdt = [dX; ddX];
end
