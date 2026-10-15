% =========================================================================
% PHASE 2: DATA SYNCHRONIZATION & KINEMATICS (DUPLICATE FIX)
% =========================================================================
clear; clc; close all;

%% 1. DEFINE MASTER TIME VECTOR (THE OVERLAP WINDOW)
dt = 0.1; % Choose a clean time step (e.g., 0.1 seconds)

% Set the bounds to the strictly overlapping timeframe!
t_start = 5;
t_end = 55;
t_master = (t_start : dt : t_end)';

%% 2. LOCATE DATASET FOLDER
% Raw datasets are in ../data/raw/; clean outputs go to ../outputs/clean/
script_dir = fileparts(mfilename('fullpath'));
dataset_dir = fullfile(script_dir, '..', 'data', 'raw');
output_dir  = fullfile(script_dir, '..', 'outputs', 'clean');

%% 3. LOAD AND INTERPOLATE 6DOF MOTION DATA
% Raw data spans 0-60s with uneven time points from WebPlotDigitizer
filenames_6dof = {'Surge_Dataset.csv', 'Sway_Dataset.csv', 'Heave_Dataset.csv', ...
                  'Roll_Dataset.csv', 'Pitch_Dataset.csv', 'Yaw_Dataset.csv'};
dof_labels = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
X_pos = zeros(length(t_master), 6);

for i = 1:6
    rawData = readmatrix(fullfile(dataset_dir, filenames_6dof{i}));
    t_raw = rawData(:,1);
    val_raw = rawData(:,2);

    % Find unique time points to prevent interp1 error from digitizer duplicates
    [t_unique, unique_idx] = unique(t_raw);
    val_unique = val_raw(unique_idx);

    % Interpolate onto the uniform master time grid
    X_pos(:, i) = interp1(t_unique, val_unique, t_master, 'linear');
end

% ANSYS AQWA outputs rotational motions in DEGREES by default.
% Convert Roll (col 4), Pitch (col 5), Yaw (col 6) to RADIANS
% so they are consistent with the ODE's M,C,K which use SI units
% (kg*m^2 for inertia, N*m/rad for stiffness, N*m*s/rad for damping).
X_pos(:, 4:6) = deg2rad(X_pos(:, 4:6));

%% 4. LOAD AND INTERPOLATE MOORING LINE TENSIONS
% Cable data spans approximately 5-55s — the AQWA mooring initialization window
filenames_mooring = {'Cable1_Tension_Dataset.csv', 'Cable2_Tension_Dataset.csv', ...
                     'Cable3_Tension_Dataset.csv', 'Cable4_Tension_Dataset.csv'};
T_mooring = zeros(length(t_master), 4);

for i = 1:4
    rawData = readmatrix(fullfile(dataset_dir, filenames_mooring{i}));
    t_raw = rawData(:,1);
    val_raw = rawData(:,2);

    [t_unique, unique_idx] = unique(t_raw);
    val_unique = val_raw(unique_idx);

    T_mooring(:, i) = interp1(t_unique, val_unique, t_master, 'linear');
end

%% 5. CALCULATE VELOCITIES & ACCELERATIONS (with low-pass filtering)
% Apply 4th-order Butterworth low-pass filter before differentiation
% to suppress digitization noise amplification (noise * 1/dt per derivative).
% Cutoff at 1 Hz — well above platform natural frequencies (0.1-0.15 Hz)
% and wave frequency (~0.1 Hz), but below digitization noise.
fc = 1.0;  % Cutoff frequency [Hz]
Fs = 1/dt;  % Sampling frequency [Hz]
[b, a] = butter(4, fc/(Fs/2), 'low');

X_pos_filt = zeros(size(X_pos));
for i = 1:6
    X_pos_filt(:, i) = filtfilt(b, a, X_pos(:, i));  % Zero-phase filtering
end

% Central difference numerical differentiation via gradient()
X_vel = zeros(size(X_pos_filt));
X_acc = zeros(size(X_pos_filt));

for i = 1:6
    X_vel(:, i) = gradient(X_pos_filt(:, i), dt);
    X_acc(:, i) = gradient(X_vel(:, i), dt);
end

%% 6. EXPORT CLEANED DATA TO CSV
fprintf('Data successfully cleaned, cropped, and synchronized to %d time steps (from %.1fs to %.1fs).\n', length(t_master), t_start, t_end);

% 6-DOF Position (low-pass filtered to suppress WebPlotDigitizer noise)
T_pos = array2table([t_master, X_pos_filt], 'VariableNames', ['Time', dof_labels]);
writetable(T_pos, fullfile(output_dir, 'Clean_6DOF_Position.csv'));
fprintf('  -> Exported Clean_6DOF_Position.csv\n');

% 6-DOF Velocity
T_vel = array2table([t_master, X_vel], 'VariableNames', ['Time', dof_labels]);
writetable(T_vel, fullfile(output_dir, 'Clean_6DOF_Velocity.csv'));
fprintf('  -> Exported Clean_6DOF_Velocity.csv\n');

% 6-DOF Acceleration
T_acc = array2table([t_master, X_acc], 'VariableNames', ['Time', dof_labels]);
writetable(T_acc, fullfile(output_dir, 'Clean_6DOF_Acceleration.csv'));
fprintf('  -> Exported Clean_6DOF_Acceleration.csv\n');

% Mooring Tensions
cable_labels = {'Line1','Line2','Line3','Line4'};
T_moor = array2table([t_master, T_mooring], 'VariableNames', ['Time', cable_labels]);
writetable(T_moor, fullfile(output_dir, 'Clean_Mooring_Tension.csv'));
fprintf('  -> Exported Clean_Mooring_Tension.csv\n');

fprintf('All clean data files written to: %s\n', output_dir);

%% 7. VERIFICATION PLOT

figure('Name', 'Synchronized Surge Data Check', 'Color', 'w');
subplot(3,1,1);
plot(t_master, X_pos(:,1), 'b', 'LineWidth', 1.5); title('Surge Position (x)'); grid on;
xlim([t_start t_end]);
subplot(3,1,2);
plot(t_master, X_vel(:,1), 'r', 'LineWidth', 1.5); title('Surge Velocity (\dot{x})'); grid on;
xlim([t_start t_end]);
subplot(3,1,3);
plot(t_master, X_acc(:,1), 'k', 'LineWidth', 1.5); title('Surge Acceleration (\ddot{x})'); grid on;
xlabel('Time [s]');
xlim([t_start t_end]);