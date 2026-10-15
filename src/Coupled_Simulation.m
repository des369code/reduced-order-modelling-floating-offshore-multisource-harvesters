% =========================================================================
% PHASE D: MASTER 6-DOF COUPLED SIMULATION (AERO + HYDRO + MOORING)
% =========================================================================
clear; clc; close all;

%% 1. LOAD SYSTEM DNA (M, C, K)
if ~isfile('../outputs/sysid/MasterSystemData.mat')
    error('MasterSystemData.mat not found. Please run your Matrix Assembly script first.');
end
load('../outputs/sysid/MasterSystemData.mat'); % Loads M_mat, C_mat, K_mat

%% 2. LOAD ENVIRONMENTAL FORCES (135-deg COLLINEAR)
fprintf('Loading force datasets...\n');
WaveTable = readtable('../outputs/forces/Clean_Wave_Force_Extracted.csv');
AeroTable = readtable('../outputs/forces/Clean_Aero_Force_Variable.csv');

% Extract time vectors (may have different sampling rates)
t_wave = WaveTable.Time;
t_aero = AeroTable.Time;
F_wave = table2array(WaveTable(:, 2:7));
F_aero = table2array(AeroTable(:, 2:7));

% --- Environmental Direction ---
% Wave direction is LOCKED to the AQWA Design 4 data (135 deg).
% The extracted F_wave was back-calculated from AQWA response at this heading.
% Rotating F_wave would be physically invalid — the hydrodynamic RAO
% (Response Amplitude Operator) is direction-dependent.
% To simulate different wave headings, new AQWA runs would be needed.
theta_wave_deg = 135;  % Wave direction [deg] — baked into Clean_Wave_Force_Extracted.csv

%% 3. INITIALIZE ODE45 SOLVER
% Start from AQWA state at t=5s to avoid artificial startup transient.
% By t=5s, the AQWA simulation has reached quasi-steady oscillatory state.
PosTable = readtable('../outputs/clean/Clean_6DOF_Position.csv');
VelTable = readtable('../outputs/clean/Clean_6DOF_Velocity.csv');
Pos0 = table2array(PosTable(1, 2:7));  % [Surge Sway Heave Roll Pitch Yaw] at t=5s
Vel0 = table2array(VelTable(1, 2:7));
Y0 = [Pos0'; Vel0'];  % 12x1 state vector
tspan = [max(t_wave(1), t_aero(1)), min(t_wave(end), t_aero(end))];

fprintf('Solving 6-DOF Equations of Motion...\n');
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);

% Run Solver
[t_out, Y_out] = ode45(@(t,Y) EoM_Coupled(t, Y, M_mat, C_mat, K_mat, t_wave, F_wave, t_aero, F_aero), tspan, Y0, options);

%% 4. POST-PROCESSING
Pos = Y_out(:, 1:6); % Extract [Surge, Sway, Heave, Roll, Pitch, Yaw]

% Convert Rotational DOFs from Radians to Degrees for easier analysis
Pos(:, 4:6) = rad2deg(Pos(:, 4:6));

%% 5. FINAL VISUALIZATION
figure('Name', '6-DOF Hybrid Platform Response', 'Color', 'w', 'Position', [50 50 1200 850]);
labels = {'Surge [m]', 'Sway [m]', 'Heave [m]', 'Roll [deg]', 'Pitch [deg]', 'Yaw [deg]'};
colors = {[0 0.45 0.74], [0.85 0.33 0.1], [0.93 0.69 0.13], [0.49 0.18 0.56], [0.47 0.67 0.19], [0.3 0.75 0.93]};

for i = 1:6
    subplot(3, 2, i);
    plot(t_out, Pos(:, i), 'LineWidth', 2, 'Color', colors{i});
    title(['\bf ' labels{i}]);
    xlabel('Time [s]'); ylabel('Amplitude');
    grid on;
    xlim(tspan);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
    th = get(gca, 'Title'); set(th, 'Color', 'k');
end

h_sg = sgtitle(sprintf('\\bf Final Coupled Response: AQWA-Extracted Waves (%d^\\circ) + VAWT Wind Load', theta_wave_deg), 'FontSize', 14);
set(h_sg, 'Color', 'k');

%% 6. EXPORT RESULTS FOR REPORT
ResultTable = array2table([t_out, Pos], 'VariableNames', {'Time','Surge','Sway','Heave','Roll','Pitch','Yaw'});
writetable(ResultTable, '../outputs/results/Final_Simulation_Results.csv');
fprintf('Success! Results saved to ../outputs/results/Final_Simulation_Results.csv.\n');

% =========================================================================
% LOCAL FUNCTION: THE DYNAMIC ENGINE
% =========================================================================
function dYdt = EoM_Coupled(t, Y, M, C, K, t_wave, Fw_ref, t_aero, Fa_ref)
    % Extract current state
    X = Y(1:6);
    dX = Y(7:12);

    % Interpolate wave and aero forces at current time (each on its own time grid)
    Fw = interp1(t_wave, Fw_ref, t, 'linear', 'extrap')';
    Fa = interp1(t_aero, Fa_ref, t, 'linear', 'extrap')';
    
    % The Equation: ddX = M^-1 * (F_env - C*dX - K*X)
    % We use matrix left-division (\) for numerical stability
    ddX = M \ (Fw + Fa - C*dX - K*X);
    
    % Return derivatives [Velocity; Acceleration]
    dYdt = [dX; ddX];
end