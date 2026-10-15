% =========================================================================
% PHASE D: 6-DOF COUPLED SIMULATION — EXTREME WIND (20 m/s)
% Usage: Coupled_Simulation_20ms('canonical')
%        Coupled_Simulation_20ms('ballast')
%        Coupled_Simulation_20ms('wide')
%        Coupled_Simulation_20ms('combined')
% =========================================================================
function Coupled_Simulation_20ms(case_name)

%% 1. LOAD SYSTEM DNA (M, C, K)
if ~isfile('../outputs/sysid/MasterSystemData.mat')
    error('MasterSystemData.mat not found. Run your matrix assembly script first.');
end
load('../outputs/sysid/MasterSystemData.mat');

%% 2. LOAD ENVIRONMENTAL FORCES
fprintf('\n=== 20 m/s COUPLED SIMULATION: %s ===\n', upper(case_name));
fprintf('Loading force datasets...\n');
WaveTable = readtable('../outputs/forces/Clean_Wave_Force_Extracted.csv');
AeroTable = readtable('../outputs/forces/Clean_Aero_Force_Variable_20ms.csv');

t_wave = WaveTable.Time;
t_aero = AeroTable.Time;
F_wave = table2array(WaveTable(:, 2:7));
F_aero = table2array(AeroTable(:, 2:7));

%% 3. INITIALIZE ODE45 SOLVER
PosTable = readtable('../outputs/clean/Clean_6DOF_Position.csv');
VelTable = readtable('../outputs/clean/Clean_6DOF_Velocity.csv');
Pos0 = table2array(PosTable(1, 2:7));
Vel0 = table2array(VelTable(1, 2:7));
Y0 = [Pos0'; Vel0'];
tspan = [max(t_wave(1), t_aero(1)), min(t_wave(end), t_aero(end))];

fprintf('Solving 6-DOF Equations of Motion...\n');
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);
[t_out, Y_out] = ode45(@(t,Y) EoM_Coupled(t, Y, M_mat, C_mat, K_mat, t_wave, F_wave, t_aero, F_aero), tspan, Y0, options);

%% 4. POST-PROCESSING
Pos = Y_out(:, 1:6);
Pos(:, 4:6) = rad2deg(Pos(:, 4:6));

%% 5. PRINT STATS
fprintf('\n--- %s RESULTS ---\n', upper(case_name));
labels = {'Surge[m]','Sway[m]','Heave[m]','Roll[deg]','Pitch[deg]','Yaw[deg]'};
for i = 1:6
    fprintf('  %s: mean=%+.3f  min=%+.3f  max=%+.3f  std=%.3f\n', ...
        labels{i}, mean(Pos(:,i)), min(Pos(:,i)), max(Pos(:,i)), std(Pos(:,i)));
end

%% 6. PLOT
figure('Name', sprintf('20 m/s — %s', case_name), 'Color', 'w', 'Position', [50 50 1200 850]);
label_names = {'Surge [m]', 'Sway [m]', 'Heave [m]', 'Roll [deg]', 'Pitch [deg]', 'Yaw [deg]'};
colors = {[0 0.45 0.74], [0.85 0.33 0.1], [0.93 0.69 0.13], [0.49 0.18 0.56], [0.47 0.67 0.19], [0.3 0.75 0.93]};
for i = 1:6
    subplot(3, 2, i);
    plot(t_out, Pos(:, i), 'LineWidth', 2, 'Color', colors{i});
    title(['\bf ' label_names{i}]); xlabel('Time [s]'); ylabel('Amplitude'); grid on; xlim(tspan);
end
sgtitle(sprintf('\\bf 20 m/s Extreme Wind — %s', case_name), 'FontSize', 14);

%% 7. SAVE
outpath = sprintf('../outputs/results/Final_Simulation_Results_20ms_%s.csv', case_name);
ResultTable = array2table([t_out, Pos], 'VariableNames', {'Time','Surge','Sway','Heave','Roll','Pitch','Yaw'});
writetable(ResultTable, outpath);
fprintf('\nSaved to: %s\n', outpath);

end

%% LOCAL: Equations of Motion
function dYdt = EoM_Coupled(t, Y, M, C, K, t_wave, Fw_ref, t_aero, Fa_ref)
    X = Y(1:6); dX = Y(7:12);
    Fw = interp1(t_wave, Fw_ref, t, 'linear', 'extrap')';
    Fa = interp1(t_aero, Fa_ref, t, 'linear', 'extrap')';
    dYdt = [dX; M \ (Fw + Fa - C*dX - K*X)];
end
