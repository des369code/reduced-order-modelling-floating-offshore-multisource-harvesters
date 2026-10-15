% =========================================================================
% VARIABLE WIND (QUASI-STEADY) DMST SIMULATION
% Features: TSR-optimized Omega, Aero-Map, Kaimal-spectrum turbulent wind.
% =========================================================================
clear; clc; close all;

%% 1. USER INPUTS
% --- Turbine Geometry ---
N_darrieus = 3;
c = 0.302;
R_darrieus = 3.765;
H_darrieus = 3.203;
BladeShape = 'TROPOSKEIN';
R_savonius = 0.3925;
H_sav_total = 2.109;
N_sav_levels = 2;
Sav_Twist = 90;
Z_eq = 1.6015;     % Midpoint of hybrid turbine system (Darrieus equator) [m]
Z_ref = 1.6015;    % Wind shear reference height [m]
alpha_w = 0.11;    % Offshore wind shear exponent
rho = 1.225;       % Air density [kg/m^3]

% --- DMST Simulation Parameters ---
Omega = 5;  % Placeholder — will be overwritten by TSR optimization below
d_theta_deg = 5;
N_streamtubes = 10;
MaxIter = 50;
Tol = 1e-4;

% --- Environmental Direction ---
theta_wind_deg = 135;        % Wind incident direction [deg, 0=North, 90=East]
                              % Default 135 deg matches AQWA Design 4, Case 1.
                              % Change this to test other wind headings.
                              % Wave direction is separately controlled in
                              % Coupled_Simulation.m (locked to AQWA data).
theta_wind = deg2rad(theta_wind_deg);

% Load Airfoil
filename = '../data/airfoil/NACA0009_360_Full_New2.csv';
T = readtable(filename);
data = table2array(T); if isnan(data(1,1)); data(1,:) = []; end
[alpha_vec, idx] = unique(data(:,1));
Cl_vec = data(idx,2);
C_d_vec = max(data(idx,3), 0.005);

% Geometry Prep
d_zeta = 1 / N_streamtubes;  % Normalized streamtube height increment d(z/H) = 1/N
zeta_nodes = linspace(-1 + d_zeta, 1 - d_zeta, N_streamtubes);
slice_height_m = H_darrieus / N_streamtubes; S_darrieus = 0;

for z = zeta_nodes
   if strcmp(BladeShape, 'STRAIGHT'); r=R_darrieus; else; r=R_darrieus*cos(asin(z)); end
   S_darrieus = S_darrieus + (2 * r * slice_height_m);
end
z_fine = linspace(-1, 1, 100); L_sum = 0;
for k = 1:99
   z1 = z_fine(k); z2 = z_fine(k+1);
   if strcmp(BladeShape, 'STRAIGHT'); r1=R_darrieus; r2=R_darrieus; else; r1=R_darrieus*cos(asin(z1)); r2=R_darrieus*cos(asin(z2)); end
   L_sum = L_sum + sqrt(((z2-z1)*(H_darrieus/2))^2 + (r2-r1)^2);
end
AR = L_sum / c;
A_swept = 2 * R_darrieus * H_darrieus;  % Projected swept area [m^2]

%% 2. TSR OPTIMIZATION — FIND OPTIMAL ROTATIONAL SPEED
% For a fixed-speed VAWT, Omega is constant and should be chosen to
% maximize Cp at the rated (mean) wind speed. This is standard design
% practice (Paraschivoiu 2002; Manwell et al., "Wind Energy Explained").
%
% Blade solidity: sigma = N*c/R = 3*0.302/3.765 = 0.241.
% At sigma = 0.24, peak Cp is expected at TSR ~ 2.5–3.5
% (Paraschivoiu 2002, Fig. 6.15). We verify by scanning the full
% DMST model (Darrieus + Savonius) at V = 6 m/s (AQWA Design 4, Case 1).
%
% Omega = 5 rad/s gives TSR = 5*3.765/6 = 3.14 — within the valid
% range of the DMST momentum theory (lambda > 1.5).

V_rated = 6;   % Rated wind speed for design-point optimization [m/s] (matches AQWA Design 4)
scan_TSR = 1.5:0.25:5.0;
scan_Cp  = zeros(size(scan_TSR));
scan_Tmean = zeros(size(scan_TSR));

fprintf('==============================================================\n');
fprintf(' TSR OPTIMIZATION AT V = %.0f m/s\n', V_rated);
fprintf(' Blade solidity sigma = N*c/R = %.3f\n', N_darrieus*c/R_darrieus);
fprintf(' Expected optimal TSR ~ 2.5–3.5 (Paraschivoiu 2002, Fig. 6.15)\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' %6s  %8s  %10s  %10s\n', 'TSR', 'Omega', 'Mean Torq', 'Cp');
fprintf('--------------------------------------------------------------\n');

% Pre-compute azimuth grids (same for all TSR values)
theta_up_opt = (-pi/2 : deg2rad(d_theta_deg) : pi/2)';
theta_dw_opt = (pi/2 : deg2rad(d_theta_deg) : 3*pi/2)';
theta_map_opt = (0:1:359)';

for k = 1:length(scan_TSR)
    tsr = scan_TSR(k);
    Omega_scan = tsr * V_rated / R_darrieus;

    [Torque_Dar_Total, ~, Torque_Sav_Total, ~, ~] = ...
        run_DMST_revolution(V_rated, Omega_scan, theta_up_opt, theta_dw_opt, theta_map_opt);

    Torque_Hybrid = Torque_Dar_Total + Torque_Sav_Total;
    T_mean = mean(Torque_Hybrid);
    P_aero = T_mean * Omega_scan;
    P_wind = 0.5 * rho * A_swept * V_rated^3;
    Cp = P_aero / P_wind;

    scan_Cp(k) = Cp;
    scan_Tmean(k) = T_mean;

    fprintf('  %5.2f   %7.2f   %9.1f   %9.4f\n', tsr, Omega_scan, T_mean, Cp);
end

% Find optimal TSR
[max_Cp, idx_opt] = max(scan_Cp);
Best_TSR = scan_TSR(idx_opt);
Omega = Best_TSR * V_rated / R_darrieus;

fprintf('--------------------------------------------------------------\n');
fprintf(' OPTIMAL: TSR = %.2f, Omega = %.2f rad/s (%.1f RPM), Cp = %.4f\n', ...
    Best_TSR, Omega, Omega*60/(2*pi), max_Cp);
fprintf(' At Omega=5 rad/s: TSR=%.2f (momentum theory invalid for TSR < 1.5)\n', 5*R_darrieus/V_rated);
fprintf('==============================================================\n\n');

% Plot Cp-TSR curve
figure('Name', 'Cp-TSR Optimization', 'Color', 'w', 'Position', [100 100 700 500]);
yyaxis left;
plot(scan_TSR, scan_Cp, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Tip-Speed Ratio \lambda'); ylabel('Power Coefficient C_p');
grid on; hold on;
plot(Best_TSR, max_Cp, 'r*', 'MarkerSize', 15);
set(gca, 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);

% Place annotation intelligently — left of peak if near right edge, else right
if Best_TSR > mean(scan_TSR([1 end]))
    x_off = -0.7; ha = 'right';
else
    x_off =  0.3; ha = 'left';
end
text(Best_TSR + x_off, max_Cp * 0.92, ...
    sprintf('Optimal: \\lambda=%.2f, C_p=%.4f\n\\Omega=%.1f rad/s', ...
    Best_TSR, max_Cp, Omega), ...
    'FontSize', 11, 'Color', 'k', 'BackgroundColor', 'w', ...
    'EdgeColor', 'k', 'HorizontalAlignment', ha);

t1 = get(gca, 'Title'); set(t1, 'Color', 'k');
title(sprintf('\\bf Cp-\\lambda Curve (V=%.0f m/s, \\sigma=%.3f)', V_rated, N_darrieus*c/R_darrieus));

yyaxis right;
plot(scan_TSR, scan_Tmean, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6);
ylabel('Mean Torque [N\cdotm]');
set(gca, 'YColor', 'k');

%% 3. BUILD THE AERO-MAP (Run DMST at 5 wind speeds using optimized Omega)
V_test_array = [2, 4, 6, 8, 10];  % covers 3-sigma turbulent range [2.3, 9.7] with headroom

% Pre-compute azimuth grids for aero-map (same as TSR optimization)
theta_up_map = (-pi/2 : deg2rad(d_theta_deg) : pi/2)';
theta_dw_map = (pi/2 : deg2rad(d_theta_deg) : 3*pi/2)';
theta_map = (0:1:359)';

Map_Torque = zeros(length(theta_map), length(V_test_array));
Map_Thrust = zeros(length(theta_map), length(V_test_array));

fprintf('Building Aerodynamic Map at Omega = %.2f rad/s...\n', Omega);
for v_idx = 1:length(V_test_array)
    V_ref = V_test_array(v_idx);
    fprintf('  -> Simulating V = %d m/s (TSR = %.2f)...\n', V_ref, Omega*R_darrieus/V_ref);

    [Torque_Dar_Total, Thrust_Dar_Total, Torque_Sav_Total, Thrust_Sav_Total, ~] = ...
        run_DMST_revolution(V_ref, Omega, theta_up_map, theta_dw_map, theta_map);

    Map_Torque(:, v_idx) = Torque_Dar_Total + Torque_Sav_Total;
    Map_Thrust(:, v_idx) = Thrust_Dar_Total + Thrust_Sav_Total;
end

% Print aero-map summary at V_rated (design point)
fprintf('\nAero-Map Summary at V=%.0f m/s (design point, TSR=%.2f):\n', V_rated, Omega*R_darrieus/V_rated);
fprintf('  Mean Darrieus Torque: %.1f N.m\n', mean(Torque_Dar_Total));
fprintf('  Mean Savonius Torque: %.1f N.m\n', mean(Torque_Sav_Total));
fprintf('  Mean Hybrid  Torque:  %.1f N.m\n', mean(Map_Torque(:,3)));
fprintf('  Mean Hybrid  Thrust:  %.1f N\n', mean(Map_Thrust(:,3)));

% Resolve design-point 6-DOF forces (matches thesis Table 4.3.1)
T_mean_dp = mean(Map_Torque(:,3));
Th_mean_dp = mean(Map_Thrust(:,3));
F_dp = zeros(1,6);
F_dp(1) = Th_mean_dp * cos(theta_wind);      % Surge
F_dp(2) = Th_mean_dp * sin(theta_wind);      % Sway
F_dp(3) = 0;                                   % Heave
F_dp(4) = -(Th_mean_dp * sin(theta_wind)) * Z_eq;  % Roll
F_dp(5) = (Th_mean_dp * cos(theta_wind)) * Z_eq;   % Pitch
F_dp(6) = -T_mean_dp;                            % Yaw

fprintf('\n--- Design-Point 6-DOF Aero Forces (V=%.0f m/s, \\theta_{wind}=%.0f°) ---\n', V_rated, theta_wind_deg);
fprintf('%-20s %10s %10s\n', 'Component', 'Mean', 'Units');
fprintf('%-20s %10.1f %10s\n', 'F_aero,surge', F_dp(1), 'N');
fprintf('%-20s %10.1f %10s\n', 'F_aero,sway',  F_dp(2), 'N');
fprintf('%-20s %10.1f %10s\n', 'F_aero,heave', F_dp(3), 'N');
fprintf('%-20s %10.1f %10s\n', 'M_aero,roll',  F_dp(4), 'N.m');
fprintf('%-20s %10.1f %10s\n', 'M_aero,pitch', F_dp(5), 'N.m');
fprintf('%-20s %10.1f %10s\n', 'M_aero,yaw',   F_dp(6), 'N.m');
fprintf('----------------------------------------------------------------------\n');

% Plot aero-map: torque and thrust vs azimuth at all wind speeds
colors = lines(length(V_test_array));
figure('Name', 'Aero-Map: Torque and Thrust vs Azimuth', 'Color', 'w', ...
    'Position', [100 100 1000 700]);

subplot(2,1,1); hold on;
for k = 1:length(V_test_array)
    plot(theta_map, Map_Torque(:,k), 'Color', colors(k,:), 'LineWidth', 1.5);
end
ylabel('Torque [N\cdotm]');
title(sprintf('\\bf Hybrid Torque vs Azimuth (\\Omega = %.1f rad/s)', Omega));
legend(cellstr(num2str(V_test_array', 'U = %d m/s')), 'Location', 'best');
grid on; xlim([0 360]);
set(gca, 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
t1 = get(gca, 'Title'); set(t1, 'Color', 'k');

subplot(2,1,2); hold on;
for k = 1:length(V_test_array)
    plot(theta_map, Map_Thrust(:,k), 'Color', colors(k,:), 'LineWidth', 1.5);
end
xlabel('Azimuth Angle [deg]'); ylabel('Thrust [N]');
title(sprintf('\\bf Hybrid Thrust vs Azimuth (\\Omega = %.1f rad/s)', Omega));
legend(cellstr(num2str(V_test_array', 'U = %d m/s')), 'Location', 'best');
grid on; xlim([0 360]);
set(gca, 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
t2 = get(gca, 'Title'); set(t2, 'Color', 'k');

%% 4. THE TIME LOOP (Variable Wind Simulation)
fprintf('\nAero-Map complete! Running Time Simulation...\n');

dt_ode = 0.001;
t_ode = (5 : dt_ode : 55)';

% Generate turbulent wind signal using IEC 61400-1 Kaimal spectrum
U_mean  = 6;     % Mean wind speed [m/s] (matches AQWA Design 4, Case 1)
I_ref   = 0.14;  % Turbulence intensity (offshore Class B)
sigma_u = I_ref * U_mean;
L1      = 0.7 * Z_eq;  % Integral length scale ~ hub height (surface layer scaling) [m]

fs   = 1 / dt_ode;
N    = length(t_ode);
f    = (0 : N-1)' * fs / N;
df   = fs / N;

S = zeros(N, 1);
nonzero = f > 0;
S(nonzero) = sigma_u^2 * (4 * L1 / U_mean) ./ ...
             (1 + 6 * f(nonzero) * L1 / U_mean).^(5/3);
S(1) = 0;

S_two_sided = S / 2;
amp = sqrt(S_two_sided * df);

phi = 2 * pi * rand(N, 1);
phi(1) = 0;

H = amp .* exp(1i * phi);
H(1) = 0;
if mod(N, 2) == 0
    H(N/2 + 1) = real(H(N/2 + 1));
    H(N/2+2 : end) = conj(flip(H(2 : N/2)));
else
    H((N+1)/2+1 : end) = conj(flip(H(2 : (N+1)/2)));
end

u_turb = ifft(H, 'symmetric') * N;
V_turbulent = U_mean + 1.2*sin(0.5 * t_ode) + real(u_turb);
V_turbulent = max(min(V_turbulent, V_test_array(end)), V_test_array(1));  % clamp to aero-map bounds

Thrust_ode = zeros(length(t_ode), 1);
Torque_ode = zeros(length(t_ode), 1);

[V_GRID, THETA_GRID] = meshgrid(V_test_array, theta_map);

for i = 1:length(t_ode)
    V_current = V_turbulent(i);
    theta_current = mod(rad2deg(Omega * t_ode(i)), 359);
    Torque_ode(i) = interp2(V_GRID, THETA_GRID, Map_Torque, V_current, theta_current);
    Thrust_ode(i) = interp2(V_GRID, THETA_GRID, Map_Thrust, V_current, theta_current);
end

% Figure 4.5: Turbulent wind field — time series + Kaimal PSD overlay
[Pxx_act, f_act] = pwelch(V_turbulent - mean(V_turbulent), ...
    floor(N/5), floor(N/10), N, fs);
S_theo = sigma_u^2 * (4 * L1 / U_mean) ./ ...
    (1 + 6 * f_act(2:end) * L1 / U_mean).^(5/3);

figure('Name', 'Turbulent Wind Field (Figure 4.5)', 'Color', 'w', ...
    'Position', [100 100 900 650]);

subplot(2,1,1);
plot(t_ode, V_turbulent, 'b', 'LineWidth', 1);
ylabel('Wind Speed [m/s]');
title(sprintf('\\bf Turbulent Wind Speed Time Series (U_{mean}=%.0f m/s, I_{ref}=%.2f)', U_mean, I_ref));
grid on; xlim([5 55]);
set(gca, 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
t1 = get(gca, 'Title'); set(t1, 'Color', 'k');

subplot(2,1,2);
loglog(f_act, Pxx_act, 'b', 'LineWidth', 1.5); hold on;
loglog(f_act(2:end), S_theo, 'r--', 'LineWidth', 2);
xlabel('Frequency [Hz]'); ylabel('PSD [m^2/s]');
title('\bf Kaimal Spectrum: Generated vs Theoretical');
legend('Generated (Welch PSD)', 'Theoretical (Kaimal)', 'Location', 'southwest');
grid on; xlim([f_act(2) f_act(end)]);
set(gca, 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
t2 = get(gca, 'Title'); set(t2, 'Color', 'k');

%% 5. EXPORT TO CSV AND PLOT
% theta_wind is defined in Section 1 — change theta_wind_deg there.
F_aero = zeros(length(t_ode), 6);
F_aero(:,1) = Thrust_ode * cos(theta_wind);
F_aero(:,2) = Thrust_ode * sin(theta_wind);
F_aero(:,3) = 0;
F_aero(:,4) = -(Thrust_ode * sin(theta_wind)) * Z_eq;
F_aero(:,5) = (Thrust_ode * cos(theta_wind)) * Z_eq;
F_aero(:,6) = -Torque_ode;

headers_Aero = {'Time', 'F_Surge', 'F_Sway', 'F_Heave', 'F_Roll', 'F_Pitch', 'F_Yaw'};
writetable(array2table([t_ode, F_aero], 'VariableNames', headers_Aero), '../outputs/forces/Clean_Aero_Force_Variable.csv');
fprintf('Successfully generated Clean_Aero_Force_Variable.csv!\n');

% Print mean 6-DOF aerodynamic forces (matches thesis Table 4.3.1)
fprintf('\n--- Time-Domain Mean 6-DOF Aero Forces (U=%.0f m/s turbulent, \\theta_{wind}=%.0f°) [VARIES PER RUN] ---\n', U_mean, theta_wind_deg);
fprintf('%-20s %10s %10s\n', 'Component', 'Mean', 'Units');
fprintf('%-20s %10.1f %10s\n', 'F_aero,surge', mean(F_aero(:,1)), 'N');
fprintf('%-20s %10.1f %10s\n', 'F_aero,sway',  mean(F_aero(:,2)), 'N');
fprintf('%-20s %10.1f %10s\n', 'F_aero,heave', mean(F_aero(:,3)), 'N');
fprintf('%-20s %10.1f %10s\n', 'M_aero,roll',  mean(F_aero(:,4)), 'N.m');
fprintf('%-20s %10.1f %10s\n', 'M_aero,pitch', mean(F_aero(:,5)), 'N.m');
fprintf('%-20s %10.1f %10s\n', 'M_aero,yaw',   mean(F_aero(:,6)), 'N.m');
fprintf('----------------------------------------------------------------\n');

% Plot
figure('Name', 'Variable Wind Response', 'Color', 'w', 'Position', [100 100 1000 750]);
subplot(3,1,1);
plot(t_ode, V_turbulent, 'b', 'LineWidth', 1);
ylabel('Wind Speed [m/s]'); title(sprintf('\\bf Turbulent Wind Input (Kaimal, U_{mean}=%.0f m/s, I_{ref}=%.2f, \\Omega=%.1f rad/s)', U_mean, I_ref, Omega));
grid on; xlim([5 55]);
set(gca, 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
t1 = get(gca, 'Title'); set(t1, 'Color', 'k');

subplot(3,1,2);
plot(t_ode, Torque_ode, 'r', 'LineWidth', 1);
ylabel('Torque [N\cdotm]'); title('\bf Aerodynamic Torque Response');
grid on; xlim([5 55]);
set(gca, 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
t2 = get(gca, 'Title'); set(t2, 'Color', 'k');

subplot(3,1,3);
plot(t_ode, Thrust_ode, 'Color', [0 0.45 0.74], 'LineWidth', 1);
ylabel('Thrust [N]'); xlabel('Time [s]'); title('\bf Aerodynamic Thrust Response');
grid on; xlim([5 55]);
set(gca, 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
t3 = get(gca, 'Title'); set(t3, 'Color', 'k');

%% HELPER: Airfoil Lookup
function [Cl, Cd] = get_aero(a, vec, Cl_v, Cd_v)
   a = mod(a + 180, 360) - 180;
   Cl = interp1(vec, Cl_v, a, 'linear', 'extrap');
   Cd = interp1(vec, Cd_v, a, 'linear', 'extrap');
   Cl = max(min(Cl, 1.4), -1.4);
end

%% HELPER: DMST One-Revolution Solver (shared by TSR optimization and aero-map)
function [Torque_Dar_Total, Thrust_Dar_Total, Torque_Sav_Total, Thrust_Sav_Total, V_sav_input] = ...
    run_DMST_revolution(V_wind, Omega_rot, theta_up, theta_dw, theta_map)
    % Runs one full DMST revolution (Darrieus upwind + downwind + Savonius)
    % at a single (V_wind, Omega_rot) operating point.
    % Called by both TSR optimization (Section 2) and aero-map build (Section 3).

    % ---- Shared parameters (from caller workspace via persistent-like access) ----
    % These are defined in the main script's Section 1 and accessed here.
    % MATLAB nested-function rules require these to be passed or in scope.
    % We access them as variables from the enclosing script workspace.

    % Assign from caller's workspace via evalin (safe — called only from this script)
    Ns  = evalin('caller', 'N_streamtubes');
    dth = evalin('caller', 'd_theta_deg');
    zn  = evalin('caller', 'zeta_nodes');
    dz  = evalin('caller', 'd_zeta');
    BS  = evalin('caller', 'BladeShape');
    Rd  = evalin('caller', 'R_darrieus');
    Hd  = evalin('caller', 'H_darrieus');
    ch  = evalin('caller', 'c');
    Nd  = evalin('caller', 'N_darrieus');
    ARg = evalin('caller', 'AR');
    Rs  = evalin('caller', 'R_savonius');
    Hs  = evalin('caller', 'H_sav_total');
    Nsl = evalin('caller', 'N_sav_levels');
    Stw = evalin('caller', 'Sav_Twist');
    Ze  = evalin('caller', 'Z_eq');
    Zr  = evalin('caller', 'Z_ref');
    aw  = evalin('caller', 'alpha_w');
    rh  = evalin('caller', 'rho');
    MI  = evalin('caller', 'MaxIter');
    Tl  = evalin('caller', 'Tol');
    av  = evalin('caller', 'alpha_vec');
    Cv  = evalin('caller', 'Cl_vec');
    Dv  = evalin('caller', 'C_d_vec');

    TS_Torque_Up = zeros(size(theta_up));
    TS_Thrust_Up = zeros(size(theta_up));
    TS_Torque_Dw = zeros(size(theta_dw));
    TS_Thrust_Dw = zeros(size(theta_dw));
    Sum_Ve = 0;

    for i = 1:Ns
        Zeta = zn(i);
        if strcmp(BS, 'STRAIGHT')
            r = Rd; eta = 1; delta = 0;
        else
            r = Rd * cos(asin(Zeta)); eta = r/Rd; delta = asin(Zeta);
        end
        Z_i = Ze + Zeta * (Hd/2);
        V_inf = (Z_i <= 0.1) * 0.001 + (Z_i > 0.1) * V_wind * (Z_i/Zr)^aw;

        % Upwind induction
        u = 1.0;
        for iter = 1:MI
            sum_f = 0;
            for j = 1:length(theta_up)
                th = theta_up(j);
                if abs(cos(th)) < 1e-5; continue; end
                V = max(u * V_inf, 0.1);
                X = (Omega_rot * r) / V;
                W = V * sqrt((X - sin(th))^2 + (cos(th)*cos(delta))^2);
                al = rad2deg(asin((cos(th)*cos(delta)) / max(W/V, 0.1)));
                [Cl, Cd_2d] = get_aero(al, av, Cv, Dv);
                Cd_i = Cd_2d + (Cl^2) / (pi * ARg * 0.9);
                Cn = Cl*cosd(al) + Cd_i*sind(al);
                Ct = Cl*sind(al) - Cd_i*cosd(al);
                N_fact = cos(th)/abs(cos(th));
                T_fact = sin(th) / (abs(cos(th))*cos(delta) + 1e-6);
                f_local = (Cn*N_fact + Ct*T_fact) * (W/max(V,0.1))^2;
                sum_f = sum_f + f_local;
            end
            f_val = (Nd*ch)/(8*pi*Rd) * sum_f * deg2rad(dth);
            u_new = (pi*eta) / (f_val + pi*eta);
            if iter > 1 && abs(u_new - u_old) < Tl; break; end
            u_old = u_new;
            u = 0.6*u + 0.4*max(u_new, 0.2);
        end

        for j = 1:length(theta_up)
            th = theta_up(j);
            if abs(cos(th)) < 1e-5; continue; end
            V = max(u * V_inf, 0.1);
            X = (Omega_rot * r) / V;
            W = V * sqrt((X - sin(th))^2 + (cos(th)*cos(delta))^2);
            al = rad2deg(asin((cos(th)*cos(delta)) / max(W/V, 0.1)));
            [Cl, Cd_2d] = get_aero(al, av, Cv, Dv);
            Cd_i = Cd_2d + (Cl^2) / (pi * ARg * 0.9);
            Ct = Cl*sind(al) - Cd_i*cosd(al);
            Cn = Cl*cosd(al) + Cd_i*sind(al);
            TS_Torque_Up(j) = TS_Torque_Up(j) + Ct*(W/V_wind)^2*(eta/max(cos(delta),1e-6))*dz;
            TS_Thrust_Up(j) = TS_Thrust_Up(j) + (Cn*cos(th) + Ct*sin(th))*(W/V_wind)^2*dz;
        end

        % Downwind
        V_e = (2*u - 1) * V_inf;
        Sum_Ve = Sum_Ve + V_e;
        u_p = u;
        for iter = 1:MI
            sum_f = 0;
            for j = 1:length(theta_dw)
                th = theta_dw(j);
                if abs(cos(th)) < 1e-5; continue; end
                V = max(u_p * V_e, 0.1);
                X = (Omega_rot * r) / V;
                W = V * sqrt((X - sin(th))^2 + (cos(th)*cos(delta))^2);
                al = rad2deg(asin((cos(th)*cos(delta)) / max(W/V, 0.1)));
                [Cl, Cd_2d] = get_aero(al, av, Cv, Dv);
                Cd_i = Cd_2d + (Cl^2) / (pi * ARg * 0.9);
                Cn = Cl*cosd(al) + Cd_i*sind(al);
                Ct = Cl*sind(al) - Cd_i*cosd(al);
                N_fact = cos(th)/abs(cos(th));
                T_fact = sin(th) / (abs(cos(th))*cos(delta) + 1e-6);
                f_local = (Cn*N_fact + Ct*T_fact) * (W/max(V,0.1))^2;
                sum_f = sum_f + f_local;
            end
            f_val = (Nd*ch)/(8*pi*Rd) * sum_f * deg2rad(dth);
            u_new = (pi*eta) / (f_val + pi*eta);
            if iter > 1 && abs(u_new - u_p_old) < Tl; break; end
            u_p_old = u_new;
            u_p = 0.6*u_p + 0.4*max(u_new, 0.2);
        end

        for j = 1:length(theta_dw)
            th = theta_dw(j);
            if abs(cos(th)) < 1e-5; continue; end
            V = max(u_p * V_e, 0.1);
            X = (Omega_rot * r) / V;
            W = V * sqrt((X - sin(th))^2 + (cos(th)*cos(delta))^2);
            al = rad2deg(asin((cos(th)*cos(delta)) / max(W/V, 0.1)));
            [Cl, Cd_2d] = get_aero(al, av, Cv, Dv);
            Cd_i = Cd_2d + (Cl^2) / (pi * ARg * 0.9);
            Ct = Cl*sind(al) - Cd_i*cosd(al);
            Cn = Cl*cosd(al) + Cd_i*sind(al);
            TS_Torque_Dw(j) = TS_Torque_Dw(j) + Ct*(W/V_wind)^2*(eta/max(cos(delta),1e-6))*dz;
            TS_Thrust_Dw(j) = TS_Thrust_Dw(j) + (Cn*cos(th) + Ct*sin(th))*(W/V_wind)^2*dz;
        end
    end

    % Merge upwind + downwind, interpolate to 1-degree grid
    theta_combined = mod(rad2deg([theta_up; theta_dw]), 360);
    [theta_sorted, unique_idx] = unique(theta_combined);
    temp_Torque = [TS_Torque_Up; TS_Torque_Dw];
    temp_Thrust = [TS_Thrust_Up; TS_Thrust_Dw];

    Torque_Dar_1B = interp1([theta_sorted; theta_sorted(1)+360], ...
        [temp_Torque(unique_idx); temp_Torque(unique_idx(1))], theta_map, 'linear');
    Thrust_Dar_1B = interp1([theta_sorted; theta_sorted(1)+360], ...
        [temp_Thrust(unique_idx); temp_Thrust(unique_idx(1))], theta_map, 'linear');

    factor_Q = 0.5*rh*V_wind^2 * (2*Rd*Hd) * Rd * ((ch*Hd)/(2*Rd*Hd));  % = 0.5*rho*V^2*c*H*R
    factor_T = 0.5*rh*V_wind^2 * (2*Rd*Hd) * ((ch*Hd)/(2*Rd*Hd));       % = 0.5*rho*V^2*c*H

    % 3-blade superposition
    Torque_Dar_Total = zeros(size(theta_map));
    Thrust_Dar_Total = zeros(size(theta_map));
    for b = 1:Nd
        shift = round((b-1) * (360/Nd));
        Torque_Dar_Total = Torque_Dar_Total + circshift(Torque_Dar_1B * factor_Q, shift);
        Thrust_Dar_Total = Thrust_Dar_Total + circshift(Thrust_Dar_1B * factor_T, shift);
    end

    % Savonius (Chauvin Pressure-Drag Model)
    V_sav_input = Sum_Ve / Ns;
    Torque_Sav_Total = zeros(size(theta_map));
    Thrust_Sav_Total = zeros(size(theta_map));
    for m = 1:length(theta_map)
        for lvl = 1:Nsl
            alpha_sav = deg2rad(theta_map(m) + (lvl-1)*Stw);
            Q_M = rh * Rs^2 * (Hs/Nsl) * 2.0 * ...
                (V_sav_input^2 + 2*Omega_rot^2*Rs^2 ...
                - (pi/2)*V_sav_input*Omega_rot*Rs*cos(alpha_sav) ...
                - 2*V_sav_input*Omega_rot*Rs*sin(alpha_sav));
            Q_D = 2 * rh * Rs^4 * (Hs/Nsl) * Omega_rot^2 * 2.0;
            Torque_Sav_Total(m) = Torque_Sav_Total(m) + (Q_M - Q_D);
            Thrust_Sav_Total(m) = Thrust_Sav_Total(m) + ((Q_M - Q_D) / Rs);
        end
    end
end
