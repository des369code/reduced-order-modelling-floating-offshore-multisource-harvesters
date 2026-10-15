% =========================================================================
% MOORING STIFFNESS EXTRACTION: FORCE-DISPLACEMENT REGRESSION
% =========================================================================
clear; clc; close all;

%% 1. LOAD DATASETS
fprintf('Loading data for stiffness analysis...\n');
Pos_Table = readtable('../outputs/clean/Clean_6DOF_Position.csv');
Moor_Table = readtable('../outputs/clean/Clean_Mooring_Tension.csv');

% Extract Time and Motion (LHS of regression)
surge = Pos_Table.Surge; 
sway  = Pos_Table.Sway;
yaw   = Pos_Table.Yaw; 

% Extract Tensions (RHS of regression)
T1 = Moor_Table.Line1;
T2 = Moor_Table.Line2;
T3 = Moor_Table.Line3;
T4 = Moor_Table.Line4;

%% 2. CALCULATE NET RESTORING FORCES
% Mooring geometry (Design 4, 4-line spread mooring):
%   Line 1: NE corner (fairlead at +x, +y)
%   Line 2: SE corner (fairlead at +x, -y)
%   Line 3: NW corner (fairlead at -x, +y)
%   Line 4: SW corner (fairlead at -x, -y)
%
% Each line tension T projects to horizontal plane with fairlead angle
% theta_v (from horizontal), then to platform axes with plan-view angle
% theta_h. For a square spread at 45 deg: theta_h = 45 deg.
%
% FAIRLEAD VERTICAL ANGLE (theta_v):
%   Assumed 45 deg from horizontal (typical for catenary mooring).
%   cos(theta_v) = cos(45) = 0.707 projects tension to horizontal plane.
%
% FAIRLEAD PLAN-VIEW ANGLE (theta_h):
%   Lines radiate at 45 deg from platform axes in plan view.
%   cos(theta_h) = sin(theta_h) = cos(45) = 0.707.
%
% FAIRLEAD DISTANCE (R_fairlead):
%   Thesis Table 3.11: connection points at (±1.28, ±1.27) → r ≈ 1.80 m.
%   Verified 2026-05-03 against the project report.
%   Directly scales K_moor_yaw proportionally.
theta_v = deg2rad(45);  % Fairlead vertical angle from horizontal [rad]
theta_h = deg2rad(45);  % Fairlead plan-view angle from x-axis [rad]
R_fairlead = 1.802;     % Fairlead distance from platform center [m] (thesis Table 3.11)

% We look at the change in tension (Delta T) to isolate the stiffness
dT1 = T1 - mean(T1);
dT2 = T2 - mean(T2);
dT3 = T3 - mean(T3);
dT4 = T4 - mean(T4);

% Net Surge Restoring Force (Fx)
% Each line: Fx = T * cos(theta_v) * cos(theta_h), with appropriate sign
Fx_net = (dT1 + dT2 - dT3 - dT4) * cos(theta_v) * cos(theta_h);

% Net Sway Restoring Force (Fy)
Fy_net = (dT1 + dT4 - dT2 - dT3) * cos(theta_v) * sin(theta_h);

% Net Yaw Restoring Moment (Mz)
% Mz = F_tangential * R_fairlead, where F_tangential is the component
% of horizontal tension perpendicular to the radial line
Mz_net = (dT1 + dT3 - dT2 - dT4) * cos(theta_v) * R_fairlead;

%% 3. PERFORM LINEAR REGRESSION (F = K * x)
% The slope of the line is the Stiffness (K)
p_surge = polyfit(surge - mean(surge), Fx_net, 1);
K_moor_surge = abs(p_surge(1));

p_sway = polyfit(sway - mean(sway), Fy_net, 1);
K_moor_sway = abs(p_sway(1));

p_yaw = polyfit(yaw - mean(yaw), Mz_net, 1);
K_moor_yaw = abs(p_yaw(1));

%% 4. DISPLAY RESULTS
fprintf('\n--- DATA-DRIVEN MOORING STIFFNESS (K_moor) ---\n');
fprintf('Surge Mooring Stiffness: %.4f N/m\n', K_moor_surge);
fprintf('Sway  Mooring Stiffness: %.4f N/m\n', K_moor_sway);
fprintf('Yaw   Mooring Stiffness: %.4f Nm/rad\n', K_moor_yaw);
fprintf('----------------------------------------------\n');

% Save for downstream scripts
save('../outputs/sysid/MooringStiffness.mat', 'K_moor_surge', 'K_moor_sway', 'K_moor_yaw');
fprintf('Mooring stiffness saved to ../outputs/sysid/MooringStiffness.mat\n');

%% 5. VISUAL VERIFICATION
figure('Name', 'Mooring Stiffness Verification', 'Color', 'w', 'Position', [100 100 900 400]);

subplot(1,2,1);
scatter(surge - mean(surge), Fx_net, 10, 'b', 'filled'); hold on;
plot(surge - mean(surge), polyval(p_surge, surge - mean(surge)), 'r', 'LineWidth', 2);
title(sprintf('Surge Stiffness (K = %.2f N/m)', K_moor_surge));
xlabel('Displacement [m]'); ylabel('Restoring Force [N]');
grid on;
set(gca, 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
t1 = get(gca, 'Title'); set(t1, 'Color', 'k');

subplot(1,2,2);
scatter(sway - mean(sway), Fy_net, 10, 'g', 'filled'); hold on;
plot(sway - mean(sway), polyval(p_sway, sway - mean(sway)), 'r', 'LineWidth', 2);
title(sprintf('Sway Stiffness (K = %.2f N/m)', K_moor_sway));
xlabel('Displacement [m]'); ylabel('Restoring Force [N]');
grid on;
set(gca, 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
t2 = get(gca, 'Title'); set(t2, 'Color', 'k');