% =========================================================================
% FIXED VITERNA-CORRIGAN EXTRAPOLATION
% No singularities at 180 degrees
% =========================================================================
clear; clc; close all;

%% 1. SETTINGS
filename = '../../data/airfoil/xf-n0009sm-il-1000000.csv'; % Source file
AR = 18.81;

%% 2. PROCESS
raw = readmatrix(filename, 'NumHeaderLines', 1);
[alpha_raw, idx] = unique(raw(:,1));
Cl_raw = raw(idx, 2);
Cd_raw = raw(idx, 3);

% Identify Stall (Max Lift)
[Cl_max, idx_stall] = max(Cl_raw);
alpha_stall = alpha_raw(idx_stall);
as_rad = deg2rad(alpha_stall);
Cd_stall = Cd_raw(idx_stall);

% Viterna Parameters
Cd_max = 1.11 + 0.018 * AR;
A1 = Cd_max / 2;
B1 = Cd_max;
A2 = (Cl_max - A1*sin(2*as_rad)) * sin(as_rad) / cos(as_rad)^2;
B2 = (Cd_stall - B1*sin(as_rad)^2) / cos(as_rad);

% Generate 360 Data
alpha_full = (-180:1:180)';
Cl_full = zeros(size(alpha_full));
Cd_full = zeros(size(alpha_full));

for i = 1:length(alpha_full)
   a = alpha_full(i);
   a_abs = abs(a);
   a_rad = deg2rad(a_abs);
  
   if a_abs <= alpha_stall
       % Linear Region
       Cl_full(i) = interp1(alpha_raw, Cl_raw, a, 'linear', 'extrap');
       Cd_full(i) = interp1(alpha_raw, Cd_raw, a, 'linear', 'extrap');
   elseif a_abs <= 90
       % Post-Stall (Viterna)
       Cl = A1*sin(2*a_rad) + A2*cos(a_rad)^2/sin(a_rad);
       Cd = B1*sin(a_rad)^2 + B2*cos(a_rad);
       Cl_full(i) = sign(a) * Cl;
       Cd_full(i) = Cd;
   else
       % Backside (90-180) - FIXED FORMULA
       % Simple Flat Plate model to avoid singularity
       rem_rad = deg2rad(180 - a_abs);
       Cl = 1.1 * sin(2 * rem_rad); % Smooth sine to 0 at 180
       Cd = Cd_max * sin(rem_rad)^2 + Cd_raw(1)*cos(rem_rad)^2;
      
       % Invert lift sign for backside
       Cl_full(i) = -sign(a) * Cl;
       Cd_full(i) = Cd;
   end
end

% Save
writematrix([alpha_full,Cl_full, Cd_full], '../../data/airfoil/NACA0009_360_Full_New2.csv');
fprintf('Success! Created ../../data/airfoil/NACA0009_360_Full_New2.csv\n');