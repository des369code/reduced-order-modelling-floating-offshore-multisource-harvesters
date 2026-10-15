% =========================================================================
% PHASE A6: FINAL MASTER SYSTEM MATRICES (WIDER PLATFORM TEST)
% =========================================================================
clear; clc;

%% 1. STRUCTURAL MASS & INERTIA
M_dry = 836 + 567;  % Platform + Turbine (kg)
Ixx_dry = 4042.45 + 6480.73;  % Total Roll Inertia (kg*m^2)
Iyy_dry = 7502.24 + 6487.92;  % Total Pitch Inertia (kg*m^2)
Izz_dry = 11473.02 + 12034.07; % Total Yaw Inertia (kg*m^2)

%% 2. ADDED MASS APPROXIMATIONS (DNV Guidelines)
M_add_surge = 0.20 * M_dry;
M_add_heave = 1.00 * M_dry;
Ixx_add = 0.20 * Ixx_dry;
Iyy_add = 0.20 * Iyy_dry;
Izz_add = 0.10 * Izz_dry;

M_total = [M_dry + M_add_surge;  % 1: Surge
           M_dry + M_add_surge;  % 2: Sway
           M_dry + M_add_heave;  % 3: Heave
           Ixx_dry + Ixx_add;    % 4: Roll
           Iyy_dry + Iyy_add;    % 5: Pitch
           Izz_dry + Izz_add];   % 6: Yaw

%% 3. DATA-DRIVEN INPUTS (Loaded from system identification scripts)
if isfile('../outputs/sysid/NaturalFrequencies.mat')
    load('../outputs/sysid/NaturalFrequencies.mat', 'wn_values');
    wn = wn_values;
    fprintf('Loaded natural frequencies from NaturalFrequencies.mat\n');
else
    warning('NaturalFrequencies.mat not found. Using cached values.');
    wn = [0.7363; 0.7363; 0.6136; 0.8590; 0.8590; 0.6136];
end

if isfile('../outputs/sysid/DampingRatios.mat')
    load('../outputs/sysid/DampingRatios.mat', 'zeta_values');
    zeta = zeta_values;
    fprintf('Loaded damping ratios from DampingRatios.mat\n');
else
    warning('DampingRatios.mat not found. Using cached values.');
    zeta = [0.2431; 0.2465; 0.7794; 0.7921; 0.8406; 0.3129];
end

%% 4. STIFFNESS ASSEMBLY (Hybrid Mooring + Hydrostatic)
if isfile('../outputs/sysid/MooringStiffness.mat')
    load('../outputs/sysid/MooringStiffness.mat', 'K_moor_surge', 'K_moor_sway', 'K_moor_yaw');
    fprintf('Loaded mooring stiffness from MooringStiffness.mat\n');
else
    warning('MooringStiffness.mat not found. Using cached values.');
    K_moor_surge = 510.6695;
    K_moor_sway  = 2374.9467;
    K_moor_yaw   = 9333.9936;
end

K_diag = zeros(6,1);
K_diag(1) = K_moor_surge;
K_diag(2) = K_moor_sway;
K_diag(3) = M_total(3) * (wn(3)^2); % Heave (Spectral)
K_diag(4) = M_total(4) * (wn(4)^2); % Roll (Spectral)
K_diag(5) = M_total(5) * (wn(5)^2); % Pitch (Spectral)
K_diag(6) = K_moor_yaw;

% Wider platform: exponentially higher roll/pitch hydrostatic stiffness
K_diag(4) = K_diag(4) * 10;
K_diag(5) = K_diag(5) * 10;

%% 5. DAMPING CALCULATION
C_diag = zeros(6,1);
for i = 1:6
    C_diag(i) = 2 * zeta(i) * M_total(i) * wn(i);
end

%% 6. ASSEMBLE AND SAVE
M_mat = diag(M_total); C_mat = diag(C_diag); K_mat = diag(K_diag);

%% 7. DISPLAY RESULTS
fprintf('=========================================================================\n');
fprintf(' WIDE PLATFORM VARIANT: 10x K_roll/K_pitch\n');
fprintf('=========================================================================\n');
dof_names = {'Surge', 'Sway ', 'Heave', 'Roll ', 'Pitch', 'Yaw  '};
fprintf('DOF \t Total Mass (M) \t Stiffness (K) \t\t Damping (C)\n');
fprintf('-------------------------------------------------------------------------\n');
for i = 1:6
    fprintf('%s \t %.2f \t\t %.2f \t\t %.2f\n', dof_names{i}, M_total(i), K_diag(i), C_diag(i));
end
fprintf('=========================================================================\n');
fprintf('* Units: M [kg, kg*m^2], K [N/m, Nm/rad], C [Ns/m, Nms/rad]\n');

save('../outputs/sysid/MasterSystemData.mat', 'M_mat', 'C_mat', 'K_mat');
fprintf('Wide platform variant saved to MasterSystemData.mat.\n');
