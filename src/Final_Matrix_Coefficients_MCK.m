% =========================================================================
% PHASE A6: FINAL MASTER SYSTEM MATRICES (M, C, K)
% =========================================================================
clear; clc;

%% 1. STRUCTURAL MASS & INERTIA (Tables 3.15 - 3.18)
M_dry = 836 + 567; % Platform + Turbine (kg)
Ixx_dry = 4042.45 + 6480.73; % Total Roll Inertia (kg*m^2)
Iyy_dry = 7502.24 + 6487.92; % Total Pitch Inertia (kg*m^2)
Izz_dry = 11473.02 + 12034.07; % Total Yaw Inertia (kg*m^2)

%% 2. ADDED MASS APPROXIMATIONS (DNV Guidelines)
% Applying hydrodynamic coefficients based on DNV-RP-C205
M_add_surge = 0.20 * M_dry;  
M_add_heave = 1.00 * M_dry;  
Ixx_add = 0.20 * Ixx_dry;    
Iyy_add = 0.20 * Iyy_dry;    
Izz_add = 0.10 * Izz_dry;    

% Final Total Mass Vector (M_total)
M_total = [M_dry + M_add_surge;  % 1: Surge
           M_dry + M_add_surge;  % 2: Sway
           M_dry + M_add_heave;  % 3: Heave
           Ixx_dry + Ixx_add;    % 4: Roll
           Iyy_dry + Iyy_add;    % 5: Pitch
           Izz_dry + Izz_add];   % 6: Yaw

%% 3. DATA-DRIVEN INPUTS (Loaded from system identification scripts)
% Load natural frequencies from PSD analysis
if isfile('../outputs/sysid/NaturalFrequencies.mat')
    load('../outputs/sysid/NaturalFrequencies.mat', 'wn_values');
    wn = wn_values;
    fprintf('Loaded natural frequencies from NaturalFrequencies.mat\n');
else
    warning('NaturalFrequencies.mat not found. Using cached values.');
    wn = [0.7363; 0.7363; 0.6136; 0.8590; 0.8590; 0.6136];
end

% Load damping ratios from half-power bandwidth analysis
if isfile('../outputs/sysid/DampingRatios.mat')
    load('../outputs/sysid/DampingRatios.mat', 'zeta_values');
    zeta = zeta_values;
    fprintf('Loaded damping ratios from DampingRatios.mat\n');
else
    warning('DampingRatios.mat not found. Using cached values.');
    zeta = [0.2431; 0.2465; 0.7794; 0.7921; 0.8406; 0.3129];
end

%% 4. STIFFNESS ASSEMBLY (Hybrid Mooring + Hydrostatic)
% Mooring stiffness from tension-displacement regression
if isfile('../outputs/sysid/MooringStiffness.mat')
    load('../outputs/sysid/MooringStiffness.mat', 'K_moor_surge', 'K_moor_sway', 'K_moor_yaw');
    fprintf('Loaded mooring stiffness from MooringStiffness.mat\n');
else
    warning('MooringStiffness.mat not found. Using cached values.');
    K_moor_surge = 510.6695;
    K_moor_sway  = 2374.9467;
    K_moor_yaw   = 9333.9936;
end

% For Heave, Roll, and Pitch, we use the values derived from 
% the natural frequency formula as they are dominated by Hydrostatics (Buoyancy).
K_diag = zeros(6,1);
K_diag(1) = K_moor_surge;
K_diag(2) = K_moor_sway;
K_diag(3) = M_total(3) * (wn(3)^2); % Heave (Spectral)
K_diag(4) = M_total(4) * (wn(4)^2); % Roll (Spectral)
K_diag(5) = M_total(5) * (wn(5)^2); % Pitch (Spectral)
K_diag(6) = K_moor_yaw;

%% 5. DAMPING CALCULATION (C)
C_diag = zeros(6,1);
for i = 1:6
    % C = 2 * zeta * M * wn
    C_diag(i) = 2 * zeta(i) * M_total(i) * wn(i);
end

%% 6. ASSEMBLE FINAL 6x6 MATRICES
M_mat = diag(M_total);
C_mat = diag(C_diag);
K_mat = diag(K_diag);

%% 7. DISPLAY RESULTS FOR FYP REPORT
fprintf('=========================================================================\n');
fprintf(' FINAL MASTER SYSTEM MATRICES (AERO-HYDRO-MOORING COUPLING)\n');
fprintf('=========================================================================\n');
fprintf('DOF \t Total Mass (M) \t Stiffness (K) \t\t Damping (C)\n');
fprintf('-------------------------------------------------------------------------\n');
dof_names = {'Surge', 'Sway ', 'Heave', 'Roll ', 'Pitch', 'Yaw  '};
for i = 1:6
    fprintf('%s \t %.2f \t\t %.2f \t\t %.2f\n', dof_names{i}, M_total(i), K_diag(i), C_diag(i));
end
fprintf('=========================================================================\n');
fprintf('* Units: M [kg, kg*m^2], K [N/m, Nm/rad], C [Ns/m, Nms/rad]\n');

%% 8. CONSISTENCY CHECK: Verify wn ≈ sqrt(K/M)
fprintf('\n--- CONSISTENCY CHECK: wn vs sqrt(K/M) ---\n');
fprintf('DOF \t PSD wn \t sqrt(K/M) \t Diff (%%)\n');
fprintf('-----------------------------------------\n');
dof_labels = {'Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'};
for i = 1:6
    wn_from_K = sqrt(K_diag(i) / M_total(i));
    diff_pct = 100 * abs(wn(i) - wn_from_K) / wn(i);
    fprintf('%s \t %.4f \t %.4f \t %.1f%%', dof_labels{i}, wn(i), wn_from_K, diff_pct);
    if diff_pct > 20
        fprintf(' *** INVESTIGATE');
    elseif diff_pct > 10
        fprintf(' *');
    end
    fprintf('\n');
end
fprintf('-----------------------------------------\n');
fprintf('>20%% difference: PSD peak may not be the true natural frequency,\n');
fprintf('  or mooring stiffness / mass may be incorrect.\n');

% Save for ODE Solver
save('../outputs/sysid/MasterSystemData.mat', 'M_mat', 'C_mat', 'K_mat');
fprintf('Successfully saved MasterSystemData.mat for final simulation.\n');