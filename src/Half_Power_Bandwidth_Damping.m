% =========================================================================
% PHASE A: HALF-POWER BANDWIDTH DAMPING EXTRACTION (FIXED)
% =========================================================================
% Fixes applied:
%  - More Welch segments (window ~N/5) for reliable PSD estimates
%  - Direct crossing-point search (replaces unique() which destroyed
%    frequency ordering by sorting PSD values)
%  - Explicit warnings when fallback zeta is used
%  - Note: half-power bandwidth on forced-response data tends to
%    overestimate damping. The extracted zeta should be treated as an
%    UPPER BOUND, not the true physical damping.
%  - Saves results to DampingRatios.mat for downstream scripts.
% =========================================================================
clear; clc; close all;

%% 1. LOAD DATA & SETUP PSD
if ~isfile('../outputs/clean/Clean_6DOF_Position.csv')
    error('Clean_6DOF_Position.csv not found. Run Clean_Dataset_Script.m first.');
end
Table_Pos = readtable('../outputs/clean/Clean_6DOF_Position.csv');
t = Table_Pos.Time;
dt = t(2) - t(1);
Fs = 1/dt;
N_total = length(t);

motions = [Table_Pos.Surge, Table_Pos.Sway, Table_Pos.Heave, ...
           Table_Pos.Roll, Table_Pos.Pitch, Table_Pos.Yaw];
dof_names = {'Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'};

% Improved PSD parameters (matching PSD_NaturalFrequency.m)
window   = floor(N_total / 5);
noverlap = floor(window / 2);
nfft     = 2^nextpow2(4 * N_total);

fprintf('========================================================\n');
fprintf(' IDENTIFIED DAMPING RATIOS (HALF-POWER BANDWIDTH — FIXED)\n');
fprintf('========================================================\n');
fprintf(' NOTE: These are UPPER BOUND estimates. Forced response to\n');
fprintf(' narrowband JONSWAP waves broadens PSD peaks, inflating zeta.\n');
fprintf('========================================================\n');
fprintf('DOF \t wn (rad/s) \t w1 \t\t w2 \t\t Zeta (Damping)\n');
fprintf('--------------------------------------------------------\n');

zeta_values = zeros(6, 1);
wn_values   = zeros(6, 1);
fallback_used = false(6, 1);
zeta_capped   = false(6, 1);

%% 2. CALCULATE ZETA FOR ALL 6 DOFS
for i = 1:6
    % Compute PSD
    [Pxx, f] = pwelch(motions(:,i) - mean(motions(:,i)), ...
                      window, noverlap, nfft, Fs);
    omega = f * 2 * pi;

    % Find peak
    [P_max, idx_peak] = max(Pxx);
    wn = omega(idx_peak);
    wn_values(i) = wn;
    P_half = P_max / 2;

    % --- Find w1 (Left of Peak): direct crossing search ---
    w1 = NaN;
    for j = idx_peak:-1:2
        if (Pxx(j-1) - P_half) * (Pxx(j) - P_half) <= 0
            % Linear interpolation at the crossing
            w1 = interp1([Pxx(j-1), Pxx(j)], [omega(j-1), omega(j)], P_half, 'linear');
            break;
        end
    end
    if isnan(w1)
        w1 = wn * 0.90;  % Wider fallback than before (was 0.95)
        fallback_used(i) = true;
        fprintf('%-5s WARNING: Left half-power crossing not found. ', dof_names{i});
        fprintf('Using w1=0.90*wn (fallback).\n');
    end

    % --- Find w2 (Right of Peak): direct crossing search ---
    w2 = NaN;
    for j = idx_peak:length(omega)-1
        if (Pxx(j) - P_half) * (Pxx(j+1) - P_half) <= 0
            w2 = interp1([Pxx(j), Pxx(j+1)], [omega(j), omega(j+1)], P_half, 'linear');
            break;
        end
    end
    if isnan(w2)
        w2 = wn * 1.10;  % Wider fallback than before (was 1.05)
        fallback_used(i) = true;
        fprintf('%-5s WARNING: Right half-power crossing not found. ', dof_names{i});
        fprintf('Using w2=1.10*wn (fallback).\n');
    end

    % --- Calculate Raw Damping Ratio ---
    zeta_raw = (w2 - w1) / (2 * wn);

    % Sanity check: must be positive
    if zeta_raw <= 0
        fprintf('%-5s ERROR: Negative damping (zeta=%.4f). Setting to NaN.\n', ...
            dof_names{i}, zeta_raw);
        zeta_raw = NaN;
    end

    % --- PHYSICAL CAP: Half-power bandwidth on forced-response data ---
    % The half-power bandwidth method assumes free decay with white-noise
    % excitation (Bendat & Piersol, 2011). When applied to JONSWAP-driven
    % forced response, the spectral peak width reflects convolution of
    % structural damping WITH the wave excitation bandwidth, not damping
    % alone (Chakrabarti, 1987, Sec. 8.4). The result is an UPPER BOUND
    % that can exceed physically possible values.
    %
    % For a floating pontoon platform without dedicated damping devices
    % (bilge keels, heave plates, etc.), typical damping ratios are:
    %
    %   DOF       Typical zeta    Source
    %   Surge     0.05–0.15       DNV-RP-C205 Sec 6.3.2 (mooring drag adds ~0.05-0.10)
    %   Sway      0.05–0.15       DNV-RP-C205 Sec 6.3.2
    %   Heave     0.03–0.08       Faltinsen (1990) Sec 6.2; radiation + viscous on pontoons
    %   Roll      0.03–0.08       DNV-RP-C205; limited by pontoon waterplane
    %   Pitch     0.03–0.08       DNV-RP-C205; same mechanism as roll
    %   Yaw       0.05–0.15       DNV-RP-C205; mooring friction + wave drift
    %
    % Caps applied to prevent non-physical overdamping in the coupled
    % simulation. Surge/Sway are trusted up to zeta=0.25 (mooring drag
    % provides genuine hydrodynamic damping beyond the structural baseline).
    % Heave/Roll/Pitch/Yaw have no such mechanism and are capped lower.
    %
    % References:
    %   [1] DNV-RP-C205 (2010), "Environmental Conditions & Loads," Sec 6.3.
    %   [2] Faltinsen, O.M. (1990), "Sea Loads on Ships & Offshore Structures."
    %   [3] Chakrabarti, S.K. (1987), "Hydrodynamics of Offshore Structures."
    %   [4] Bendat, J.S. & Piersol, A.G. (2011), "Random Data."

    switch dof_names{i}
        case {'Surge', 'Sway'}
            zeta_max = 0.25;  % Mooring drag provides real damping
            zeta_lit = 0.15;  % Conservative literature value for cap
        case {'Heave', 'Roll', 'Pitch'}
            zeta_max = 0.30;  % Above this: definitively non-physical
            zeta_lit = 0.05;  % Faltinsen (1990); DNV-RP-C205
        case {'Yaw'}
            zeta_max = 0.30;
            zeta_lit = 0.10;  % Mooring friction + wave drift
    end

    if zeta_raw > zeta_max || isnan(zeta_raw)
        zeta = zeta_lit;
        fprintf('%-5s CAPPED: zeta_raw=%.4f > %.2f -> using literature zeta=%.2f [%s]\n', ...
            dof_names{i}, zeta_raw, zeta_max, zeta_lit, ...
            'DNV-RP-C205; Faltinsen (1990)');
        zeta_capped(i) = true;
    else
        zeta = zeta_raw;
        zeta_capped(i) = false;
    end

    zeta_values(i) = zeta;

    fprintf('%s \t %.4f \t %.4f \t %.4f \t %.4f (%.1f%%)', ...
        dof_names{i}, wn, w1, w2, zeta_raw, zeta_raw*100);
    if zeta_capped(i)
        fprintf(' -> capped to %.4f (%.1f%%)', zeta, zeta*100);
    end
    if fallback_used(i); fprintf(' [XING FALLBACK]'); end
    fprintf('\n');
end
fprintf('========================================================\n');

n_capped = sum(zeta_capped);
if n_capped > 0
    fprintf('PHYSICS CAP: %d DOF(s) capped to literature damping values.\n', n_capped);
    fprintf('  Rationale: Half-power bandwidth on forced JONSWAP response\n');
    fprintf('  overestimates zeta. See DNV-RP-C205 Sec 6.3.2, Faltinsen (1990).\n');
    fprintf('  Capped DOFs: ');
    for i = 1:6
        if zeta_capped(i)
            fprintf('%s (%.1f%%) ', dof_names{i}, zeta_values(i)*100);
        end
    end
    fprintf('\n');
end
if any(fallback_used)
    fprintf('WARNING: %d DOF(s) used fallback crossing frequencies.\n', sum(fallback_used));
end

%% 3. SAVE RESULTS
save('../outputs/sysid/DampingRatios.mat', 'zeta_values', 'wn_values', 'dof_names', 'fallback_used', 'zeta_capped');
fprintf('\nDamping ratios saved to ../outputs/sysid/DampingRatios.mat\n');
