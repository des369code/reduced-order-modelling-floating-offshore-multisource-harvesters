% =========================================================================
% PHASE A: 6-DOF POWER SPECTRAL DENSITY (PSD) & NATURAL FREQUENCIES
% =========================================================================
% Improved PSD estimation: more Welch segments (window ~N/5 = 100 samples,
% ~9 segments) vs the old default of 3 segments (window ~N/2).
% Uses findpeaks() with prominence threshold instead of naive max().
% Uses parabolic peak fitting for sub-bin frequency refinement.
% Saves results to NaturalFrequencies.mat for downstream scripts.
% =========================================================================
clear; clc; close all;

%% 1. LOAD CLEAN DATA
fprintf('Loading Clean Checkpoint Data...\n');
if ~isfile('../outputs/clean/Clean_6DOF_Position.csv')
    error('Clean_6DOF_Position.csv not found. Run Clean_Dataset_Script.m first.');
end
Table_Pos = readtable('../outputs/clean/Clean_6DOF_Position.csv');
t = Table_Pos.Time;
dt = t(2) - t(1);
Fs = 1/dt;
N_total = length(t);

% Extract all 6 DOFs (Roll, Pitch, Yaw are already in radians)
motions = [Table_Pos.Surge, Table_Pos.Sway, Table_Pos.Heave, ...
           Table_Pos.Roll, Table_Pos.Pitch, Table_Pos.Yaw];
dof_names = {'Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'};
units = {'m^2/s', 'm^2/s', 'm^2/s', 'rad^2/s', 'rad^2/s', 'rad^2/s'};

%% 2. PSD PARAMETERS (Improved: ~9 Welch segments for reliable estimates)
window   = floor(N_total / 5);   % ~100 samples = 10s window
noverlap = floor(window / 2);    % 50% overlap
nfft     = 2^nextpow2(4 * N_total);  % ~2048 points for fine frequency grid
df       = Fs / nfft;            % Frequency resolution

fprintf('PSD settings: N=%d, window=%d, noverlap=%d, nfft=%d, df=%.4f Hz\n', ...
    N_total, window, noverlap, nfft, df);
fprintf('Number of Welch segments: ~%d (vs 3 previously)\n\n', ...
    floor((N_total - window) / (window - noverlap)) + 1);

wn_values = zeros(6, 1);

figure('Name', '6-DOF Power Spectral Density (Improved)', ...
       'Color', 'w', 'Position', [50 50 1400 900]);

for i = 1:6
    % Compute PSD
    current_motion = motions(:, i);
    [Pxx, f] = pwelch(current_motion - mean(current_motion), ...
                      window, noverlap, nfft, Fs);
    omega = f * 2 * pi; % Convert Hz to rad/s

    % Find dominant peak using findpeaks with prominence threshold
    min_prominence = max(Pxx) * 0.1;  % 10% of max PSD as minimum prominence
    [pks, locs, w_prom] = findpeaks(Pxx, 'MinPeakProminence', min_prominence, ...
                                    'SortStr', 'descend', 'NPeaks', 3);

    if isempty(pks)
        % Fallback to max() if findpeaks finds nothing
        [max_Pxx, idx_peak] = max(Pxx);
        wn = omega(idx_peak);
        fprintf('%-5s WARNING: findpeaks found no peaks; using max() fallback\n', dof_names{i});
    else
        % Parabolic fit around the highest peak for sub-bin refinement
        idx_peak = locs(1);
        if idx_peak > 1 && idx_peak < length(omega)
            % Fit parabola to 3 points around peak: y = a*f^2 + b*f + c
            f_vals = f(idx_peak-1 : idx_peak+1);
            p_vals = Pxx(idx_peak-1 : idx_peak+1);
            p_fit = polyfit(f_vals, p_vals, 2);
            f_peak = -p_fit(2) / (2 * p_fit(1));  % Vertex of parabola
            wn = f_peak * 2 * pi;
        else
            wn = omega(idx_peak);
        end
        max_Pxx = pks(1);

        % Report prominence of dominant peak
        if length(pks) > 1
            fprintf('%-5s Peak 1: %.4f rad/s (prominence=%.1f%%), Peak 2: %.4f rad/s\n', ...
                dof_names{i}, wn, 100*w_prom(1)/max(Pxx), omega(locs(2)));
        else
            fprintf('%-5s Peak 1: %.4f rad/s (prominence=%.1f%%)\n', ...
                dof_names{i}, wn, 100*w_prom(1)/max(Pxx));
        end
    end

    wn_values(i) = wn;
    fprintf('%-5s Natural Frequency : %.4f rad/s  (%.4f Hz)\n', ...
        dof_names{i}, wn, wn/(2*pi));

    % Plot
    subplot(2, 3, i);
    plot(omega, Pxx, 'b', 'LineWidth', 1.5); hold on;
    plot(wn, max_Pxx, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    % Also mark secondary peaks
    if length(locs) > 1
        for k = 2:min(length(locs), 3)
            plot(omega(locs(k)), pks(k), 'kv', 'MarkerSize', 6);
        end
    end
    title(sprintf('%s (\\omega_n = %.3f rad/s)', dof_names{i}, wn));
    xlabel('Frequency \omega (rad/s)');
    ylabel(sprintf('Energy (%s)', units{i}));
    grid on;
    xlim([0 3]);
    set(gca, 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
    title_handle = get(gca, 'Title');
    set(title_handle, 'Color', 'k');
end

h_sg = sgtitle('\bf 6-DOF Power Spectral Density — Improved Welch Estimate', 'FontSize', 14);
set(h_sg, 'Color', 'k');

%% 3. SAVE RESULTS
save('../outputs/sysid/NaturalFrequencies.mat', 'wn_values', 'dof_names');
fprintf('\nNatural frequencies saved to ../outputs/sysid/NaturalFrequencies.mat\n');
fprintf('Verification: sqrt(K_moor/M_total) should be close to wn for surge/sway.\n');
fprintf('  If they differ significantly (>20%%), investigate PSD peaks or mooring K.\n');
