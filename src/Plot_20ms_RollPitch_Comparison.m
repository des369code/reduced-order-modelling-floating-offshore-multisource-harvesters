% Plot_20ms_RollPitch_Comparison.m
% Compares Roll and Pitch response at 20 m/s extreme wind
% Canonical platform vs Combined variant (ballast + wide)

clear; clc; close all;

% Load both CSVs
canon = readtable('../outputs/results/Final_Simulation_Results_20ms_canonical.csv');
comb  = readtable('../outputs/results/Final_Simulation_Results_20ms_combined.csv');

% Time vectors
t_canon = canon.Time;
t_comb  = comb.Time;

% Roll and Pitch
roll_canon = canon.Roll;
roll_comb  = comb.Roll;
pitch_canon = canon.Pitch;
pitch_comb  = comb.Pitch;

% Brighter, more visible colours
blue = [0.00, 0.45, 0.74];   % MATLAB default blue
red  = [0.85, 0.33, 0.10];   % MATLAB default red

% Create figure with forced white everywhere
figure('Position', [100, 100, 900, 600], 'Color', 'w', ...
       'InvertHardcopy', 'off');

% --- Roll ---
subplot(2,1,1);
hold on;
h1 = plot(t_canon, roll_canon, '-',  'Color', blue, 'LineWidth', 1.8);
h2 = plot(t_comb,  roll_comb,  '-',  'Color', red,  'LineWidth', 1.8);
hold off;
ylabel('Roll (rad)', 'FontSize', 12, 'Color', 'k');
title('Roll Response at 20 m/s Extreme Wind — Canonical vs Combined', ...
      'FontSize', 13, 'FontWeight', 'bold', 'Color', 'k');
h_leg1 = legend([h1, h2], 'Canonical', 'Combined', 'Location', 'best');
set(h_leg1, 'TextColor', 'k', 'Color', 'w');
grid on;
set(gca, 'Color', 'w', 'FontSize', 11, 'XColor', 'k', 'YColor', 'k', ...
         'GridColor', [0.15, 0.15, 0.15], 'GridAlpha', 0.15);

% --- Pitch ---
subplot(2,1,2);
hold on;
h3 = plot(t_canon, pitch_canon, '-',  'Color', blue, 'LineWidth', 1.8);
h4 = plot(t_comb,  pitch_comb,  '-',  'Color', red,  'LineWidth', 1.8);
hold off;
xlabel('Time (s)', 'FontSize', 12, 'Color', 'k');
ylabel('Pitch (rad)', 'FontSize', 12, 'Color', 'k');
title('Pitch Response at 20 m/s Extreme Wind — Canonical vs Combined', ...
      'FontSize', 13, 'FontWeight', 'bold', 'Color', 'k');
h_leg2 = legend([h3, h4], 'Canonical', 'Combined', 'Location', 'best');
set(h_leg2, 'TextColor', 'k', 'Color', 'w');
grid on;
set(gca, 'Color', 'w', 'FontSize', 11, 'XColor', 'k', 'YColor', 'k', ...
         'GridColor', [0.15, 0.15, 0.15], 'GridAlpha', 0.15);

% Print stats
fprintf('=== Roll Stats (20 m/s) ===\n');
fprintf('Canonical: max=%.4f rad (%.2f deg), RMS=%.4f rad (%.2f deg)\n', ...
    max(abs(roll_canon)), max(abs(roll_canon))*180/pi, rms(roll_canon), rms(roll_canon)*180/pi);
fprintf('Combined : max=%.4f rad (%.2f deg), RMS=%.4f rad (%.2f deg)\n', ...
    max(abs(roll_comb)), max(abs(roll_comb))*180/pi, rms(roll_comb), rms(roll_comb)*180/pi);

fprintf('\n=== Pitch Stats (20 m/s) ===\n');
fprintf('Canonical: max=%.4f rad (%.2f deg), RMS=%.4f rad (%.2f deg)\n', ...
    max(abs(pitch_canon)), max(abs(pitch_canon))*180/pi, rms(pitch_canon), rms(pitch_canon)*180/pi);
fprintf('Combined : max=%.4f rad (%.2f deg), RMS=%.4f rad (%.2f deg)\n', ...
    max(abs(pitch_comb)), max(abs(pitch_comb))*180/pi, rms(pitch_comb), rms(pitch_comb)*180/pi);

fprintf('\n=== Reduction (Combined vs Canonical) ===\n');
fprintf('Roll  max reduction: %.1f%%\n', (1 - max(abs(roll_comb))/max(abs(roll_canon)))*100);
fprintf('Pitch max reduction: %.1f%%\n', (1 - max(abs(pitch_comb))/max(abs(pitch_canon)))*100);

% Save figure with white background
exportgraphics(gcf, '../outputs/figures/RollPitch_20ms_Canonical_vs_Combined.png', ...
    'Resolution', 200, 'BackgroundColor', 'w');
fprintf('\nFigure saved to ../outputs/figures/RollPitch_20ms_Canonical_vs_Combined.png\n');
