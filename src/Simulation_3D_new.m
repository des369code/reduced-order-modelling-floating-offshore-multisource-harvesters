function Simulation_3D(W, L, H_mast, csv_filename)
    % =====================================================================
    % 3D DYNAMIC STABILITY VISUALIZER (SMART WIND & WAVE HUD)
    % =====================================================================
    
    % Set defaults if no arguments are provided
    if nargin < 4, csv_filename = '../outputs/results/Final_Simulation_Results.csv'; end
    if nargin < 3, H_mast = 5; end
    if nargin < 2, L = 4; end
    if nargin < 1, W = 4; end

    clearvars -except W L H_mast csv_filename; clc; close all;
    
    %% 1. LOAD DATA
    if ~isfile(csv_filename)
        error('Cannot find %s. Please run the simulation first.', csv_filename);
    end
    Res = readtable(csv_filename);
    
    % Extract Data Vectors
    t = Res.Time;
    x = Res.Surge;  y = Res.Sway;  z = Res.Heave;
    rx = deg2rad(Res.Roll);  
    ry = deg2rad(Res.Pitch); 
    rz = deg2rad(Res.Yaw);
    
    %% 2. GEOMETRY DEFINITION
    platform_points = [-L/2 -W/2 0; L/2 -W/2 0; L/2 W/2 0; -L/2 W/2 0; 0 0 0; 0 0 H_mast]; 
    
    %% 3. SETUP FIGURE & UI
    fig = figure('Name', '3D Stability Analysis (Interactive)', 'Color', 'w', 'Position', [100 100 1000 700]);
    ax = axes('Parent', fig, 'Units', 'pixels', 'Position', [100 150 800 500]);
    view(3); grid on; axis equal; hold on;
    
    % Create Sea Surface
    [X_sea, Y_sea] = meshgrid(-20:5:60, -20:5:20);
    surf(ax, X_sea, Y_sea, zeros(size(X_sea)), 'FaceAlpha', 0.1, 'EdgeColor', [0.7 0.7 1], 'FaceColor', 'b');
    
    % Initialize Platform Graphics
    h_beams = plot3(ax, 0,0,0, 'Color', [0.3 0.3 0.3], 'LineWidth', 3); 
    h_mast = plot3(ax, 0,0,0, 'k', 'LineWidth', 5);
    h_pontoons = scatter3(ax, 0,0,0, 300, 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
    h_blades = plot3(ax, 0,0,0, 'k', 'LineWidth', 2.5); % Turbine blades
    
    xlabel('Surge (X) [m]'); ylabel('Sway (Y) [m]'); zlabel('Heave (Z) [m]');
    
    %% 3b. SMART ENVIRONMENTAL HUD
    % --- Detect WIND Direction ---
    if isfile('../outputs/forces/Clean_Aero_Force.csv')
        Aero = readtable('../outputs/forces/Clean_Aero_Force.csv');
        theta_wind = atan2(mean(Aero.F_Sway), mean(Aero.F_Surge)); 
    else
        theta_wind = deg2rad(135); 
    end
    theta_wind_deg = mod(round(rad2deg(theta_wind)), 360); 
    
    % --- Detect WAVE Direction ---
    if isfile('../outputs/forces/Clean_Wave_Force.csv')
        Wave = readtable('../outputs/forces/Clean_Wave_Force.csv');
        force_mag = sqrt(Wave.F_Surge.^2 + Wave.F_Sway.^2);
        [~, peak_idx] = max(force_mag);
        theta_wave = atan2(Wave.F_Sway(peak_idx), Wave.F_Surge(peak_idx));
        
        dot_product = cos(theta_wind)*cos(theta_wave) + sin(theta_wind)*sin(theta_wave);
        if dot_product < 0
            theta_wave = theta_wave + pi; 
        end
    else
        theta_wave = deg2rad(135);
    end
    theta_wave_deg = mod(round(rad2deg(theta_wave)), 360); 
    
    % --- Draw HUD Objects ---
    h_wind = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color', 'c', 'LineWidth', 3, 'MaxHeadSize', 1.5);
    h_text_wind = text(ax, 0,0,0, sprintf(' WIND (%d^\\circ)', theta_wind_deg), 'Color', 'c', 'FontWeight', 'bold', 'FontSize', 10);
    
    h_wave = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color', 'm', 'LineWidth', 4, 'MaxHeadSize', 1.5);
    h_text_wave = text(ax, 0,0,0, sprintf(' WAVE (%d^\\circ)', theta_wave_deg), 'Color', 'm', 'FontWeight', 'bold', 'FontSize', 10);
    
    h_ref = quiver3(ax, 0,0,0, 0,0,0, 0, 'Color', 'w', 'LineStyle', '-', 'LineWidth', 2, 'MaxHeadSize', 1.0);
    h_text_ref = text(ax, 0,0,0, ' 0^\circ (+X Axis)', 'Color', 'w', 'FontWeight', 'bold', 'FontSize', 9);
    
    %% 4. SLIDER UI
    uicontrol('Style', 'text', 'Position', [400 60 200 20], 'String', 'Drag to Scrub Through Time', 'BackgroundColor', 'w');
    
    sld = uicontrol('Style', 'slider', 'Min', 1, 'Max', length(t), 'Value', 1, ...
        'Position', [200 40 600 20], 'Callback', @(src, event) updatePlot(round(src.Value), rx, ry, rz, x, y, z, platform_points, h_beams, h_mast, h_pontoons, h_blades, h_wind, h_wave, h_ref, h_text_wind, h_text_wave, h_text_ref, theta_wind, theta_wave, H_mast, ax, t));
    
    updatePlot(1, rx, ry, rz, x, y, z, platform_points, h_beams, h_mast, h_pontoons, h_blades, h_wind, h_wave, h_ref, h_text_wind, h_text_wave, h_text_ref, theta_wind, theta_wave, H_mast, ax, t);
end 

%% 5. LOCAL UPDATE FUNCTION
function updatePlot(idx, rx, ry, rz, x, y, z, platform_points, h_beams, h_mast, h_pontoons, h_blades, h_wind, h_wave, h_ref, h_text_wind, h_text_wave, h_text_ref, theta_wind, theta_wave, H_mast, ax, t)
    
    % --- PLATFORM KINEMATICS ---
    Rx = [1 0 0; 0 cos(rx(idx)) -sin(rx(idx)); 0 sin(rx(idx)) cos(rx(idx))];
    Ry = [cos(ry(idx)) 0 sin(ry(idx)); 0 1 0; -sin(ry(idx)) 0 cos(ry(idx))];
    Rz = [cos(rz(idx)) -sin(rz(idx)) 0; sin(rz(idx)) cos(rz(idx)) 0; 0 0 1];
    R = Rz * Ry * Rx;
    
    pts = (R * platform_points')';
    pts(:,1) = pts(:,1) + x(idx); 
    pts(:,2) = pts(:,2) + y(idx); 
    pts(:,3) = pts(:,3) + z(idx); 
    
    base_indices = [1 2 3 4 1];
    set(h_beams, 'XData', pts(base_indices, 1), 'YData', pts(base_indices, 2), 'ZData', pts(base_indices, 3));
    set(h_mast, 'XData', pts(5:6, 1), 'YData', pts(5:6, 2), 'ZData', pts(5:6, 3));
    set(h_pontoons, 'XData', pts(1:4, 1), 'YData', pts(1:4, 2), 'ZData', pts(1:4, 3));
    
    % --- TURBINE SPINNING ---
    R_rotor = 2.5; omega = 4; theta_rotor = omega * t(idx);
    a1 = theta_rotor; a2 = theta_rotor + (2*pi/3); a3 = theta_rotor + (4*pi/3);
    hub = [0, 0, 0];
    tip1 = [0, R_rotor*cos(a1), R_rotor*sin(a1)];
    tip2 = [0, R_rotor*cos(a2), R_rotor*sin(a2)];
    tip3 = [0, R_rotor*cos(a3), R_rotor*sin(a3)];
    local_blades = [hub; tip1; hub; tip2; hub; tip3];
    
    R_yaw_wind = [cos(theta_wind) -sin(theta_wind) 0; sin(theta_wind) cos(theta_wind) 0; 0 0 1];
    aligned_blades = (R_yaw_wind * local_blades')';
    aligned_blades(:,3) = aligned_blades(:,3) + H_mast;
    
    pts_blades = (R * aligned_blades')';
    pts_blades(:,1) = pts_blades(:,1) + x(idx); 
    pts_blades(:,2) = pts_blades(:,2) + y(idx); 
    pts_blades(:,3) = pts_blades(:,3) + z(idx); 
    set(h_blades, 'XData', pts_blades(:,1), 'YData', pts_blades(:,2), 'ZData', pts_blades(:,3));

    % --- HUD ARROWS ---
    arrow_len = 5; 
    hud_x = x(idx) - 8; hud_y = y(idx) - 8;
    
    set(h_wind, 'XData', hud_x, 'YData', hud_y, 'ZData', 8, 'UData', arrow_len*cos(theta_wind), 'VData', arrow_len*sin(theta_wind), 'WData', 0);
    set(h_text_wind, 'Position', [hud_x, hud_y, 9]);
    
    set(h_ref, 'XData', hud_x, 'YData', hud_y, 'ZData', 4, 'UData', arrow_len*cos(0), 'VData', arrow_len*sin(0), 'WData', 0);
    set(h_text_ref, 'Position', [hud_x, hud_y, 5]);

    set(h_wave, 'XData', hud_x, 'YData', hud_y, 'ZData', 0, 'UData', arrow_len*cos(theta_wave), 'VData', arrow_len*sin(theta_wave), 'WData', 0);
    set(h_text_wave, 'Position', [hud_x, hud_y, 1]);
    
    % --- CAMERA & TITLE ---
    xlim(ax, [x(idx)-12, x(idx)+12]); ylim(ax, [y(idx)-12, y(idx)+12]); zlim(ax, [-10, 10]);
    
    if abs(rad2deg(rx(idx))) > 90 || abs(rad2deg(ry(idx))) > 90
        status_color = 'r'; status_text = ' [CAPSIZED!]';
    else
        status_color = 'k'; status_text = '';
    end

    t_obj = title(ax, sprintf('TIME: %.2f s | SURGE: %.1f m | ROLL: %.1f^o | PITCH: %.1f^o%s', ...
        t(idx), x(idx), rad2deg(rx(idx)), rad2deg(ry(idx)), status_text));
    set(t_obj, 'Color', status_color);
end