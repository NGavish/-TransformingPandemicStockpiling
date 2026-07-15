function Figure1_Japan_Scenario
% --- Setup and Data Initialization ---
close all;
figure('Color', 'w', 'Position', [100, 100, 1000, 500]);
t = tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Japan Population estimate (approximate for standardization)
N_Japan = 123.8e6; 

% Japan Scenario Targets and Calibrated Values (From Table 2)
beta = 0.72; 
gammaI = 0.5; 
gammaH = 0.25; 
clinicalAttackRate = 0.25; % 25% clinical attack rate
hosp_rate = 0.078;         % 7.8% Hospitalization rate of symptomatic cases (Standardized)
p=0.0405;
duration_weeks = 8;        % 8 weeks duration
duration_days = duration_weeks * 7;

% --- Panel A: Static Model (Lower Bound) ---
nexttile
lg = plotPanelA_Japan(N_Japan, clinicalAttackRate, hosp_rate, duration_days, gammaH);
text(0.05, 0.95, 'A', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');
box on;

% --- Panel B: Dynamic Model (SIHR) ---
nexttile
plotPanelB_Japan(N_Japan, beta, hosp_rate, gammaI, gammaH, clinicalAttackRate, duration_days, duration_weeks);
text(0.05, 0.95, 'B', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');
box on;

% --- Shared Legend ---
lg.Layout.Tile = 'South';
lg.Orientation = 'horizontal';

printGraph('./graphs/Figure_Japan_Scenario')

end

function lg = plotPanelA_Japan(N, clinicalAR, hospRate, duration_days, gammaH)
    % Calculate flat distribution (Static Model) for Japan
    total_clinical_cases = N * clinicalAR;
    total_hospitalizations = total_clinical_cases * hospRate;
    length_of_stay = 1 / gammaH;
    total_bed_days = total_hospitalizations * length_of_stay;
    
    val_unmit = total_bed_days / duration_days; % Flat daily demand
    
    % Antiviral mitigation bounds (Assuming 30% to 60% effectiveness)
    av_min = 0.4 * val_unmit; % 60% reduction
    av_max = 0.7 * val_unmit; % 30% reduction
    val_mit = (av_min + av_max) / 2;
    err_l = val_mit - av_min;
    err_h = av_max - val_mit;

    hold on;
    % Scaled Background Zones for Japan (proportional to ~60k bed capacity threshold)
    x_limits = [0 2 2 0 0]; 
    % p6 = fill(x_limits, [2e4 2e4 4e4 4e4 2e4], [1 1 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    % p7 = fill(x_limits, [4e4 4e4 6e4 6e4 4e4], [1 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    p_line = yline(180000, 'r--', 'LineWidth', 2);
    
    % Bars
    b1 = bar(0.85, val_unmit, 0.25);
    b2 = bar(1.15, val_mit, 0.25);
    
    % Style
    b1.FaceColor = [0.9, 0.9, 0.9]; b1.EdgeColor = [55 112 184]/255; b1.LineWidth = 1.5; b1.FaceAlpha = 0.7;
    b2.FaceColor = [0.4, 0.4, 0.4]; b2.EdgeColor = 'k'; b2.LineWidth = 1; b2.FaceAlpha = 0.7;

    % Errorbar
    errorbar(1.15, val_mit, err_l, err_h, 'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

    % Formatting
    ylabel('Daily Bed Demand', 'FontSize', 12);
    title('Static Model (Japan Lower Bound)', 'FontSize', 12);
    
    % Scale Y-axis to accommodate Japan's peak
    ylim([0 5e5]); 
    set(gca, 'YTick', 0:1e5:6e5);
    ytickformat('%.0fk');
    labels = arrayfun(@(v) sprintf('%.0fk', v/1e3), get(gca, 'YTick'), 'UniformOutput', false);
    set(gca, 'YTickLabel', labels);
    
    xlim([0.5 1.5]);
    set(gca, 'XTick', []); 
    grid on; ax = gca; ax.YGrid = 'on'; ax.XGrid = 'off';
    
    % Legend handle
    lg = legend([b1, b2,  p_line], ...
        {'Unmitigated', 'With Antivirals',  'Collapse Threshold (~180k)'}, ...
        'FontSize', 10);
    hold off;
end

function plotPanelB_Japan(N, beta, p, gammaI, gammaH, clinicalAttackRate, duration_days, duration_weeks)
    % Initial state
    I0 = 1e-3; 
    y0 = [1-I0; I0; 0; 0; 0];
    
    % Unmitigated Simulation
    [t, y] = ode45(@(t,y)SIHR(t,y,beta,gammaI,gammaH,p), [0 duration_days], y0);
    
    Z = y(end,5); % Cumulative transitions
    
    % Mitigated Simulation
    [t_mit, y_mit] = ode45(@(t,y)SIHR(t,y,beta,gammaI,gammaH,p*clinicalAttackRate/Z), [0 duration_days], y0);
    
    % Scale time to weeks
    t_weeks = t_mit/7; 
    H = y_mit(:,3); % Hospitalized compartment
    
    hold on;
    % Backgrounds (Scaled to Japan ~60k Collapse Threshold)
    x_lims = [0 duration_weeks duration_weeks 0 0];
    % fill(x_lims, [2e4 2e4 4e4 4e4 2e4], [1 1 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    % fill(x_lims, [4e4 4e4 6e4 6e4 4e4], [1 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    yline(180000, 'r--', 'LineWidth', 2);
    text(0.5, 182000, {'Free bed capacity','          (~180k)'}, 'Color', 'r', 'FontSize', 9, 'VerticalAlignment', 'middle');

    % Curves
    h_lower = N * H * 0.4; % Best case (60% efficacy)
    h_upper = N * H * 0.7; % Worst case (30% efficacy)
    fill([t_weeks; flipud(t_weeks)], [h_upper; flipud(h_lower)], [0.4 0.4 0.4], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    plot(t_weeks, N*H, 'Color', [55 112 184]/255, 'LineWidth', 2);

    % Formatting
    title('Dynamic Model (Japan Scenario)', 'FontSize', 12);
    ylabel('Daily Bed Demand');
    xlabel('Weeks into wave');
    
    % Scale Y-axis to accommodate Japan's peak
    ylim([0 5e5]);
    set(gca, 'YTick', 0:1e5:6e5);
    ytickformat('%.0fk');
    labels = arrayfun(@(v) sprintf('%.0fk', v/1e3), get(gca, 'YTick'), 'UniformOutput', false);
    set(gca, 'YTickLabel', labels);
    
    xlim([0 duration_weeks]);
    grid on;
end

% SIHR Differential Equations
function dy=SIHR(t,y,beta,gammaI,gammaH,pHospitalization)
    S=y(1); I=y(2); H=y(3); R=y(4); ZI=y(5);
    dy=[
        -beta*I*S; 
        beta*I*S - gammaI*I; 
        pHospitalization*gammaI*I - gammaH*H; 
        (1-pHospitalization)*gammaI*I + gammaH*H; 
        gammaI*I
    ];
end