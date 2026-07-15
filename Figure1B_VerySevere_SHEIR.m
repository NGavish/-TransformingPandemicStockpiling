function Figure1_VerySevere_Revised_SEIHR
% --- Setup and Data Initialization ---
close all;
figure('Color', 'w', 'Position', [100, 100, 1400, 500]);
t = tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel A: Static Model ---
nexttile
lg = plotPanelA;
text(0.05, 0.95, 'A', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');
box on;

% --- Panel B: Dynamic Model (SEIHR) ---
nexttile
% Base parameters
gammaI=1/2; 
gammaH=1/4; 
pHospitalization=0.12; 
duration=90; 
clinicalAttackRate=0.25; 
N=320e6;

% NEW: Incubation period parameter (sigma)
% Assumes a 3-day incubation period (1/3). Adjust as needed.
sigma = 1/3; 

plotPanelB(pHospitalization, duration, clinicalAttackRate, gammaI, gammaH, N, sigma);
text(0.05, 0.95, 'B', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');
box on;

% --- Panel C: Historical Reconstruction ---
nexttile
plotPanelC();
text(0.05, 0.95, 'C', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');
box on;

% --- Shared Legend ---
lg.Layout.Tile = 'South';
lg.Orientation = 'horizontal';

% Save (Optional)
printGraph('./graphs/Figure1_VerySevere_Revised_SHEIR');

end

function lg = plotPanelA
    % Data for Very Severe (Index 3 from original code logic)
    peak_no_av = [0.8+1.2, 3.8+5.8, 7.7+11.5]/2/42*4*1e6; 
    av_min = 0.4*peak_no_av;
    av_max = 0.7*peak_no_av;
    av_mid = (av_min + av_max) / 2;
    err_low = av_mid - av_min;
    err_high = av_max - av_mid;
    
    val_unmit = peak_no_av(3);
    val_mit = av_mid(3);
    err_l = err_low(3);
    err_h = err_high(3);

    hold on;
    % Background Zones
    x_limits = [0 2 2 0 0]; 
    p6 = fill(x_limits, [5e4 5e4 1e5 1e5 5e4], [1 1 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    p7 = fill(x_limits, [1e5 1e5 1.5e5 1.5e5 1e5], [1 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    p_line = yline(1.5e5, 'r--', 'LineWidth', 2);
    
    % Bars (Centered closer together)
    b1 = bar(0.85, val_unmit, 0.25);
    b2 = bar(1.15, val_mit, 0.25);
    
    % Style
    b1.FaceColor = [0.9, 0.9, 0.9]; b1.EdgeColor = [55 112 184]/255; b1.LineWidth = 1.5; b1.FaceAlpha = 0.7;
    b2.FaceColor = [0.4, 0.4, 0.4]; b2.EdgeColor = 'k'; b2.LineWidth = 1; b2.FaceAlpha = 0.7;

    % Errorbar
    errorbar(1.15, val_mit, err_l, err_h, 'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

    % Formatting
    ylabel('Daily Bed Demand', 'FontSize', 12);
    title('Static Model (Lower Bound)', 'FontSize', 12);
    
    ylim([0 2.5e6]); 
    set(gca, 'YTick', 0:5e5:2e6);
    ytickformat('%.1fM');
    labels = arrayfun(@(v) sprintf('%.1fM', v/1e6), get(gca, 'YTick'), 'UniformOutput', false);
    set(gca, 'YTickLabel', labels);
    
    xlim([0.5 1.5]);
    set(gca, 'XTick', []); 
    grid on; ax = gca; ax.YGrid = 'on'; ax.XGrid = 'off';
    
    % Legend handle
    lg = legend([b1, b2, p6, p7, p_line], ...
        {'Unmitigated', 'With Antivirals', 'Severe Strain Zone', 'Crisis Zone', 'Collapse Threshold'}, ...
        'FontSize', 10);
    hold off;
end

function plotPanelB(pHospitalization, duration, clinicalAttackRate, gammaI, gammaH, N, sigma)
    betaParm = 2*gammaI;
    % Pass sigma into the findR0 function
    beta = fzero(@(betaParm)findR0(betaParm,sigma,gammaI,gammaH,pHospitalization), betaParm);
    
    % Initial state: S=1-I0, E=0, I=I0, H=0, R=0, ZI=0
    I0 = 1e-3; 
    y0 = [1-I0; 0; I0; 0; 0; 0];
    
    [t, y] = ode45(@(t,y)SEIHR(t,y,beta,sigma,gammaI,gammaH,pHospitalization), [0 duration], y0);
    
    Z = y(end,6); % ZI is now at index 6
    [t_mit, y_mit] = ode45(@(t,y)SEIHR(t,y,beta,sigma,gammaI,gammaH,pHospitalization*clinicalAttackRate/Z), [0 duration], y0);
    t = t_mit/7; 
    H = y_mit(:,4); % H is now at index 4
    
    hold on;
    % Backgrounds
    x_lims = [0 8 8 0 0];
    fill(x_lims, [5e4 5e4 1e5 1e5 5e4], [1 1 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill(x_lims, [1e5 1e5 1.5e5 1.5e5 1e5], [1 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    yline(1.5e5, 'r--', 'LineWidth', 2);
    text(0.2, 1.6e5, 'Collapse Threshold (~150k)', 'Color', 'r', 'FontSize', 9, 'VerticalAlignment', 'bottom');

    % Curves
    h_lower = N * H * 0.4; % Best case (60% efficacy)
    h_upper = N * H * 0.7; % Worst case (30% efficacy)
    fill([t; flipud(t)], [h_upper; flipud(h_lower)], [0.4 0.4 0.4], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    plot(t, N*H, 'Color', [55 112 184]/255, 'LineWidth', 2);
    
    % Annotate Gap ("Unmet Demand")
    [max_val, idx] = max(h_lower);
    t_peak = t(idx);
    plot([t_peak t_peak], [1.5e5 max_val], 'k-', 'LineWidth', 1.5);
    plot(t_peak, 1.5e5, 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    plot(t_peak, max_val, 'k^', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    
    text(t_peak+0.1, (1.5e5 + max_val)/2, ' Unmet Demand', 'FontSize', 10, 'FontWeight', 'bold');
    text(t_peak+0.1, (1.5e5 + max_val)/2 - 1.5e5, ' (~5x Capacity)', 'FontSize', 9);

    % Formatting
    title('Dynamic Model (SEIHR, CDC Very Severe)', 'FontSize', 12);
    ylabel('Daily Bed Demand');
    xlabel('Weeks into wave');
    ylim([0 2.5e6]);
    set(gca, 'YTick', 0:5e5:2e6);
    ytickformat('%.1fM');
    labels = arrayfun(@(v) sprintf('%.1fM', v/1e6), get(gca, 'YTick'), 'UniformOutput', false);
    set(gca, 'YTickLabel', labels);
    xlim([0 8]);
    grid on;
end

function plotPanelC()
    % Setup Data
    pop_old = 310e3; pop_new = 320e3;
    scale_factor = pop_new / pop_old/7/0.24;
    west_coast=0.24;
    data_1918_west_orig = [3 20e3; 4 40e3; 5 120e3; 6 220e3; 7 310e3; 8 240e3; 9 140e3; 10 70e3; 11 40e3; 12 25e3; 13 20e3];
    data_1918_west_orig(:,2)=data_1918_west_orig(:,2)*west_coast;
    hosp_1918_west_orig = mortality_to_hospitalization(data_1918_west_orig);
    hosp_1918_west = hosp_1918_west_orig; 
    hosp_1918_west(:,2) = hosp_1918_west(:,2) * scale_factor;

    hold on;
    % Background
    x_lims = [0 10 10 0 0];
    fill(x_lims, [5e4 5e4 1e5 1e5 5e4], [1 1 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill(x_lims, [1e5 1e5 1.5e5 1.5e5 1e5], [1 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    yline(1.5e5, 'r--', 'LineWidth', 2);

    % Curves
    x = hosp_1918_west(:,1)/7-2;
    y = hosp_1918_west(:,2);
    h_lower = y * 0.4;
    h_upper = y * 0.7;
    
    fill([x; flipud(x)], [h_upper; flipud(h_lower)], [0.4 0.4 0.4], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    plot(x, y, 'Color', [55 112 184]/255, 'LineWidth', 2);
    
    title('1918 Historical Reconstruction', 'FontSize', 12);
    ylabel('Daily Bed Demand');
    xlabel('Weeks into wave');
    ylim([0 2.5e6]);
    set(gca, 'YTick', 0:5e5:2e6);
    ytickformat('%.1fM');
    labels = arrayfun(@(v) sprintf('%.1fM', v/1e6), get(gca, 'YTick'), 'UniformOutput', false);
    set(gca, 'YTickLabel', labels);
    xlim([0 10]);
    grid on;
end

% --- Helpers ---
function hosp_data = mortality_to_hospitalization(weekly_mortality)
    if isempty(weekly_mortality); hosp_data = []; return; end
    first_data_week = weekly_mortality(1,1);
    prepended_points = [];
    if first_data_week > 1; prepended_points = [first_data_week - 1, 0]; end
    if first_data_week > 0
       if isempty(prepended_points) || prepended_points(1,1) > 0; prepended_points = [0 0; prepended_points]; end
    end
    weekly_mortality = [prepended_points; weekly_mortality];
    [~, unique_indices] = unique(weekly_mortality(:,1));
    weekly_mortality = weekly_mortality(unique_indices,:);

    weeks = weekly_mortality(:,1); deaths = weekly_mortality(:,2);
    days_sparse = weeks * 7;
    if ~isempty(days_sparse) && any(days_sparse == 0); days_sparse(days_sparse == 0) = 1; end
    if length(days_sparse) == 1; days_sparse = [days_sparse; days_sparse + 1]; deaths = [deaths; deaths]; end
    days_continuous = min(days_sparse):max(days_sparse);
    daily_deaths = interp1(days_sparse, deaths, days_continuous, 'pchip');
    daily_deaths(daily_deaths < 0) = 0;
    hospitalization_series = 10 * conv(daily_deaths, ones(1, 4), 'full');
    daily_hosps = hospitalization_series(1:length(days_continuous));
    hosp_data = [days_continuous', daily_hosps'];
end

% NEW: Added sigma to parameters
function err=findR0(betaParm,sigma,gammaI,gammaH,pHospitalization)
    I0=1e-3; 
    y0=[1-I0;0;I0;0;0;0]; % E starts at 0
    [t,y]=ode45(@(t,y)SEIHR(t,y,betaParm,sigma,gammaI,gammaH,pHospitalization),[0 3000],y0,odeset('RelTol',1e-10,'Events', @(t,y) peakEventsFcn(t,y,betaParm,sigma,gammaI,gammaH,pHospitalization)));
    err=t(end)-21; % Anchor peak to day 21
end

% NEW: Added sigma to parameters and adjusted index
function [position,isterminal,direction] = peakEventsFcn(t,y,beta,sigma,gammaI,gammaH,pHospitalization)
    dydt = SEIHR(t,y,beta,sigma,gammaI,gammaH,pHospitalization);
    % We want the peak of Infectious compartment. I is now at index 3.
    position = dydt(3); 
    isterminal = 1; 
    direction = -1;
end

% NEW: The actual SEIHR differential equation
function dy=SEIHR(t,y,beta,sigma,gammaI,gammaH,pHospitalization)
    S=y(1); E=y(2); I=y(3); H=y(4); R=y(5); ZI=y(6);
    
    dy = [
        -beta*I*S;                                      % dS
         beta*I*S - sigma*E;                            % dE (Exposed)
         sigma*E - gammaI*I;                            % dI (Infectious)
         pHospitalization*gammaI*I - gammaH*H;          % dH (Hospitalized)
         (1-pHospitalization)*gammaI*I + gammaH*H;      % dR (Recovered)
         gammaI*I                                       % dZI (Cumulative transitions out of I)
    ];
end

% function printGraph(name)
%     print(name,'-dpng','-r300');
% end