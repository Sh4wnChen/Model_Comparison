%% Batch run Yates09 & ShoreFor models (3-year segments) -- Pier-calibrated only
% 10 three-year periods, calibrated using pier data only
% Pier-calibrated model output is compared against both pier obs and FT mean
% to show that model output aligns better with FT mean (cusp noise is not fitted)
%
% Outputs:
%   - Model_<mdl>_pier_<period>.csv : daily model output
%   - model_performance_3yr.csv     : R^2/RMSE summary table
%   - Time series comparison plots (5x2) : one figure per model
%   - RMSE_pier vs RMSE_FT scatter plot  : key evidence figure

clear all;

%% ========== Global configuration ==========
periods_start = datetime(1990:3:2017, 1, 1)';
periods_end   = datetime(1993:3:2020, 1, 1)';
n_periods = length(periods_start);

shore_pier_file = 'Hasaki_Shore_ref1.40.csv';
shore_ft_file   = 'FT_mean_shoreline.csv';
wave_file       = 'Hasaki_Wave_JRA55.csv';

first_day = datetime(1990,1,1);

%% Read wave data
wave_data = readtable(wave_file, 'DatetimeType', 'datetime');
wave_tt_all = table2timetable(wave_data, 'RowTimes', 'time');
wave_tt_daily_all = retime(wave_tt_all, 'daily', 'mean');

w = 0.15; % fall velocity
wave_tt_daily_all.Omega = wave_tt_daily_all.Hs ./ (w * wave_tt_daily_all.Tp);
wave_tt_daily_all.F = ones(height(wave_tt_daily_all), 1);
wave_dt_seconds = seconds(median(diff(wave_tt_daily_all.time)));

%% Read shoreline data (pier for calibration; FT mean for validation)
shore_pier = table2timetable(readtable(shore_pier_file, 'DatetimeType', 'datetime'), ...
    'RowTimes', 'time');
shore_ft   = table2timetable(readtable(shore_ft_file,   'DatetimeType', 'datetime'), ...
    'RowTimes', 'time');

%% ========== R^2 summary table ==========
% vs_Pier : model (pier-calibrated) vs pier observations
% vs_FT   : model (pier-calibrated) vs FT mean observations
period_labels = cell(n_periods, 1);
for p = 1:n_periods
    period_labels{p} = [datestr(periods_start(p),'yyyy') '-' datestr(periods_end(p),'yyyy')];
end

R2_table = table(period_labels, ...
    NaN(n_periods,1), NaN(n_periods,1), NaN(n_periods,1), NaN(n_periods,1), ...
    NaN(n_periods,1), NaN(n_periods,1), NaN(n_periods,1), NaN(n_periods,1), ...
    NaN(n_periods,1), NaN(n_periods,1), ...
    'VariableNames', {'Period', ...
    'Yates09_R2_vs_Pier',  'Yates09_R2_vs_FT', ...
    'Yates09_RMSE_vs_Pier','Yates09_RMSE_vs_FT', ...
    'ShoreFor_R2_vs_Pier', 'ShoreFor_R2_vs_FT', ...
    'ShoreFor_RMSE_vs_Pier','ShoreFor_RMSE_vs_FT', ...
    'N_Pier_obs', 'N_FT_obs'});

%% ==================== Yates09 (pier-calibrated) ====================
fprintf('\n========== Yates09 (Pier-calibrated) ==========\n');

for p = 1:n_periods
    t_start = periods_start(p);
    t_end   = periods_end(p);
    fprintf('\n--- Yates09 | %s ---\n', period_labels{p});

    cal_tr = timerange(t_start, t_end, 'closed');
    pier_cal = shore_pier(cal_tr, :);
    ft_cal   = shore_ft(cal_tr, :);
    wave_cal = wave_tt_all(cal_tr, :);
    wave_cal.E = wave_cal.Hs.^2;
    wave_daily = retime(wave_cal, 'daily', 'mean');
    dt_days = 1;

    if height(pier_cal) < 3
        fprintf('  Insufficient pier observations, skipping\n');
        continue;
    end

    % --- Optimization (using pier data) ---
    initialYini = mean(pier_cal.y);
    x0 = [-0.005, 1, -1, -2, initialYini];
    lb = [-0.1, -30, -50, -10, min(pier_cal.y) - 10];
    ub = [0, 10, 0, 0, max(pier_cal.y) + 10];

    opts = optimoptions('fmincon', 'Display', 'off', ...
        'MaxIterations', 1000, 'MaxFunctionEvaluations', 10000, ...
        'OptimalityTolerance', 1e-6, 'Algorithm', 'sqp');

    objFunc = @(params) yates09_objective(params, wave_daily, dt_days, pier_cal);
    problem = createOptimProblem('fmincon', 'objective', objFunc, ...
        'x0', x0, 'lb', lb, 'ub', ub, 'options', opts);

    ms = MultiStart('Display', 'off', 'UseParallel', true, 'MaxTime', 600);
    [optimParams, ~] = run(ms, problem, 5);

    % --- Daily model output ---
    [Y_opt, ~] = yates09(wave_daily.E, dt_days, ...
        optimParams(1), optimParams(2), optimParams(3), optimParams(4), optimParams(5));

    % --- vs pier observations ---
    [r2_pier, rmse_pier, n_pier] = compute_metrics(wave_daily.time, Y_opt, pier_cal);
    R2_table.Yates09_R2_vs_Pier(p)   = r2_pier;
    R2_table.Yates09_RMSE_vs_Pier(p) = rmse_pier;
    R2_table.N_Pier_obs(p) = n_pier;

    % --- vs FT mean ---
    [r2_ft, rmse_ft, n_ft] = compute_metrics(wave_daily.time, Y_opt, ft_cal);
    R2_table.Yates09_R2_vs_FT(p)   = r2_ft;
    R2_table.Yates09_RMSE_vs_FT(p) = rmse_ft;
    R2_table.N_FT_obs(p) = n_ft;

    % --- Save model time series ---
    out_tt = timetable(wave_daily.time, Y_opt, 'VariableNames', {'y_model'});
    out_tt.Properties.DimensionNames{1} = 'time';
    fname = ['Model_Yates09_pier_' period_labels{p} '.csv'];
    writetable(timetable2table(out_tt), fname);

    fprintf('  vs Pier: R^2=%.3f RMSE=%.2f (n=%d) | vs FT: R^2=%.3f RMSE=%.2f (n=%d)\n', ...
        r2_pier, rmse_pier, n_pier, r2_ft, rmse_ft, n_ft);
end

%% ==================== ShoreFor (pier-calibrated) ====================
fprintf('\n========== ShoreFor (Pier-calibrated) ==========\n');

phi_test = [5:5:50 60:10:100 150:50:1000];
nn = 0.5;

for p = 1:n_periods
    t_start = periods_start(p);
    t_end   = periods_end(p);
    fprintf('\n--- ShoreFor | %s ---\n', period_labels{p});

    cal_tr = timerange(t_start, t_end, 'closed');
    pier_cal = shore_pier(cal_tr, :);
    ft_cal   = shore_ft(cal_tr, :);
    wave_cal = wave_tt_daily_all(cal_tr, :);
    cal_indices = ismember(wave_tt_daily_all.time, wave_cal.time);

    [common_times, wave_idx, shore_idx] = intersect(wave_cal.time, pier_cal.time);
    if length(common_times) < 3
        fprintf('  Insufficient overlapping data points, skipping\n');
        continue;
    end

    best_error = Inf;
    best_phi = NaN; best_c1 = NaN(3,1); best_c2 = NaN(3,1);

    for i = 1:length(phi_test)
        D = 2 * phi_test(i);
        omega_eq = WS85FilterConv(wave_tt_daily_all.Omega, D, phi_test(i), wave_dt_seconds);
        omega_eq = omega_eq(cal_indices);

        ie = omega_eq - wave_cal.Omega <= 0;
        ia = omega_eq - wave_cal.Omega > 0;

        DeltaX = wave_cal.F .^ nn .* (omega_eq - wave_cal.Omega);
        DeltaX_acc = DeltaX; DeltaX_acc(ie) = 0;
        DeltaX_ero = DeltaX; DeltaX_ero(ia) = 0;

        cum_dx_acc = cumsum(DeltaX_acc);
        cum_dx_ero = cumsum(DeltaX_ero);

        cum_acc_c = cum_dx_acc(wave_idx);
        cum_ero_c = cum_dx_ero(wave_idx);
        t_day = days(common_times - first_day);
        y_obs = pier_cal.y(shore_idx);

        % Two-step correction
        y_obs_dt = remove_linear_trend(pier_cal.y, pier_cal.time);
        y_obs_dt = y_obs_dt(shore_idx);

        A1 = [cum_acc_c, cum_ero_c, ones(length(t_day),1)];
        c1 = A1 \ y_obs_dt;
        y_basic = c1(1)*cum_acc_c + c1(2)*cum_ero_c + c1(3);

        A2 = [y_basic, t_day, ones(length(t_day),1)];
        c2 = A2 \ y_obs;

        y_model_obs = c2(1)*y_basic + c2(2)*t_day + c2(3);
        err = sqrt(mean((y_model_obs - y_obs).^2));

        if err < best_error
            best_error = err;
            best_phi = phi_test(i);
            best_c1 = c1; best_c2 = c2;
        end
    end

    % --- Compute daily model output with best phi ---
    omega_eq = WS85FilterConv(wave_tt_daily_all.Omega, 2*best_phi, best_phi, wave_dt_seconds);
    omega_eq = omega_eq(cal_indices);
    ie = omega_eq - wave_cal.Omega <= 0;
    ia = omega_eq - wave_cal.Omega > 0;
    DeltaX = wave_cal.F .^ nn .* (omega_eq - wave_cal.Omega);
    DeltaX_acc = DeltaX; DeltaX_acc(ie) = 0;
    DeltaX_ero = DeltaX; DeltaX_ero(ia) = 0;
    cum_dx_acc = cumsum(DeltaX_acc);
    cum_dx_ero = cumsum(DeltaX_ero);

    t_day_daily = days(wave_cal.time - first_day);
    y_basic_daily = best_c1(1)*cum_dx_acc + best_c1(2)*cum_dx_ero + best_c1(3);
    y_daily_model = best_c2(1)*y_basic_daily + best_c2(2)*t_day_daily + best_c2(3);

    % --- vs pier observations ---
    [r2_pier, rmse_pier, n_pier] = compute_metrics(wave_cal.time, y_daily_model, pier_cal);
    R2_table.ShoreFor_R2_vs_Pier(p)   = r2_pier;
    R2_table.ShoreFor_RMSE_vs_Pier(p) = rmse_pier;

    % --- vs FT mean ---
    [r2_ft, rmse_ft, n_ft] = compute_metrics(wave_cal.time, y_daily_model, ft_cal);
    R2_table.ShoreFor_R2_vs_FT(p)   = r2_ft;
    R2_table.ShoreFor_RMSE_vs_FT(p) = rmse_ft;

    % --- Save ---
    out_tt = timetable(wave_cal.time, y_daily_model, 'VariableNames', {'y_model'});
    out_tt.Properties.DimensionNames{1} = 'time';
    fname = ['Model_ShoreFor_pier_' period_labels{p} '.csv'];
    writetable(timetable2table(out_tt), fname);

    fprintf('  vs Pier: R^2=%.3f RMSE=%.2f (n=%d) | vs FT: R^2=%.3f RMSE=%.2f (n=%d) | phi=%d\n', ...
        r2_pier, rmse_pier, n_pier, r2_ft, rmse_ft, n_ft, best_phi);
end

%% ==================== Export summary table ====================
writetable(R2_table, 'model_performance_3yr.csv');
fprintf('\n========== Summary table saved: model_performance_3yr.csv ==========\n');
disp(R2_table);

%% ==================== Time series comparison plots (5x2) ====================
c_obs_pier = [0.75 0.75 0.75];
c_edge     = [0.5  0.5  0.5];
c_obs_ft   = [50,  50,  50] / 255;
c_model    = [72, 128, 184] / 255;

mdl_names = {'Yates09', 'ShoreFor'};
for m = 1:length(mdl_names)
    mdl = mdl_names{m};

    fig = figure('Position', [50, 50, 900, 1000], 'Color', 'w');
    tl = tiledlayout(5, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, [mdl ' Model'], ...
        'FontSize', 16, 'FontWeight', 'bold');

    for p = 1:n_periods
        ax = nexttile;
        hold on;
        tr = timerange(periods_start(p), periods_end(p), 'closed');

        ax.Color = [0.94 0.94 0.94];
        ax.Box = 'off';
        ax.FontSize = 11;
        ax.TickDir = 'out';
        ax.XAxis.LineWidth = 0.8;
        ax.YAxis.LineWidth = 0.8;
        ax.GridColor = [1 1 1];
        ax.GridAlpha = 1;
        ax.Layer = 'bottom';

        % Pier obs
        obs_p = shore_pier(tr, :);
        h_obs_pier = [];
        if height(obs_p) > 0
            h_obs_pier = scatter(obs_p.time, obs_p.y, 15, ...
                'MarkerFaceColor', c_obs_pier, 'MarkerEdgeColor', c_edge, ...
                'Marker', 'o', 'DisplayName', 'Obs (Pier)');
        end

        % FT obs
        obs_f = shore_ft(tr, :);
        h_obs_ft = [];
        if height(obs_f) > 0
            h_obs_ft = scatter(obs_f.time, obs_f.y, 35, ...
                'MarkerFaceColor', c_obs_ft, 'MarkerEdgeColor', c_edge, ...
                'Marker', 's', 'DisplayName', 'Obs (MT)');
        end

        % Model
        fname = ['Model_' mdl '_pier_' period_labels{p} '.csv'];
        h_mod = [];
        if exist(fname, 'file')
            tmp = readtable(fname, 'DatetimeType', 'datetime');
            mod_tt = table2timetable(tmp, 'RowTimes', 'time');
            h_mod = plot(mod_tt.time, mod_tt.y_model, '-', ...
                'Color', c_model, 'LineWidth', 1.4, ...
                'DisplayName', [mdl ' (Pier-cal)']);
        end

        xlim([periods_start(p) periods_end(p)]);
        ylim([-40 40]);
        years_vec = year(periods_start(p)):year(periods_end(p));
        xticks(datetime(years_vec, 1, 1));
        xtickformat('yyyy');

        % Period label
        text(0.5, 0.95, period_labels{p}, ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', 10, 'FontAngle', 'italic', 'Color', [0.45 0.45 0.45]);

        % Metrics: vs Pier (top line) and vs FT (bottom line)
        if strcmp(mdl, 'Yates09')
            r2p = R2_table.Yates09_R2_vs_Pier(p);
            rmp = R2_table.Yates09_RMSE_vs_Pier(p);
            r2f = R2_table.Yates09_R2_vs_FT(p);
            rmf = R2_table.Yates09_RMSE_vs_FT(p);
        else
            r2p = R2_table.ShoreFor_R2_vs_Pier(p);
            rmp = R2_table.ShoreFor_RMSE_vs_Pier(p);
            r2f = R2_table.ShoreFor_R2_vs_FT(p);
            rmf = R2_table.ShoreFor_RMSE_vs_FT(p);
        end

        if ~isnan(r2p)
            text(0.05, 0.13, sprintf('vs Pier: R^2=%.2f RMSE=%.1f', r2p, rmp), ...
                'Units', 'normalized', 'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'bottom', 'FontSize', 9, ...
                'FontAngle', 'italic', 'Color', [0.3 0.3 0.3]);
        end
        if ~isnan(r2f)
            text(0.05, 0.05, sprintf('vs MT:   R^2=%.2f RMSE=%.1f', r2f, rmf), ...
                'Units', 'normalized', 'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'bottom', 'FontSize', 9, ...
                'FontAngle', 'italic', 'Color', [0.3 0.3 0.3]);
        end

        grid on;
        if mod(p,2) == 1, ylabel('Position (m)'); end
        if p >= 9, xlabel('Time'); end

        if p == 2 && ~isempty(h_obs_pier) && ~isempty(h_obs_ft) && ~isempty(h_mod)
            legend([h_obs_pier, h_obs_ft, h_mod], ...
                'Obs (Pier)', 'Obs (MT)', [mdl ' (Pier-cal)'], ...
                'NumColumns', 2, 'FontSize', 7, 'Box', 'off', ...
                'Location', 'southeast');
        end

        hold off;
    end
end



%% ==================== Local functions ====================

function [r2, rmse, n] = compute_metrics(t_model, y_model, obs_tt)
    % Interpolate model time series to observation times, compute R^2 and RMSE
    r2 = NaN; rmse = NaN; n = 0;
    if height(obs_tt) < 3, return; end
    y_at_obs = interp1(datenum(t_model), y_model, datenum(obs_tt.time), 'linear', NaN);
    valid = ~isnan(y_at_obs) & ~isnan(obs_tt.y);
    n = sum(valid);
    if n < 3, return; end
    y_m = y_at_obs(valid); y_o = obs_tt.y(valid);
    ss_res = sum((y_o - y_m).^2);
    ss_tot = sum((y_o - mean(y_o)).^2);
    r2 = 1 - ss_res / ss_tot;
    rmse = sqrt(mean((y_m - y_o).^2));
end

function [Y, Seq] = yates09(E, dt, a, b, cacr, cero, Yini)
    Seq = (E - b) / a;
    Y = zeros(size(E));
    Y(1) = Yini;
    for i = 1:(length(E)-1)
        if Y(i) < Seq(i+1)
            Y(i+1) = ((Y(i)-Seq(i+1)) * exp(-1*a*cacr*(E(i+1)^0.5)*dt)) + Seq(i+1);
        else
            Y(i+1) = ((Y(i)-Seq(i+1)) * exp(-1*a*cero*(E(i+1)^0.5)*dt)) + Seq(i+1);
        end
    end
end

function err = yates09_objective(params, E_tt, dt_days, Yobs_tt)
    a = params(1); b = params(2);
    cacr = params(3); cero = params(4); Yini = params(5);
    [Y_model, ~] = yates09(E_tt.E, dt_days, a, b, cacr, cero, Yini);
    Y_model_at_obs = interp1(datenum(E_tt.time), Y_model, ...
        datenum(Yobs_tt.time), 'linear', 'extrap');
    valid = ~isnan(Y_model_at_obs) & ~isnan(Yobs_tt.y);
    if sum(valid) < 2, err = 1e6; return; end
    err = sqrt(mean((Y_model_at_obs(valid) - Yobs_tt.y(valid)).^2));
end

function omegaFiltered = WS85FilterConv(omega, D, phi, dt)
    % Equilibrium omega via Wright & Short (1985) exponential filter
    % Uses convolution for speed (Davidson, 2012)
    dt = dt / (3600 * 24);
    D = round(D / dt);
    phi = round(phi / dt);
    meanOmega = mean(omega);
    omega = omega - meanOmega;
    ii = 0:(D-1);
    padding = zeros(D-1, 1);
    filterCoeff = 10.^(-abs(ii) / phi);
    filterCoeff = [padding(:); filterCoeff(:)];
    window = filterCoeff / sum(filterCoeff);
    omegaFiltered = conv(omega, window, 'same');
    omegaFiltered = omegaFiltered + meanOmega;
end

function detrended = remove_linear_trend(s_obs, t_obs)
    t_num = datenum(t_obs);
    p = polyfit(t_num, s_obs, 1);
    detrended = s_obs - polyval(p, t_num);
end