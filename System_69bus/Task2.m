%% Task 2: Distributed RE Integration Analysis (L-index > 0.04)
clc; clear; close all;
% addpath(genpath('C:\Program Files\Matlab\matpower8.1'));
% savepath;
define_constants;

mpc = case69;

nb = size(mpc.bus, 1); % Get system size

%% 1. Base Case Analysis
mpopt = mpoption('pf.alg','NR', 'verbose',0, 'out.all',0);
res_base = runpf(mpc, mpopt);
V_base = res_base.bus(:, VM);
L_base_all = compute_L_index(mpc, res_base);

% System Losses
p_loss_base = -sum(res_base.branch(:, PF) + res_base.branch(:, PT));
q_loss_base = -sum(res_base.branch(:, QF) + res_base.branch(:, QT));

% Power Balance Mismatch
dp_base = abs(sum(res_base.gen(:, PG)) - (sum(res_base.bus(:, PD)) + abs(p_loss_base)));
dq_base = abs(sum(res_base.gen(:, QG)) - (sum(res_base.bus(:, QD)) + abs(q_loss_base)));

%% 2. Dynamic Targeting for Distributed RE Placement
% SCENARIO A: STABILITY-FOCUSED PLACEMENT

% Set L-target based on the system
if nb == 118, L_target = 0.05;
elseif nb == 69, L_target = 0.1; % DIY target for 69-bus
elseif nb == 33, L_target = 0.085; % DIY target for 33-bus
else, L_target = 0.05; % Default fallback
end
% Identify weak buses
weak_buses = find(L_base_all > L_target);
num_weak = length(weak_buses);

% Avoid division by zero if no buses meet the target
if num_weak > 0
    P_total = 10.0; % Total RE capacity in MW
    P_per_bus = P_total / num_weak;
    fprintf('Identified %d weak buses for (L > %.3f)\n', num_weak, L_target);
    fprintf('Buses: %s\n', num2str(weak_buses'));
    fprintf('Injecting %.2f MW at each weak bus.\n\n', P_per_bus);
else
    error('No buses found with L-index > %.3f. Try lowering L_target.\n', L_target);
end

% Insert Generator at all chosen weak busses based on L_target
mpc_re = mpc;
for i = 1:num_weak
    b_idx = weak_buses(i);
    mpc_re.bus(b_idx, BUS_TYPE) = PV; 
    new_gen = zeros(1, size(mpc.gen, 2));
    new_gen(GEN_BUS) = b_idx; 
    new_gen(PG) = P_per_bus;
    new_gen(QMAX) = 10; 
    new_gen(QMIN) = -10; 
    new_gen(VG) = 1.02; 
    new_gen(GEN_STATUS) = 1;
    mpc_re.gen = [mpc_re.gen; new_gen];
    mpc_re.gencost = [mpc_re.gencost; zeros(1, size(mpc.gencost, 2))];
end

% SCENARIO B: LOSS-FOCUSED PLACEMENT
[opt_bus, min_loss_val, sensitivity_data] = run_loss_sensitivity_scan(mpc, mpopt);

fprintf('Optimal Location found at Bus %d.\n', opt_bus);

mpc_loss = mpc;
mpc_loss.bus(opt_bus, PD) = mpc_loss.bus(opt_bus, PD) - 1.0; 

%% 3. Post-Integration Analysis
% SCENARIO A: STABILITY-FOCUSED PLACEMENT
res_re = runpf(mpc_re, mpopt);
V_re = res_re.bus(:, VM);
L_re_all = compute_L_index(mpc_re, res_re);

% RE Case Losses
p_loss_re = -sum(res_re.branch(:, PF) + res_re.branch(:, PT));
q_loss_re = -sum(res_re.branch(:, QF) + res_re.branch(:, QT));

% RE Case Mismatch
dp_re = abs(sum(res_re.gen(:, PG)) - (sum(res_re.bus(:, PD)) + abs(p_loss_re)));
dq_re = abs(sum(res_re.gen(:, QG)) - (sum(res_re.bus(:, QD)) + abs(q_loss_re)));

% SCENARIO B: LOSS-FOCUSED PLACEMENT
res_loss = runpf(mpc_loss, mpopt);
V_lossOpt = res_loss.bus(:, VM); 
L_lossOpt_all = compute_L_index(mpc_loss, res_loss);

% Loss Case Mismatch
p_loss_lossOpt = -sum(res_loss.branch(:, PF) + res_loss.branch(:, PT));
q_loss_lossOpt = -sum(res_loss.branch(:, QF) + res_loss.branch(:, QT)); 

% Loss Case Mismatch
dp_lossOpt = abs(sum(res_loss.gen(:, PG)) - (sum(res_loss.bus(:, PD)) + abs(p_loss_lossOpt)));
dq_lossOpt = abs(sum(res_loss.gen(:, QG)) - (sum(res_loss.bus(:, QD)) + abs(q_loss_lossOpt)));

%% 4. Export for Report
% Automatically set the filename based on the number of buses
if nb == 118
    sys_name = 'IEEE118';
elseif nb == 69
    sys_name = 'IEEE69';
elseif nb == 33
    sys_name = 'IEEE33';
else
    sys_name = ['Custom', num2str(nb)];
end

filename = sprintf('%s_task2.xlsx', sys_name); % Generates IEEE118_task1.xlsx etc.

if exist(filename, 'file'), delete(filename); end


% Summary Table
T_sum = table(["Base Case"; "Scenario A (Stability)"; "Scenario B (Loss Opt)"], ...
    [p_loss_base; p_loss_re; p_loss_lossOpt], ...
    [q_loss_base; q_loss_re; q_loss_lossOpt], ...
    [max(L_base_all); max(L_re_all); max(L_lossOpt_all)], ...
    [dp_base; dp_re; dp_lossOpt], ...
    [dq_base; dq_re; dq_lossOpt], ...
    'VariableNames', {'Case','P_Loss_MW','Q_Loss_Mvar','Max_L_Index','Delta_P','Delta_Q'});

writetable(T_sum, filename, 'Sheet', 'Stability_Summary');

T_profile = table((1:size(mpc.bus,1))', V_base, V_re, V_lossOpt, L_base_all, L_re_all, L_lossOpt_all, ...
    'VariableNames', {'Bus','V_Base','V_ScenA','V_ScenB','L_Base','L_ScenA','L_ScenB'});
writetable(T_profile, filename, 'Sheet', 'Detailed_Data');

fprintf('Task 2 Complete. Results exported to: %s\n', filename);

%% 5. Visualization
% --- PLOT 1: Voltage Profile Comparison ---
fig1 = figure('Color', 'w', 'Name', 'Voltage Comparison');
ax1 = axes(fig1);
plot(ax1, V_base, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Base Case'); hold on;
plot(ax1, V_re, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Distributed RE');
yline(ax1, 0.95, 'k:', 'Lower Limit', 'LineWidth', 1);
grid on; ylabel(ax1, 'Voltage Magnitude (p.u.)'); xlabel(ax1, 'Bus Number');
title(ax1, ['Voltage Profile Comparison: Base vs RE (', sys_name, ')']);
legend(ax1, 'Location', 'southoutside', 'Orientation', 'horizontal');
ax1.Toolbar.Visible = 'off'; % Clean export
exportgraphics(fig1, sprintf('%s_Task2_VP.png', sys_name), 'Resolution', 300);

% --- PLOT 2: L-Index Stability Profile ---
fig2 = figure('Color', 'w', 'Name', 'Stability Profile');
ax2 = axes(fig2);
plot(ax2, L_base_all, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Base Case'); hold on;
plot(ax2, L_re_all, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Distributed RE');
grid on; ylabel(ax2, 'L-index (Stability Metric)'); xlabel(ax2, 'Bus Number');
title(ax2, ['Voltage Stability Improvement (', sys_name, ')']);
legend(ax2, 'Location', 'southoutside', 'Orientation', 'horizontal');
ax2.Toolbar.Visible = 'off';
exportgraphics(fig2, sprintf('%s_Task2_Lidx.png', sys_name), 'Resolution', 300);

% --- PLOT 3: Voltage Rise for Lidx (Impact Analysis) ---
% This plot shows exactly how much the voltage improved at each bus
fig3 = figure('Color', 'w', 'Name', 'Voltage Rise');
ax3 = axes(fig3);
stem(ax3, (V_re - V_base), 'Color', [0 0.5 0], 'Marker', 'none', 'LineWidth', 1);
grid on; ylabel(ax3, 'Delta Voltage (p.u.)'); xlabel(ax3, 'Bus Number');
title(ax3, ['Nodal Voltage Rise due to RE Integration (', sys_name, ')']);
ax3.Toolbar.Visible = 'off';
exportgraphics(fig3, sprintf('%s_Task2_VR after Lidx.png', sys_name), 'Resolution', 300);

% --- PLOT 4: Loss Sensitivity Scan (Scenario B Justification) ---
fig4 = figure('Color', 'w', 'Name', 'Loss Sensitivity B');
ax4 = axes(fig4);
plot(ax4, sensitivity_data(:,1), sensitivity_data(:,2), 'b-o', 'MarkerSize', 4);
grid on; ylabel(ax4, 'Total System Loss (MW)'); xlabel(ax4, 'Bus Number');
title(ax4, ['Loss Sensitivity Analysis: Scenario B (', sys_name, ')']);
xline(ax4, opt_bus, 'r--', ['Optimal Bus ' num2str(opt_bus)], 'LabelVerticalAlignment', 'bottom');
ax4.Toolbar.Visible = 'off';
exportgraphics(fig4, sprintf('%s_Task2_LossCurve_ScenB.png', sys_name), 'Resolution', 300);

% --- PLOT 5: Voltage Profile Comparison (Scenario B) ---
fig5 = figure('Color', 'w', 'Name', 'Voltage Comparison B');
ax5 = axes(fig5);
plot(ax5, V_base, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Base Case'); hold on;
plot(ax5, V_lossOpt, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Scenario B (Loss Opt)');
yline(ax5, 0.95, 'k:', 'Lower Limit', 'LineWidth', 1);
grid on; ylabel(ax5, 'Voltage Magnitude (p.u.)'); xlabel(ax5, 'Bus Number');
title(ax5, ['Voltage Profile: Scenario B (', sys_name, ')']);
legend(ax5, 'Location', 'southoutside', 'Orientation', 'horizontal');
ax5.Toolbar.Visible = 'off'; 
exportgraphics(fig5, sprintf('%s_Task2_VP_ScenB.png', sys_name), 'Resolution', 300);

fprintf('Task 2 individual plots saved for %s.\n', sys_name);

%% --- L-index Function ---
function L = compute_L_index(mpc, res)
    define_constants;
    Ybus = makeYbus(mpc);
    nb = size(res.bus, 1);
    gen_buses = res.gen(:, GEN_BUS);
    load_buses = setdiff((1:nb)', gen_buses);
    YLL = Ybus(load_buses, load_buses);
    YLG = Ybus(load_buses, gen_buses);
    F = -YLL \ YLG;
    V = res.bus(:, VM) .* exp(1j * res.bus(:, VA) * pi/180);
    L = zeros(nb, 1);
    for k = 1:length(load_buses)
        idx = load_buses(k);
        Li = 1 - sum(F(k,:) .* (V(gen_buses).' ./ V(idx)));
        L(idx) = abs(Li);
    end
end

%% --- System Losses ---
function [optimal_bus, min_loss, results] = run_loss_sensitivity_scan(mpc, mpopt)
    define_constants;
    nb = size(mpc.bus, 1);
    results = zeros(nb, 2); % [Bus_ID, Loss_Value]
    P_inject = 1.0; % Fixed size 1 MW
    
    for i = 1:nb
        mpc_test = mpc;
        bus_id = mpc.bus(i, BUS_I); 
        mpc_test.bus(i, PD) = mpc_test.bus(i, PD) - P_inject;
        res_test = runpf(mpc_test, mpopt);
        if res_test.success == 1
            loss = sum(res_test.branch(:, PF) + res_test.branch(:, PT));
            results(i, :) = [bus_id, loss];
        else
            results(i, :) = [bus_id, 9999]; 
        end
    end
    [min_loss, idx] = min(results(:, 2));
    optimal_bus = results(idx, 1);
end