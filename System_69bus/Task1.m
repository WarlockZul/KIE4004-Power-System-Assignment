%% Task 1: System Modeling, Analysis, and Excel Export
clc; clear; close all;

% 1. Setup and Load Data
% addpath(genpath('C:\Program Files\Matlab\matpower8.1'));
% savepath;
define_constants;

mpc = case69;

nb = size(mpc.bus, 1); % Get system size

%% 2. Run Power Flow Methods
% Newton-Raphson (NR)
mpopt_nr = mpoption('pf.alg', 'NR', 'verbose', 0, 'out.all', 0);
tic; 
res_nr = runpf(mpc, mpopt_nr);
time_nr = toc;

% Fast Decoupled (FDLF)
mpopt_fd = mpoption('pf.alg', 'FDXB', 'verbose', 0, 'out.all', 0);
tic; 
res_fd = runpf(mpc, mpopt_fd);
time_fd = toc;

%% 3. Numerical Accuracy Calculations
% Calculate difference between NR and FD voltage magnitudes
v_diff = res_nr.bus(:, VM) - res_fd.bus(:, VM);

max_v_err = max(abs(v_diff));   % Maximum absolute error
rms_v_err = sqrt(mean(v_diff.^2)); % Root Mean Square error


% Loss Function: P_loss = -(P_from + P_to), Q_loss = -(Q_from + Q_to)
compute_losses = @(branch) deal( ...
    -sum(branch(:, PF) + branch(:, PT)), ...   % P loss MW
    -sum(branch(:, QF) + branch(:, QT)) ...    % Q loss Mvar
);

% Calculate NR Losses
[p_loss_nr, q_loss_nr] = compute_losses(res_nr.branch);
[p_loss_fd, q_loss_fd] = compute_losses(res_fd.branch);

% Calculate Power Balance: Error = |Generation - (Load + |Loss|)|
dp_nr = abs(sum(res_nr.gen(:, PG)) - (sum(res_nr.bus(:, PD)) + abs(p_loss_nr)));
dp_fd = abs(sum(res_fd.gen(:, PG)) - (sum(res_fd.bus(:, PD)) + abs(p_loss_fd)));
dq_nr = abs(sum(res_nr.gen(:, QG)) - (sum(res_nr.bus(:, QD)) + abs(q_loss_nr)));
dq_fd = abs(sum(res_fd.gen(:, QG)) - (sum(res_fd.bus(:, QD)) + abs(q_loss_fd)));

%% 4. Excel to Export
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

filename = sprintf('%s_task1.xlsx', sys_name); % Generates IEEE118_task1.xlsx etc.

if exist(filename, 'file'), delete(filename); end

% Sheet 1: Voltage Bus Data
T_voltage = table((1:nb)', ...
    res_nr.bus(:, VM), res_nr.bus(:, VA), ...
    res_fd.bus(:, VM), res_fd.bus(:, VA), ...
    abs(v_diff),...
    'VariableNames', {'Bus_No', 'NR_Vmag_pu', 'NR_Angle_deg', 'FD_Vmag_pu', 'FD_Angle_deg', 'Abs_Error'});
writetable(T_voltage, filename, 'Sheet', 'Voltage_Bus');

% Sheet 2: Summary
T_summary = table(["Newton-Raphson"; "Fast-Decoupled"], ...
    [res_nr.iterations; res_fd.iterations], ...
    [time_nr; time_fd], ...
    [p_loss_nr; p_loss_fd], [q_loss_nr; q_loss_fd], ...
    [dp_nr; dp_fd], ...
    [dq_nr; dq_fd], ...
    [0; max_v_err], [0; rms_v_err], ...
    'VariableNames', {'Method', 'Iterations', 'Runtime_s', 'P_Loss_MW', 'Q_Loss_Mvar', 'Delta_P', 'Delta_Q', 'Max_V_Error', 'RMS_V_Error'});
writetable(T_summary, filename, 'Sheet', 'Summary');

fprintf('Successfully exported results for %s to: %s\n', sys_name, filename);

%% 5. Visualization and Auto-Save Voltage Profile
fig = figure('Color', 'w', 'Name', [sys_name, ' Voltage Profile']);
ax = axes(fig); % Define the axes handle

% Plotting the Voltage Profile
plot(ax, res_nr.bus(:, BUS_I), res_nr.bus(:, VM), 'b-s', 'MarkerSize', 4, 'LineWidth', 1);
grid on; hold on;

% Remove the toolbar from the saved image
ax.Toolbar.Visible = 'off'; 

% Add statutory limit lines
yline(ax, 0.95, 'r--', 'Lower Limit (0.95 p.u.)', 'LabelVerticalAlignment', 'bottom');
yline(ax, 1.05, 'r--', 'Upper Limit (1.05 p.u.)');

% Formatting the plot
title(ax, ['Task 1: Base Case Voltage Profile (', sys_name, ')']);
xlabel(ax, 'Bus Number');
ylabel(ax, 'Voltage Magnitude (p.u.)');
ylim(ax, [min(res_nr.bus(:, VM))-0.05, 1.1]); 

% Save the plot
plot_filename = sprintf('%s_Task1_VP.png', sys_name);
exportgraphics(fig, plot_filename, 'Resolution', 300); % Higher quality than saveas

fprintf('Voltage profile plot saved as: %s\n', plot_filename);