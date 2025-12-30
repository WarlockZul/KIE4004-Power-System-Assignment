%% Task 1: System Modeling, Analysis, and Excel Export
clc; clear; close all;

% 1. Setup and Load Data
addpath(genpath('C:\Program Files\Matlab\matpower8.1')); 
savepath;
mpc = loadcase('case118');
define_constants; % Required for naming columns like VM, VA, PD, etc.

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

% System Losses
% % Power Balance Error Equation: |Gen - (Load + Loss)|
% calc_dP = @(res, p_loss) abs(sum(res.gen(:, PG)) - (sum(res.bus(:, PD)) + p_loss));
% calc_dQ = @(res, q_loss) abs(sum(res.gen(:, QG)) - (sum(res.bus(:, QD)) + q_loss));
% 
% p_loss_nr = sum(real(get_losses(res_nr)));
% q_loss_nr = sum(imag(get_losses(res_nr)));
% p_loss_fd = sum(real(get_losses(res_fd)));
% q_loss_fd = sum(imag(get_losses(res_fd)));
% 
% dp_nr = calc_dP(res_nr, p_loss_nr);
% dq_nr = calc_dQ(res_nr, q_loss_nr);
% dp_fd = calc_dP(res_fd, p_loss_fd);
% dq_fd = calc_dQ(res_fd, q_loss_fd);

% Loss Function: P_loss = -(P_from + P_to), Q_loss = -(Q_from + Q_to)
% This automatically includes line charging and tap changers
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
filename = 'IEEE118_Task1_Results.xlsx';

if exist(filename, 'file') % Delete old file if needed
    delete(filename);
end

% --- Sheet 1: Voltage Bus Data ---
T_voltage = table((1:118)', ...
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

fprintf('Excel export complete: %s\n', filename);

%% 5. Visualization
figure;
plot(res_nr.bus(:, BUS_I), res_nr.bus(:, VM), 'b-s', 'MarkerSize', 4);
grid on; hold on;
yline(0.95, 'r--', 'Lower Limit');
title('Task 1: Base Case Voltage Profile (IEEE 118-Bus)');
xlabel('Bus Number'); ylabel('Voltage Magnitude (p.u.)');