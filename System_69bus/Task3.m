%% Task 3: Base Case Unsymmetrical Fault Analysis
clc; clear; close all;
% addpath(genpath('C:\Program Files\Matlab\matpower8.1'));
% savepath;
define_constants;

%% 1. Load Base Case and Run Power Flowy
mpc = case69;

nb = size(mpc.bus, 1); % Get system size'

mpopt = mpoption('pf.alg','NR','verbose',0,'out.all',0);
res_base = runpf(mpc, mpopt);

%% 2. Extract Prefault Complex Voltages
Vpref = res_base.bus(:, VM) .* exp(1j * deg2rad(res_base.bus(:, VA)));

%% 3. Build Sequence Impedance Matrices (Z-bus)
[Ybus, ~, ~] = makeYbus(mpc);

% --- Ground the Slack Bus ---
% The slack bus is determined by the inbuilt find function
slack_idx = find(mpc.bus(:, BUS_TYPE) == 3);
fprintf('Slack bus is at %d\n', slack_idx);

Ybus_grounded = Ybus;
% We add a realistic source impedance at slack bus to represent the connection to the stiff utility grid/substation.
Y_source = 1 / (0.01 + 1j*0.1); % use a realistic source impedance (e.g., Z = 0.01 + j0.1) (Eyad Zaim Check this)
Ybus_grounded(slack_idx, slack_idx) = Ybus_grounded(slack_idx, slack_idx) + Y_source;

Zbus1 = inv(full(Ybus_grounded)); 
Zbus2 = Zbus1; 
Zbus0 = 3 * Zbus1;

%% 4. Preallocate Result Arrays
SLG_IF = zeros(nb,1); 
LL_IF = zeros(nb,1); 
DLG_IF = zeros(nb,1);

%% 5. Fault Calculation Loop (Symmetrical Components)
for k = 1:nb
    Vf = Vpref(k);
    Z1 = Zbus1(k,k); 
    Z2 = Zbus2(k,k); 
    Z0 = Zbus0(k,k);
    
    % Single Line-to-Ground (SLG)
    SLG_IF(k) = abs(3 * Vf / (Z1 + Z2 + Z0));
    
    % Line-to-Line (LL)
    LL_IF(k)  = abs(sqrt(3) * Vf / (Z1 + Z2));
    
    % Double Line-to-Ground (DLG)
    Zeq = (Z2 * Z0) / (Z2 + Z0);
    DLG_IF(k) = abs(3 * Vf / (Z1 + Zeq));
end

%% 6. ANALYTICAL VALIDATION (SLG, LL, DLG)
% Cross check in excel Sheet 1: Fault Currents
f_bus = 10; % Selection for validation
Vf_v = Vpref(f_bus);
Z1_v = Zbus1(f_bus, f_bus); 
Z2_v = Zbus2(f_bus, f_bus); 
Z0_v = Zbus0(f_bus, f_bus);

% --- SLG Validation ---    
Ia0_slg = Vf_v / (Z1_v + Z2_v + Z0_v);
IF_SLG_v = abs(3 * Ia0_slg);

% --- LL Validation ---
Ia1_ll = Vf_v / (Z1_v + Z2_v);
IF_LL_v = abs(sqrt(3) * Ia1_ll);

% --- DLG Validation ---
Zeq_v = (Z2_v * Z0_v) / (Z2_v + Z0_v); 
Ia1_dlg = Vf_v / (Z1_v + Zeq_v); 
IF_DLG_v = abs(3 * Vf_v / (Z1_v + Zeq_v)); % Total Fault Current flowing to ground (3 * Ia0)

%% 7. POST-FAULT VOLTAGES (During SLG Fault at Bus f_bus)
V_post_fault = zeros(nb, 1);
for i = 1:nb
    Va1 = Vpref(i) - Zbus1(i, f_bus) * Ia0_slg;
    Va2 = 0 - Zbus2(i, f_bus) * Ia0_slg;
    Va0 = 0 - Zbus0(i, f_bus) * Ia0_slg;
    V_post_fault(i) = abs(Va1 + Va2 + Va0);
end

%% 8. Export to Excel
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

filename = sprintf('%s_task3.xlsx', sys_name); % Generates IEEE118_task1.xlsx etc.

if exist(filename, 'file'), delete(filename); end

% Sheet 1: Fault Currents
T_currents = table((1:nb)', abs(Vpref), SLG_IF, LL_IF, DLG_IF, ...
    'VariableNames', {'Bus', 'Prefault_V', 'SLG_Current', 'LL_Current', 'DLG_Current'});
writetable(T_currents, filename, 'Sheet', 'Fault_Current_Results');

% Sheet 2: Post-Fault Voltages
sheet_name = sprintf('PostFault_V_Bus%d', f_bus); % Records the fault bus in the sheet name
T_voltages = table((1:nb)', V_post_fault, ...
    'VariableNames', {'Bus', 'Post_Fault_Voltage_during_SLG_at_BusX'});
writetable(T_voltages, filename, 'Sheet', sheet_name);

%% 9. Console Output & Visualization
fprintf('\n--- ANALYTICAL VALIDATION AT BUS %d ---\n', f_bus);
fprintf('SLG Fault Current: %.4f pu\n', IF_SLG_v);
fprintf('LL  Fault Current: %.4f pu\n', IF_LL_v);
fprintf('DLG Fault Current: %.4f pu\n', IF_DLG_v);
fprintf('Post-Fault Voltage (Phase A) at Bus %d: %.4f pu\n', f_bus, V_post_fault(f_bus));

%% 10. Visualization
fig = figure('Color', 'w', 'Name', [sys_name, ' Fault Analysis']);
ax = axes(fig);
hold(ax, 'on');

% Plotting all three fault types on one plot
plot(ax, 1:nb, SLG_IF, 'r-', 'LineWidth', 1.5, 'DisplayName', 'SLG Fault');
plot(ax, 1:nb, LL_IF, 'b--', 'LineWidth', 1.5, 'DisplayName', 'LL Fault');
plot(ax, 1:nb, DLG_IF, 'k:', 'LineWidth', 1.5, 'DisplayName', 'DLG Fault');

grid(ax, 'on');
xlabel(ax, 'Bus Number');
ylabel(ax, 'Fault Current (p.u.)');
title(ax, ['Unsymmetrical Fault Currents Across ', sys_name, ' System']);
legend(ax, 'Location', 'best');

% Auto-save the plot
ax.Toolbar.Visible = 'off';
plot_filename = sprintf('%s_Task3_Fault.png', sys_name);
exportgraphics(fig, plot_filename, 'Resolution', 300);

fprintf('\nFault plot saved as: %s\n', plot_filename);
fprintf('Task 3 Complete. Results exported to: %s\n', filename);
