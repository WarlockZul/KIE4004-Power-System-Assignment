clc; clear; close all;

% addpath(genpath('C:\Users\User\Downloads\matpower assignment\matpower8.1\matpower8.1'));
% savepath;

define_constants

%% Load IEEE 69-bus system (MATPOWER case)
mpc = case69;

%% ---------- Newton–Raphson ---------- 
mpopt_NR = mpoption('pf.alg', 'NR', ...
                    'verbose', 0, ...
                    'out.all', 0);

%% ---------- Fast Decoupled ----------
mpopt_FD = mpoption('pf.alg', 'FDXB', ...
                    'verbose', 0, ...
                    'out.all', 0);

tic;
res_NR = runpf(mpc, mpopt_NR);
time_NR = toc;

tic;
res_FD = runpf(mpc, mpopt_FD);
time_FD = toc;

fprintf('NR computation time   : %.6f s\n', time_NR);
fprintf('FDLF computation time: %.6f s\n', time_FD);

%% ---------- Voltage Profile ----------
bus = mpc.bus(:, BUS_I);

figure;
plot(bus, res_NR.bus(:, VM), 'b-o','LineWidth',1.5); hold on;
plot(bus, res_FD.bus(:, VM), 'r--s','LineWidth',1.5);
xlabel('Bus Number');
ylabel('Voltage Magnitude (p.u.)');
title('IEEE 69-Bus Voltage Profile (MATPOWER)');
legend('Newton–Raphson','Fast Decoupled');
grid on;

%% ---------- System Losses ----------
P_loss_NR = sum(res_NR.branch(:, PF) + res_NR.branch(:, PT));
P_loss_FD = sum(res_FD.branch(:, PF) + res_FD.branch(:, PT));

fprintf('NR Total Real Power Loss  : %.4f MW\n', P_loss_NR);
fprintf('FDLF Total Real Power Loss: %.4f MW\n', P_loss_FD);

%% ---------- Total System Power Balance ----------
Pgen = sum(res_NR.gen(:, PG));
Pload = sum(res_NR.bus(:, PD));
Ploss = Pgen - Pload;

fprintf('Total Generation : %.4f MW\n', Pgen);
fprintf('Total Load       : %.4f MW\n', Pload);
fprintf('Power Balance    : %.4f MW (loss)\n', Ploss);

%% ---------- Accuracy Comparison ----------
V_diff = abs(res_NR.bus(:, VM) - res_FD.bus(:, VM));
fprintf('Max voltage difference NR vs FDLF: %.6f p.u.\n', max(V_diff));


%% ================= TASK 2: SOLAR PV INTEGRATION =================

%% ---------- Base Case Results (NR) ----------
res_base = res_NR;   % reuse NR base-case result

%% ---------- Solar PV Parameters ----------
PV_bus  = 61;    % location
PV_size = 1.0;   % MW injected

mpc_PV = mpc;

% Model PV as negative load
mpc_PV.bus(PV_bus, PD) = mpc_PV.bus(PV_bus, PD) - PV_size;

%% ---------- Load Flow with PV ----------
res_PV = runpf(mpc_PV, mpopt_NR);

%% ---------- Voltage Profile Comparison ----------
figure;
plot(bus, res_base.bus(:, VM), 'k-o','LineWidth',1.5); hold on;
plot(bus, res_PV.bus(:, VM), 'g-s','LineWidth',1.5);
xlabel('Bus Number');
ylabel('Voltage Magnitude (p.u.)');
title('Impact of Solar PV Integration on IEEE 69-Bus System');
legend('Base Case','With 1 MW Solar PV @ Bus 61');
grid on;

%% ---------- Loss Comparison ----------
loss_base = sum(res_base.branch(:, PF) + res_base.branch(:, PT));
loss_PV   = sum(res_PV.branch(:, PF)   + res_PV.branch(:, PT));

fprintf('\n--- Solar PV Integration Results ---\n');
fprintf('Base Case Loss        : %.4f MW\n', loss_base);
fprintf('With Solar PV Loss    : %.4f MW\n', loss_PV);
fprintf('Loss Reduction        : %.2f %%\n', ...
        (loss_base - loss_PV)/loss_base*100);

%% ---------- Voltage Improvement at Weakest Bus ----------
fprintf('Voltage at Bus %d (Base Case): %.4f p.u.\n', ...
        PV_bus, res_base.bus(PV_bus, VM));
fprintf('Voltage at Bus %d (With PV)  : %.4f p.u.\n', ...
        PV_bus, res_PV.bus(PV_bus, VM));

%% --- OPTIONAL: Algorithm to Find Best Location (Task 2 Justification) ---
results_matrix = zeros(69, 2); % Store [Bus_ID, Total_Loss]

for i = 1:69
    % 1. Reset system to base
    mpc_test = mpc; 
    
    % 2. Place 1 MW PV at Bus 'i'
    mpc_test.bus(i, PD) = mpc_test.bus(i, PD) - 1.0; 
    
    % 3. Run Power Flow
    res_test = runpf(mpc_test, mpopt_NR);
    
    % 4. Calculate and Store Loss
    total_loss = sum(res_test.branch(:, PF) + res_test.branch(:, PT));
    results_matrix(i, :) = [i, total_loss];
end

% Find the minimum loss
[min_loss, index] = min(results_matrix(:, 2));
optimal_bus = results_matrix(index, 1);

fprintf('Automatic Search Complete:\n');
fprintf('The optimal location for 1 MW PV is Bus %d with loss %.4f MW\n', optimal_bus, min_loss);

% Plot the "Sensitivity" Graph
figure;
plot(results_matrix(:,1), results_matrix(:,2), 'b-o');
xlabel('Bus Number');
ylabel('Total System Loss (MW)');
title('Loss vs. PV Location (Justification for Placement)');
grid on;

%% TASK 3: Unsymmetrical Fault Analysis using Symmetrical Components
% Based on Lecture Notes L6 & L7
clear; clc;

%% 1. Load System Data (IEEE 69-Bus)
% Ensure MATPOWER is in path
define_constants;
mpc = case69;

% --- ENGINEERING ASSUMPTIONS (Crucial for Task 3) ---
% L6/L7 notes typically assume Unloaded System for hand calcs (V_pre = 1.0),
% but for a "Project", using the actual pre-fault voltage from Task 1 is better.
USE_ACTUAL_VOLTAGE = true; 

if USE_ACTUAL_VOLTAGE
    % Run base case to get pre-fault voltages
    mpopt = mpoption('pf.alg', 'NR', 'verbose', 0, 'out.all', 0);
    res = runpf(mpc, mpopt);
    V_pre = res.bus(:, VM) .* exp(1j * deg2rad(res.bus(:, VA)));
else
    V_pre = ones(69, 1); % Flat start assumption (1.0 p.u.)
end

%% 2. Construct Sequence Impedance Matrices (Z_bus)
% PROBLEM FIX: We must "ground" the Slack Bus (Bus 1) to create a return path
% for the current. In Thevenin equivalent, voltage sources are shorted to ground.
% We do this by adding a Large Admittance (1e8) to the Slack Bus diagonal.

slack_bus_idx = 1; % IEEE 69-Bus Slack is Bus 1
large_admittance = 1e8; % Represents connection to stiff external grid

% --- Positive Sequence (Z1) ---
[Y_bus1, ~, ~] = makeYbus(mpc);
% "Ground" the slack bus so the matrix is invertible
Y_bus1(slack_bus_idx, slack_bus_idx) = Y_bus1(slack_bus_idx, slack_bus_idx) + large_admittance;
Z_bus1 = inv(full(Y_bus1)); 

% --- Negative Sequence (Z2) ---
% Assumption: Z2 = Z1 for lines
Z_bus2 = Z_bus1;

% --- Zero Sequence (Z0) ---
% Assumption: R0 = 3*R1, X0 = 3*X1
mpc_zero = mpc;
mpc_zero.branch(:, BR_R) = mpc_zero.branch(:, BR_R) * 3; 
mpc_zero.branch(:, BR_X) = mpc_zero.branch(:, BR_X) * 3; 
[Y_bus0, ~, ~] = makeYbus(mpc_zero);
% "Ground" the slack bus in Zero Sequence (Assuming Substation Transformer is Grounded Wye)
Y_bus0(slack_bus_idx, slack_bus_idx) = Y_bus0(slack_bus_idx, slack_bus_idx) + large_admittance;
Z_bus0 = inv(full(Y_bus0));

%% 3. Select Fault Location
Fault_Bus = 61; % You can change this to any bus (e.g., 33, 50)
k = Fault_Bus;  % Index for easy typing

% Get Thevenin Impedances for this specific bus (Diagonal elements)
Z1_th = Z_bus1(k, k);
Z2_th = Z_bus2(k, k);
Z0_th = Z_bus0(k, k);

% Pre-fault Voltage at Fault Bus
V_f = V_pre(k);

% Fault Impedance (Z_f) - Usually 0 for "bolted" faults
Z_f = 0; 

% Symmetrical Component Transformation Matrix (A) [L6 Slide]
a = exp(1j * 120 * pi/180); % 1 at 120 degrees
A_matrix = [1 1 1; 1 a^2 a; 1 a a^2];

fprintf('=================================================\n');
fprintf('TASK 3: UNSYMMETRICAL FAULT ANALYSIS (Bus %d)\n', Fault_Bus);
fprintf('Pre-fault Voltage: %.4f < %.2f deg p.u.\n', abs(V_f), rad2deg(angle(V_f)));
fprintf('=================================================\n\n');

%% 4. Perform Fault Calculations (Formulas from L7)

% --- CASE A: Single Line-to-Ground (SLG) Fault (Phase A) ---
% Formula: Networks in Series
% I0 = I1 = I2 = Vf / (Z1 + Z2 + Z0 + 3*Zf)
I_seq_SLG = V_f / (Z1_th + Z2_th + Z0_th + 3*Z_f);
I012_SLG = [I_seq_SLG; I_seq_SLG; I_seq_SLG];

% Convert to Phase Currents: I_abc = A * I_012
I_abc_SLG = A_matrix * I012_SLG;
If_SLG = abs(I_abc_SLG(1)); % Fault current is Ia

fprintf('--- Type 1: Single Line-to-Ground (SLG) ---\n');
fprintf('Sequence Current (I0=I1=I2): %.4f p.u.\n', abs(I_seq_SLG));
fprintf('Fault Current (Ia)         : %.4f p.u.\n', If_SLG);
fprintf('Phase B Current            : %.4f p.u. (Should be 0)\n', abs(I_abc_SLG(2)));
fprintf('\n');

% --- CASE B: Line-to-Line (LL) Fault (Phase B to C) ---
% Formula: Positive & Negative in Parallel. Zero is inactive.
% I1 = -I2 = Vf / (Z1 + Z2 + Zf)
% I0 = 0
I1_LL = V_f / (Z1_th + Z2_th + Z_f);
I2_LL = -I1_LL;
I0_LL = 0;
I012_LL = [I0_LL; I1_LL; I2_LL];

% Convert to Phase Currents
I_abc_LL = A_matrix * I012_LL;
If_LL = abs(I_abc_LL(2)); % Fault current flows B -> C

fprintf('--- Type 2: Line-to-Line (LL) ---\n');
fprintf('Positive Seq Current (I1)  : %.4f p.u.\n', abs(I1_LL));
fprintf('Fault Current (|Ib|)       : %.4f p.u.\n', If_LL);
fprintf('Phase A Current            : %.4f p.u. (Should be 0)\n', abs(I_abc_LL(1)));
fprintf('\n');

% --- CASE C: Double Line-to-Ground (DLG) Fault (Phase B-C-G) ---
% Formula: Networks in Parallel
% I1 = Vf / (Z1 + (Z2 || (Z0 + 3Zf)))
Z_eq_parallel = (Z2_th * (Z0_th + 3*Z_f)) / (Z2_th + Z0_th + 3*Z_f);
I1_DLG = V_f / (Z1_th + Z_eq_parallel);

% Current Division for I2 and I0
I2_DLG = -I1_DLG * ((Z0_th + 3*Z_f) / (Z2_th + Z0_th + 3*Z_f));
I0_DLG = -I1_DLG * (Z2_th / (Z2_th + Z0_th + 3*Z_f));
I012_DLG = [I0_DLG; I1_DLG; I2_DLG];

% Convert to Phase Currents
I_abc_DLG = A_matrix * I012_DLG;
If_DLG_Total = abs(I_abc_DLG(1) + I_abc_DLG(2) + I_abc_DLG(3)); % Ground current (3I0)

fprintf('--- Type 3: Double Line-to-Ground (DLG) ---\n');
fprintf('Fault Current in Phase B   : %.4f p.u.\n', abs(I_abc_DLG(2)));
fprintf('Fault Current in Phase C   : %.4f p.u.\n', abs(I_abc_DLG(3)));
fprintf('Total Ground Current (If)  : %.4f p.u.\n', abs(3*I0_DLG));
%% ================= TASK 2: SOLAR PV INTEGRATION =================

%% ---------- Base Case Results (NR) ----------
res_base = res_NR;   % reuse NR base-case result

%% ---------- Solar PV Parameters ----------
PV_bus  = 61;    % location
PV_size = 1.0;   % MW injected

mpc_PV = mpc;

% Model PV as negative load
mpc_PV.bus(PV_bus, PD) = mpc_PV.bus(PV_bus, PD) - PV_size;

%% ---------- Load Flow with PV ----------
res_PV = runpf(mpc_PV, mpopt_NR);

%% ---------- Voltage Profile Comparison ----------
figure;
plot(bus, res_base.bus(:, VM), 'k-o','LineWidth',1.5); hold on;
plot(bus, res_PV.bus(:, VM), 'g-s','LineWidth',1.5);
xlabel('Bus Number');
ylabel('Voltage Magnitude (p.u.)');
title('Impact of Solar PV Integration on IEEE 69-Bus System');
legend('Base Case','With 1 MW Solar PV @ Bus 61');
grid on;

%% ---------- Loss Comparison ----------
loss_base = sum(res_base.branch(:, PF) + res_base.branch(:, PT));
loss_PV   = sum(res_PV.branch(:, PF)   + res_PV.branch(:, PT));

fprintf('\n--- Solar PV Integration Results ---\n');
fprintf('Base Case Loss        : %.4f MW\n', loss_base);
fprintf('With Solar PV Loss    : %.4f MW\n', loss_PV);
fprintf('Loss Reduction        : %.2f %%\n', ...
        (loss_base - loss_PV)/loss_base*100);

%% ---------- Voltage Improvement at Weakest Bus ----------
fprintf('Voltage at Bus %d (Base Case): %.4f p.u.\n', ...
        PV_bus, res_base.bus(PV_bus, VM));
fprintf('Voltage at Bus %d (With PV)  : %.4f p.u.\n', ...
        PV_bus, res_PV.bus(PV_bus, VM));

%% --- OPTIONAL: Algorithm to Find Best Location (Task 2 Justification) ---
results_matrix = zeros(69, 2); % Store [Bus_ID, Total_Loss]

for i = 1:69
    % 1. Reset system to base
    mpc_test = mpc; 
    
    % 2. Place 1 MW PV at Bus 'i'
    mpc_test.bus(i, PD) = mpc_test.bus(i, PD) - 1.0; 
    
    % 3. Run Power Flow
    res_test = runpf(mpc_test, mpopt_NR);
    
    % 4. Calculate and Store Loss
    total_loss = sum(res_test.branch(:, PF) + res_test.branch(:, PT));
    results_matrix(i, :) = [i, total_loss];
end

% Find the minimum loss
[min_loss, index] = min(results_matrix(:, 2));
optimal_bus = results_matrix(index, 1);

fprintf('Automatic Search Complete:\n');
fprintf('The optimal location for 1 MW PV is Bus %d with loss %.4f MW\n', optimal_bus, min_loss);

% Plot the "Sensitivity" Graph
figure;
plot(results_matrix(:,1), results_matrix(:,2), 'b-o');
xlabel('Bus Number');
ylabel('Total System Loss (MW)');
title('Loss vs. PV Location (Justification for Placement)');
grid on;

%% TASK 3: Unsymmetrical Fault Analysis using Symmetrical Components
% Based on Lecture Notes L6 & L7
clear; clc;

%% 1. Load System Data (IEEE 69-Bus)
% Ensure MATPOWER is in path
define_constants;
mpc = case69;

% --- ENGINEERING ASSUMPTIONS (Crucial for Task 3) ---
% L6/L7 notes typically assume Unloaded System for hand calcs (V_pre = 1.0),
% but for a "Project", using the actual pre-fault voltage from Task 1 is better.
USE_ACTUAL_VOLTAGE = true; 

if USE_ACTUAL_VOLTAGE
    % Run base case to get pre-fault voltages
    mpopt = mpoption('pf.alg', 'NR', 'verbose', 0, 'out.all', 0);
    res = runpf(mpc, mpopt);
    V_pre = res.bus(:, VM) .* exp(1j * deg2rad(res.bus(:, VA)));
else
    V_pre = ones(69, 1); % Flat start assumption (1.0 p.u.)
end

%% 2. Construct Sequence Impedance Matrices (Z_bus)
% PROBLEM FIX: We must "ground" the Slack Bus (Bus 1) to create a return path
% for the current. In Thevenin equivalent, voltage sources are shorted to ground.
% We do this by adding a Large Admittance (1e8) to the Slack Bus diagonal.

slack_bus_idx = 1; % IEEE 69-Bus Slack is Bus 1
large_admittance = 1e8; % Represents connection to stiff external grid

% --- Positive Sequence (Z1) ---
[Y_bus1, ~, ~] = makeYbus(mpc);
% "Ground" the slack bus so the matrix is invertible
Y_bus1(slack_bus_idx, slack_bus_idx) = Y_bus1(slack_bus_idx, slack_bus_idx) + large_admittance;
Z_bus1 = inv(full(Y_bus1)); 

% --- Negative Sequence (Z2) ---
% Assumption: Z2 = Z1 for lines
Z_bus2 = Z_bus1;

% --- Zero Sequence (Z0) ---
% Assumption: R0 = 3*R1, X0 = 3*X1
mpc_zero = mpc;
mpc_zero.branch(:, BR_R) = mpc_zero.branch(:, BR_R) * 3; 
mpc_zero.branch(:, BR_X) = mpc_zero.branch(:, BR_X) * 3; 
[Y_bus0, ~, ~] = makeYbus(mpc_zero);
% "Ground" the slack bus in Zero Sequence (Assuming Substation Transformer is Grounded Wye)
Y_bus0(slack_bus_idx, slack_bus_idx) = Y_bus0(slack_bus_idx, slack_bus_idx) + large_admittance;
Z_bus0 = inv(full(Y_bus0));

%% 3. Select Fault Location
Fault_Bus = 61; % You can change this to any bus (e.g., 33, 50)
k = Fault_Bus;  % Index for easy typing

% Get Thevenin Impedances for this specific bus (Diagonal elements)
Z1_th = Z_bus1(k, k);
Z2_th = Z_bus2(k, k);
Z0_th = Z_bus0(k, k);

% Pre-fault Voltage at Fault Bus
V_f = V_pre(k);

% Fault Impedance (Z_f) - Usually 0 for "bolted" faults
Z_f = 0; 

% Symmetrical Component Transformation Matrix (A) [L6 Slide]
a = exp(1j * 120 * pi/180); % 1 at 120 degrees
A_matrix = [1 1 1; 1 a^2 a; 1 a a^2];

fprintf('=================================================\n');
fprintf('TASK 3: UNSYMMETRICAL FAULT ANALYSIS (Bus %d)\n', Fault_Bus);
fprintf('Pre-fault Voltage: %.4f < %.2f deg p.u.\n', abs(V_f), rad2deg(angle(V_f)));
fprintf('=================================================\n\n');

%% 4. Perform Fault Calculations (Formulas from L7)

% --- CASE A: Single Line-to-Ground (SLG) Fault (Phase A) ---
% Formula: Networks in Series
% I0 = I1 = I2 = Vf / (Z1 + Z2 + Z0 + 3*Zf)
I_seq_SLG = V_f / (Z1_th + Z2_th + Z0_th + 3*Z_f);
I012_SLG = [I_seq_SLG; I_seq_SLG; I_seq_SLG];

% Convert to Phase Currents: I_abc = A * I_012
I_abc_SLG = A_matrix * I012_SLG;
If_SLG = abs(I_abc_SLG(1)); % Fault current is Ia

fprintf('--- Type 1: Single Line-to-Ground (SLG) ---\n');
fprintf('Sequence Current (I0=I1=I2): %.4f p.u.\n', abs(I_seq_SLG));
fprintf('Fault Current (Ia)         : %.4f p.u.\n', If_SLG);
fprintf('Phase B Current            : %.4f p.u. (Should be 0)\n', abs(I_abc_SLG(2)));
fprintf('\n');

% --- CASE B: Line-to-Line (LL) Fault (Phase B to C) ---
% Formula: Positive & Negative in Parallel. Zero is inactive.
% I1 = -I2 = Vf / (Z1 + Z2 + Zf)
% I0 = 0
I1_LL = V_f / (Z1_th + Z2_th + Z_f);
I2_LL = -I1_LL;
I0_LL = 0;
I012_LL = [I0_LL; I1_LL; I2_LL];

% Convert to Phase Currents
I_abc_LL = A_matrix * I012_LL;
If_LL = abs(I_abc_LL(2)); % Fault current flows B -> C

fprintf('--- Type 2: Line-to-Line (LL) ---\n');
fprintf('Positive Seq Current (I1)  : %.4f p.u.\n', abs(I1_LL));
fprintf('Fault Current (|Ib|)       : %.4f p.u.\n', If_LL);
fprintf('Phase A Current            : %.4f p.u. (Should be 0)\n', abs(I_abc_LL(1)));
fprintf('\n');

% --- CASE C: Double Line-to-Ground (DLG) Fault (Phase B-C-G) ---
% Formula: Networks in Parallel
% I1 = Vf / (Z1 + (Z2 || (Z0 + 3Zf)))
Z_eq_parallel = (Z2_th * (Z0_th + 3*Z_f)) / (Z2_th + Z0_th + 3*Z_f);
I1_DLG = V_f / (Z1_th + Z_eq_parallel);

% Current Division for I2 and I0
I2_DLG = -I1_DLG * ((Z0_th + 3*Z_f) / (Z2_th + Z0_th + 3*Z_f));
I0_DLG = -I1_DLG * (Z2_th / (Z2_th + Z0_th + 3*Z_f));
I012_DLG = [I0_DLG; I1_DLG; I2_DLG];

% Convert to Phase Currents
I_abc_DLG = A_matrix * I012_DLG;
If_DLG_Total = abs(I_abc_DLG(1) + I_abc_DLG(2) + I_abc_DLG(3)); % Ground current (3I0)

fprintf('--- Type 3: Double Line-to-Ground (DLG) ---\n');
fprintf('Fault Current in Phase B   : %.4f p.u.\n', abs(I_abc_DLG(2)));
fprintf('Fault Current in Phase C   : %.4f p.u.\n', abs(I_abc_DLG(3)));
fprintf('Total Ground Current (I