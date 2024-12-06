%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:       Giovanni Fereoli
% DATE:         December 1, 2024
% CLASS:        ASEN 6015: Aerospace Vehicle Guidance and Control
% INSTRUCTOR:   Prof. Jay W. McMahon
% ASSIGNMENT:   Project
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Soft-Landing, 25143 Itokawa
clear; clc; close;

% Small Body Parameters
mu = 2.36 * 1e-9;  % Gravitational parameter (km^3/s^2)
omega = [0; 0; 1.439 * 1e-4];  % Rotation rate (rad/s)

% Initial Spacecraft State
pos = [0; 0; 0.5];  % Position (km)
vel = [-0.01; 0; 0];  % Velocity (km/s)
target_pos = [0; 0; 0.1222];  
target_vel = [0; 0; 0];  
tgo = 50;  

% Simulation Parameters
dt_sim = 1;  
dt_guidance = 1;  
u_max = 0.75 * 1e-3;  % Max acceleration (km/s^2)
time = 0; 
theta_glide = deg2rad(15); 

% Initialize arrays for storage
pos_history = [];
vel_history = [];
a_cmd_history = [];
ZEM_history = [];
ZEV_history = [];
kr_history = [];
kv_history = [];
tgo_history = [];
time_history = [];
cpu_time = [];
a_cmd = [0; 0; 0];  

while tgo > 0
    % Guidance Update
    if mod(time, dt_guidance) == 0 || tgo <= dt_guidance
         % Compute commanded acceleration using predictor-corrector
        t_start = cputime;
        [a_cmd, ZEM, ZEV, tgo, kr, kv] = predictor_corrector_zem_zev(pos,...
            vel, target_pos, target_vel, mu, tgo, u_max, omega, dt_sim, theta_glide);
        t_end = cputime - t_start;
        
        % Print results
        fprintf('Time: %.1f s, Commanded Acceleration: [%f, %f, %f] km/s^2\n', ...
                time, a_cmd(1), a_cmd(2), a_cmd(3));
    end

    % Integrate Dynamics 
    [pos, vel] = integrate_dynamics_ode113(pos, vel, a_cmd, dt_sim, omega, mu);

    % Store Results
    pos_history = [pos_history; pos'];
    vel_history = [vel_history; vel'];
    a_cmd_history = [a_cmd_history; a_cmd'];
    ZEM_history = [ZEM_history; ZEM'];
    ZEV_history = [ZEV_history; ZEV'];
    kr_history = [kr_history; kr];
    kv_history = [kv_history; kv];
    tgo_history = [tgo_history, tgo];
    time_history = [time_history; time];
    cpu_time = [cpu_time; t_end];

    % Update Time
    time = time + dt_sim;
    tgo = max(tgo - dt_sim, 0);  
end

%% Check CPU time and glide-slope constraint

% Print the results of CPU Time
mean_cpu_time = mean(cpu_time);
std_cpu_time = std(cpu_time);
fprintf('Mean CPU time: %.6f seconds\n', mean_cpu_time);
fprintf('Standard deviation of CPU time: %.6f seconds\n', std_cpu_time);

 % Compute glide-slope constraint over the trajectory
r_rel = pos_history - target_pos'; 
r_perp = vecnorm(r_rel(:, 1:2), 2, 2);  % Perpendicular distance 
r_parallel = r_rel(:, 3);               % Parallel distance 
h_glideslope = r_perp - tan(theta_glide) * r_parallel;  
glide_constraints = max(h_glideslope, [], 2);  

%% Plots

% Plot spacecraft trajectory
gca1 = figure(1);
hold on;
[vertices, faces] = readObj('Itokawa.obj');
trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), ...
    'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha',...
    0.8, 'DisplayName', '25142 Itokawa');
plot3(pos_history(:, 1), pos_history(:, 2), pos_history(:, 3),...
    'LineWidth', 3, 'Color', 'k', 'DisplayName', 'Trajectory');
quiver3(pos_history(:, 1), pos_history(:, 2), pos_history(:, 3), ...
    -a_cmd_history(:, 1), -a_cmd_history(:, 2), -a_cmd_history(:, 3), ...
    0.5, 'r', 'LineWidth', 1.2, 'DisplayName', 'Control Direction');

cone_resolution = 50;  
cone_height = max(pos_history(:, 3));  
cone_radius = cone_height * tan(theta_glide);  
[theta, z] = meshgrid(linspace(0, 2 * pi, cone_resolution), linspace(0,...
    cone_height, cone_resolution));
x = z * tan(theta_glide) .* cos(theta);
y = z * tan(theta_glide) .* sin(theta);
x = x + target_pos(1);
y = y + target_pos(2);
z = z + target_pos(3);
surf(x, y, z, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor',...
    [0.1, 0.5, 0.8], 'DisplayName', 'Glide-Slope Cone');

xlabel('$x$ (km)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y$ (km)', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('$z$ (km)', 'Interpreter', 'latex', 'FontSize', 14);
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
grid on;
axis equal;
lighting gouraud;
camlight;
set(gca, 'FontSize', 12);
view(0,0);
exportgraphics(gca1, 'P_1.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot commanded acceleration
gca2 = figure(2);
hold on;
plot(time_history, a_cmd_history(:, 1), 'LineWidth', 1.5, 'Color',...
    [0.8500 0.3250 0.0980], 'DisplayName', '$u_x$');
plot(time_history, a_cmd_history(:, 2), 'LineWidth', 1.5, 'Color',...
    [0.9290 0.6940 0.1250], 'DisplayName', '$u_y$');
plot(time_history, a_cmd_history(:, 3), 'LineWidth', 1.5, 'Color',...
    [0.4940 0.1840 0.5560], 'DisplayName', '$u_z$');
yline(u_max, '--', 'LineWidth', 1.5, 'Color', 'k', 'DisplayName',...
    '$u_{\max}$');
yline(-u_max, '--', 'LineWidth', 1.5, 'Color', 'k', 'DisplayName',...
    '$-u_{\max}$');
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Acceleration $\times 1e-4$ (km/s$^2$)', 'Interpreter', 'latex', 'FontSize', 14);
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);
exportgraphics(gca2, 'P_2.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot state components convergence
gca3 = figure(3);
subplot(2, 1, 1);
hold on;
plot(time_history, pos_history(:, 1), 'LineWidth', 1.5, 'Color',...
    [0.8500 0.3250 0.0980], 'DisplayName', '$x$');
plot(time_history, pos_history(:, 2), 'LineWidth', 1.5, 'Color',...
    [0.9290 0.6940 0.1250], 'DisplayName', '$y$');
plot(time_history, pos_history(:, 3), 'LineWidth', 1.5, 'Color',...
    [0.4940 0.1840 0.5560], 'DisplayName', '$z$');
yline(target_pos(1), '--', 'LineWidth', 1.2, 'Color',...
    [0.8500 0.3250 0.0980], 'DisplayName', '$x_{target}$');
yline(target_pos(2), '--', 'LineWidth', 1.2, 'Color',...
    [0.9290 0.6940 0.1250], 'DisplayName', '$y_{target}$');
yline(target_pos(3), '--', 'LineWidth', 1.2, 'Color',...
    [0.4940 0.1840 0.5560], 'DisplayName', '$z_{target}$');
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Position (km)', 'Interpreter', 'latex', 'FontSize', 14);
title('Position Convergence', 'Interpreter', 'latex', 'FontSize', 16);
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);
subplot(2, 1, 2);
hold on;
plot(time_history, vel_history(:, 1), 'LineWidth', 1.5, 'Color', ...
    [0.8500 0.3250 0.0980], 'DisplayName', '$v_x$');
plot(time_history, vel_history(:, 2), 'LineWidth', 1.5, 'Color', ...
    [0.9290 0.6940 0.1250], 'DisplayName', '$v_y$');
plot(time_history, vel_history(:, 3), 'LineWidth', 1.5, 'Color', ...
    [0.4940 0.1840 0.5560], 'DisplayName', '$v_z$');
yline(target_vel(1), '--', 'LineWidth', 1.2, 'Color', ...
    [0.8500 0.3250 0.0980], 'DisplayName', '$v_{x,target}$');
yline(target_vel(2), '--', 'LineWidth', 1.2, 'Color', ...
    [0.9290 0.6940 0.1250], 'DisplayName', '$v_{y,target}$');
yline(target_vel(3), '--', 'LineWidth', 1.2, 'Color', ...
    [0.4940 0.1840 0.5560], 'DisplayName', '$v_{z,target}$');
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Velocity (km/s)', 'Interpreter', 'latex', 'FontSize', 14);
title('Velocity Convergence', 'Interpreter', 'latex', 'FontSize', 16);
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);
exportgraphics(gca3, 'P_3.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot ZEM/ZEV quantities
gca4 = figure(4);
subplot(2, 1, 1);
hold on;
plot(time_history, ZEM_history(:, 1), 'LineWidth', 1.5, 'Color', ...
    [0.8500 0.3250 0.0980], 'DisplayName', '$ZEM_x$');
plot(time_history, ZEM_history(:, 2), 'LineWidth', 1.5, 'Color', ...
    [0.9290 0.6940 0.1250], 'DisplayName', '$ZEM_y$');
plot(time_history, ZEM_history(:, 3), 'LineWidth', 1.5, 'Color', ...
    [0.4940 0.1840 0.5560], 'DisplayName', '$ZEM_z$');
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('ZEM (km)', 'Interpreter', 'latex', 'FontSize', 14);
title('ZEM Evolution', 'Interpreter', 'latex', 'FontSize', 16);
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);
subplot(2, 1, 2);
hold on;
plot(time_history, ZEV_history(:, 1), 'LineWidth', 1.5, 'Color',...
    [0.8500 0.3250 0.0980], 'DisplayName', '$ZEV_x$');
plot(time_history, ZEV_history(:, 2), 'LineWidth', 1.5, 'Color', ...
    [0.9290 0.6940 0.1250], 'DisplayName', '$ZEV_y$');
plot(time_history, ZEV_history(:, 3), 'LineWidth', 1.5, 'Color', ...
    [0.4940 0.1840 0.5560], 'DisplayName', '$ZEV_z$');
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('ZEV (km/s)', 'Interpreter', 'latex', 'FontSize', 14);
title('ZEV Evolution', 'Interpreter', 'latex', 'FontSize', 16);
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);
exportgraphics(gca4, 'P_4.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot kr and kv histories
gca5 = figure(5);
subplot(2, 1, 1);
hold on;
plot(time_history, kr_history, 'LineWidth', 1.5, 'Color', ...
    [0.8500 0.3250 0.0980], 'DisplayName', '$k_r$');
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$k_r$ (-)', 'Interpreter', 'latex', 'FontSize', 14);
title('$k_r$ History', 'Interpreter', 'latex', 'FontSize', 16);
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
yticks(floor(min(kr_history)):1:ceil(max(kr_history)));
grid on;
set(gca, 'FontSize', 12);
subplot(2, 1, 2);
hold on;
plot(time_history, kv_history, 'LineWidth', 1.5, 'Color', ...
    [0.4940 0.1840 0.5560], 'DisplayName', '$k_v$');
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$k_v$ (-)', 'Interpreter', 'latex', 'FontSize', 14);
title('$k_v$ History', 'Interpreter', 'latex', 'FontSize', 16);
legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
yticks(floor(min(kv_history)):1:ceil(max(kv_history)));
grid on;
set(gca, 'FontSize', 12);
exportgraphics(gca5, 'P_5.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Monte Carlo Simulation Parameters
clear; clc; close;

% Small Body Parameters
mu = 2.36 * 1e-9;  % Gravitational parameter (km^3/s^2)
omega = [0; 0; 1.439 * 1e-4];  % Rotation rate (rad/s)

% Initial Spacecraft State
pos = [0; 0; 0.5];  % Position (km)
vel = [-0.01; 0; 0];  % Velocity (km/s)
target_pos = [0; 0; 0.1222];  
target_vel = [0; 0; 0];  
tgo = 50;  

% Simulation Parameters
dt_sim = 1;  
dt_guidance = 1;  
u_max = 0.75 * 1e-3;  % Max acceleration (km/s^2)
tgo_min = 0; 
time = 0; 
theta_glide = deg2rad(15); 
num_simulations = 50;  % Number of Monte Carlo runs
pos_errors = zeros(num_simulations, 3);  % Final position errors
vel_errors = zeros(num_simulations, 3);  % Final velocity errors
all_trajectories = cell(num_simulations, 1);  % Store trajectories

rng(1);  % Set random seed for reproducibility

for i = 1:num_simulations
    % Perturb initial conditions
    pos_perturb = pos + norm(pos) / 100 * randn(3, 1);  % Add random noise to position
    vel_perturb = vel + norm(vel) / 100 * randn(3, 1); % Add random noise to velocity

    % Initialize variables
    pos_mc = pos_perturb;
    vel_mc = vel_perturb;
    tgo_mc = tgo;
    time_mc = 0;
    trajectory = [];
    a_cmd = [0; 0; 0];
    
    while tgo_mc > 0
        % Guidance Update
        if mod(time_mc, dt_guidance) == 0 || tgo_mc <= dt_guidance
            [a_cmd, ~, ~, tgo_mc, ~, ~] = predictor_corrector_zem_zev( ...
                pos_mc, vel_mc, target_pos, target_vel, mu, tgo_mc, u_max, omega, dt_sim, theta_glide);
        end

        % Apply 3% perturbation to a_cmd
        perturbation_factor = 1 + 0.03 * (2 * rand(3, 1) - 1);
        a_cmd_perturbed = a_cmd .* perturbation_factor;

        % Integrate Dynamics
        [pos_mc, vel_mc] = integrate_dynamics_ode113(pos_mc, vel_mc, a_cmd, dt_sim, omega, mu);

        % Store trajectory
        trajectory = [trajectory; pos_mc'];

        % Update Time
        time_mc = time_mc + dt_sim;
        tgo_mc = max(tgo_mc - dt_sim, 0);
    end

    % Store final position and velocity errors
    pos_errors(i, :) = pos_mc - target_pos;
    vel_errors(i, :) = vel_mc - target_vel;

    % Store trajectory
    all_trajectories{i} = trajectory;
end

%% Plots MC

% Compute norms of position and velocity errors
pos_norms = vecnorm(pos_errors, 2, 2);  % Norm of position errors
vel_norms = vecnorm(vel_errors, 2, 2);  % Norm of velocity errors

% Plot Trajectories
gca10 = figure(10);
hold on;
[vertices, faces] = readObj('Itokawa.obj');
trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), ...
    'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha',...
    0.8);

cone_resolution = 50;  
cone_height = max(all_trajectories{1}(:, 3));  
cone_radius = cone_height * tan(theta_glide);  
[theta, z] = meshgrid(linspace(0, 2 * pi, cone_resolution), linspace(0,...
    cone_height, cone_resolution));
x = z * tan(theta_glide) .* cos(theta);
y = z * tan(theta_glide) .* sin(theta);
x = x + target_pos(1);
y = y + target_pos(2);
z = z + target_pos(3);
surf(x, y, z, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor',...
    [0.1, 0.5, 0.8], 'DisplayName', 'Glide-Slope Cone');
for i = 1:num_simulations
    traj = all_trajectories{i};
    plot3(traj(:, 1), traj(:, 2), traj(:, 3), 'LineWidth', 2, 'Color', ...
        rand(1, 3));  % Random colors for each trajectory    
end
xlabel('$x$ (km)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y$ (km)', 'Interpreter', 'latex', 'FontSize', 14);
zlabel('$z$ (km)', 'Interpreter', 'latex', 'FontSize', 14);
legend({'25142 Itokawa', 'Glide-Slope Cone', 'Trajectory'}, 'Interpreter', ...
    'latex', 'Location', 'best', 'FontSize', 12);
grid on;
axis equal;
lighting gouraud;
camlight;
set(gca, 'FontSize', 12);
view(0, 0);
exportgraphics(gca10, 'P_10.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Position Error Norm Histogram
gca11 = figure(11);
hold on;
histogram(pos_norms, num_simulations / 2, 'Normalization', 'pdf', 'FaceColor', ...
    [0.2, 0.6, 0.8], 'DisplayName', 'Norm of Position Errors');
pos_mean_norm = mean(pos_norms);
pos_std_norm = std(pos_norms);
x_range = linspace(min(pos_norms), max(pos_norms), 100);
y_pdf = normpdf(x_range, pos_mean_norm, pos_std_norm);
plot(x_range, y_pdf, 'LineWidth', 2, 'Color', [0.8, 0.2, 0.2], ...
    'DisplayName', sprintf('Gaussian Fit\n$\\mu=%.2e$, $\\sigma=%.2e$', ...
    pos_mean_norm, pos_std_norm));
xlabel('$||\delta \mathbf{r}_f|| \times 1e-8$ (km)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('PDF $\times 1e8$  (-)', 'Interpreter', 'latex', 'FontSize', 14);
legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast');
grid on;
set(gca, 'FontSize', 12);
exportgraphics(gca11, 'P_11.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Velocity Error Norm Histogram
gca12 = figure(12);
hold on;
histogram(vel_norms, num_simulations / 2, 'Normalization', 'pdf', 'FaceColor', ...
    [0.3, 0.8, 0.4], 'DisplayName', 'Norm of Velocity Errors');
vel_mean_norm = mean(vel_norms);
vel_std_norm = std(vel_norms);
x_range = linspace(min(vel_norms), max(vel_norms), 100);
y_pdf = normpdf(x_range, vel_mean_norm, vel_std_norm);
plot(x_range, y_pdf, 'LineWidth', 2, 'Color', [0.8, 0.2, 0.2], ...
    'DisplayName', sprintf('Gaussian Fit\n$\\mu=%.2e$, $\\sigma=%.2e$', ...
    vel_mean_norm, vel_std_norm));
xlabel('$||\delta \mathbf{v}_f|| \times 1e-8$ (km/s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('PDF $\times 1e7$  (-)', 'Interpreter', 'latex', 'FontSize', 14);
legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast');
grid on;
set(gca, 'FontSize', 12);
exportgraphics(gca12, 'P_12.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Functions

function [a_cmd, ZEM, ZEV, tf_final, kr_final, kv_final] = predictor_corrector_zem_zev(r0, ...
    v0, r_target, v_target, mu, tgo_initial, u_max, omega, dt_sim, theta_glide)
    % Initial guess for optimization [kr, kv, tf]
    x0 = [1; 1; tgo_initial];
    
    % Optimization bounds
    lb = [-100; -100; 0.9 * tgo_initial]; 
    ub = [100; 100; 1.1 * tgo_initial];  

    % Extended Optimization options for fmincon
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'final', ... 
        'MaxIterations', 50, ...
        'StepTolerance', 1e-14, ...
        'ConstraintTolerance', 1e-3,...  % 'OptimalityTolerance', 1e-6, ...
        'EnableFeasibilityMode', true, ...
        'SubproblemAlgorithm', 'cg',...
        'MaxFunctionEvaluations', 500,...
        'UseParallel', true);

    try
        % Perform a single iteration with fmincon
        [x_opt, J_opt] = fmincon(@(x) cost_function_zem_zev(x, r0, v0, ...
            r_target, v_target, mu, omega, dt_sim), ...
                                 x0, [], [], [], [], lb, ub, ...
                                 @(x) constraints_zem_zev(x, r0, v0, ...
                                 r_target, v_target, mu, omega, u_max, ...
                                 dt_sim, theta_glide), ...
                                 options);

        % Extract optimized parameters
        kr = x_opt(1);
        kv = x_opt(2);
        tf = x_opt(3);

        % Compute ZEM, ZEV, and commanded acceleration
        [ZEM, ZEV, a_cmd] = compute_zem_zev(r0, v0, r_target, v_target, ...
            kr, kv, tf, mu, omega);

        % Clipping acceleration, if needed
        a_cmd = max(min(a_cmd, u_max), -u_max);

        % Output results
        fprintf('Optimization complete. Cost: %.6e\n', J_opt);
        fprintf('kr: %.4f, kv: %.4f, tf: %.4f\n', kr, kv, tf);
        
    catch ME
        % If optimization fails, return the initial guess
        warning('Optimization failed. Returning initial guess.');
        kr = x0(1);
        kv = x0(2);
        tf = x0(3);

        % Compute ZEM, ZEV, and commanded acceleration using initial guess
        [ZEM, ZEV, a_cmd] = compute_zem_zev(r0, v0, r_target, v_target, ...
            kr, kv, tf, mu, omega);
    end

    % Final outputs
    kr_final = kr;
    kv_final = kv;
    tf_final = tf;
end

% Cost Function
function J = cost_function_zem_zev(x, r0, v0, r_target, v_target, mu, ...
    omega, dt_sim)
    % Extract optimization variables
    kr = x(1);
    kv = x(2);
    tf = x(3);

    % Initialize state
    time = 0;
    tgo = tf;
    pos = r0;
    vel = v0;

    % Initialize total commanded acceleration sum
    total_a_cmd = 0;

    % Simulate dynamics over the trajectory
    while tgo > 0
        % Compute commanded acceleration
        [~, ~, a_cmd] = compute_zem_zev(pos, vel, r_target, v_target, ...
            kr, kv, tgo, mu, omega);

        % Accumulate the norm of commanded acceleration
        total_a_cmd = total_a_cmd + norm(a_cmd) * dt_sim;

        % Integrate dynamics
        [pos, vel] = integrate_dynamics_ode113(pos, vel, a_cmd, dt_sim, ...
            omega, mu);

        % Update time and time-to-go
        time = time + dt_sim;
        tgo = max(tgo - dt_sim, 0);
    end

    % Cost function is the total sum of commanded accelerations
    J = total_a_cmd;
end

function [c, ceq] = constraints_zem_zev(x, r0, v0, r_target, v_target, mu, ...
    omega, u_max, dt_sim, theta_glide)
    % Extract optimization variables
    kr = x(1);
    kv = x(2);
    tf = x(3);

    % Initialize state
    time = 0;
    tgo = tf;
    pos = r0;
    vel = v0;

    % Initialize arrays to track constraints
    a_cmd_components = [];
    pos_history = [];
    vel_history = [];

    % Simulate dynamics over the trajectory
    while tgo > 0
        % Compute commanded acceleration
        [~, ~, a_cmd] = compute_zem_zev(pos, vel, r_target, v_target, kr,...
            kv, tgo, mu, omega);

        % Store the magnitude of commanded acceleration
        a_cmd_components = [a_cmd_components, a_cmd];

        % Integrate dynamics
        [pos, vel] = integrate_dynamics_ode113(pos, vel, a_cmd, dt_sim,...
            omega, mu);

        % Update time and time-to-go
        time = time + dt_sim;
        tgo = max(tgo - dt_sim, 0);

        % Store position history
        pos_history = [pos_history, pos];
        vel_history  = [vel_history, vel];
    end

    % Final state from the simulation
    r_final = pos;
    v_final = vel;

    % Compute glide-slope constraint over the trajectory
    r_rel = pos_history - r_target; 
    r_perp = vecnorm(r_rel(1:2, :), 2, 1);  % Perpendicular distance 
    r_parallel = r_rel(3, :);               % Parallel distance 
    h_glideslope = r_perp - tan(theta_glide) * r_parallel;  
    glide_constraints = max(h_glideslope, [], 2);  

    % Check acceleration component constraints
    a_cmd_x = a_cmd_components(1, :);
    a_cmd_y = a_cmd_components(2, :);
    a_cmd_z = a_cmd_components(3, :);

    % Inequality constraints (c <= 0)
    c = [
         max(abs(a_cmd_x)) - u_max;  
         max(abs(a_cmd_y)) - u_max;  
         max(abs(a_cmd_z)) - u_max;  
         glide_constraints; 
    ];

    % Equality constraints (ceq = 0)
    ceq = [
        r_final - r_target; 
        v_final - v_target   
    ];
end

% Commanded Acceleration Calculation (ZEM/ZEV simplified for speed)
function [ZEM, ZEV, a_cmd] = compute_zem_zev(r0, v0, r_target, v_target, ...
    kr, kv, tf, mu, omega)
    % Initialization
    tgo = tf;
    a0 = dynamics_rotating_frame(0, [r0; v0], omega, mu, [0; 0; 0]);
    a0 = a0(4:6);

    % Compute Zero Effort Miss (ZEM)
    ZEM = r_target - r0 - v0 * tgo - 0.5 * a0 * tgo^2;

    % Compute Zero Effort Velocity (ZEV)
    ZEV = v_target - v0 - a0 * tgo;

    % Compute commanded acceleration
    a_cmd = (6 * kr / tgo^2) * ZEM - (2 * kv / tgo) * ZEV;
end

function [pos, vel] = integrate_dynamics_ode113(pos, vel, a_cmd, dt, omega, mu)
    % Initial state vector [position; velocity]
    state_init = [pos; vel];

    % ODE solver options
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

    % Integrate over the timestep
    [~, state] = ode113(@(t, state) dynamics_rotating_frame(t, state, ...
        omega, mu, a_cmd), [0, dt], state_init, options);

    % Extract updated position and velocity
    pos = state(end, 1:3)';
    vel = state(end, 4:6)';
end

function dstatedt = dynamics_rotating_frame(~, state, omega, mu, a_cmd)
    % Extract position and velocity
    pos = state(1:3);
    vel = state(4:6);

    % Compute gravitational acceleration
    r = norm(pos);
    acc_gravity = -mu / r^3 * pos;

    % Compute Coriolis and centrifugal accelerations
    acc_coriolis = - 2 * cross(omega, vel);
    acc_centrifugal = - cross(omega, cross(omega, pos));

    % Total acceleration
    acc_total = acc_gravity + acc_coriolis + acc_centrifugal + a_cmd;

    % Time derivative of the state vector
    dstatedt = [vel; acc_total];
end

function [vertices, faces] = readObj(filename)
    vertices = [];
    faces = [];
    fid = fopen(filename, 'r');
    if fid == -1, error('Cannot open the file.'); end
    while ~feof(fid)
        line = fgetl(fid);
        if startsWith(line, 'v ')
            vertex = sscanf(line, 'v %f %f %f');
            vertices = [vertices; vertex'];
        elseif startsWith(line, 'f ')
            face = sscanf(line, 'f %d %d %d');
            faces = [faces; face'];
        end
    end
    fclose(fid);
end


