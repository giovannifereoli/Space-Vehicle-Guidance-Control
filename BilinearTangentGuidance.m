%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:       Giovanni Fereoli
% DATE:         November 11, 2024
% CLASS:        ASEN 6015: Aerospace Vehicle Guidance and Control
% INSTRUCTOR:   Prof. Jay W. McMahon
% ASSIGNMENT:   Homework 4 - Pt2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Controlled Trajectory
clear; clc; close;

% Constants
mu = 398600;  % [km^3/s^2]
t0 = 0;  
tf_guess = 23;  
thrust_acc = 30 / 1000;  % [km/s^2]
guidance_frequency = 10;  % [Hz]

% Target orbital elements
a_target = 10000;      
e_target = 0.0001;     
i_target = deg2rad(0);  
Omega_target = deg2rad(0);
omega_target = deg2rad(35); 
nu_target = deg2rad(0);    

% Initial state vector (inertial frame)
r0 = [8276; 5612; 5];  
v0 = [-3.142; 4.672; 0]; 
state0 = [r0; v0];

% Compute final target state 
[r_f, v_f] = orbital_elements_to_cartesian(mu, a_target, e_target,...
    i_target, Omega_target, omega_target, nu_target);

% Initialize guidance parameters (LVLH frame)
A = [1; 0; 0]; 
B = [0; 0; 0];
tf = tf_guess;

% Initialize storage 
orbital_elements_history = [];
guidance_history = [];
time_history = [];
error_history = [];

% Online BLT Guidance 
dt_guidance = 1 / guidance_frequency;  % Time step for guidance updates
time = t0;
state = state0;
trajectory = [];
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

while time < tf
    % Integrate for a small time step (up to the next guidance update)
    [t_segment, state_segment] = ode113(@(t, y) equations_of_motion(t,...
        y, A, B, thrust_acc, mu), [time, min(time + dt_guidance,...
        tf)], state, options);
    
    % Append the trajectory
    trajectory = [trajectory; state_segment];
    time = t_segment(end);  % Update time
    state = state_segment(end, :)';  % Update state

    % Update t0 
    t0 = time;  

    % Compute current orbital elements
    r_current = state(1:3);
    v_current = state(4:6);
    [a, e, i, omega, nu] = cartesian_to_orbital_elements(r_current,...
        v_current, mu);
    orbital_elements_history = [orbital_elements_history; a, e,...
        rad2deg(i), rad2deg(omega)];
    time_history = [time_history; time];
    
    % Update guidance parameters
    [A, B, tf, ~, error_final] = predictor_corrector(state, A, B,...
        thrust_acc, mu, a_target, e_target, t0, tf, omega_target);
    guidance_history = [guidance_history; A', B'];
    error_history = [error_history, error_final];
end

% Print final true anomaly and final burn time
fprintf('Final True Anomaly (nu): %.14f deg\n', rad2deg(nu));
fprintf('Final Burn Time (tf): %.14f seconds\n', tf);

% Extract trajectory
r_controlled = trajectory(:, 1:3);

% Integrate initial state without control
t_span1 = [0, 23];  
[~, state_free1] = ode113(@(t, y) two_body_eom(t, y, mu),...
    t_span1, [r0; v0], options);
t_span2 = [0, -23]; 
[~, state_free2] = ode113(@(t, y) two_body_eom(t, y, mu),...
    t_span2, [r0; v0], options);
state_free = [flipud(state_free2); state_free1];

% Integrate target orbit without control from (r_f, v_f)
[~, state_target1] = ode113(@(t, y) two_body_eom(t, y, mu), t_span1,...
    [r_f; v_f], options);
[~, state_target2] = ode113(@(t, y) two_body_eom(t, y, mu), t_span2,...
    [r_f; v_f], options);
state_target = [flipud(state_target2); state_target1];

% Plot Trajectories
gca1 = figure(1);
plot3(state_free(:, 1), state_free(:, 2), state_free(:, 3),...
    'm--', 'LineWidth', 1.5, 'DisplayName', 'Initial Orbit');
hold on;
plot3(state_target(:, 1), state_target(:, 2), state_target(:, 3),...
    'g--', 'LineWidth', 1.5, 'DisplayName', 'Target Orbit');
plot3(r_controlled(:, 1), r_controlled(:, 2), r_controlled(:, 3), 'r',...
    'LineWidth', 2, 'DisplayName', 'Controlled Trajectory');
scatter3(r0(1), r0(2), r0(3), 'bo', 'DisplayName',...
    'Initial Position ($\mathbf{r}_0$)', 'LineWidth', 1.5);
scatter3(r_controlled(end, 1), r_controlled(end, 2),...
    r_controlled(end, 3), 'ro', 'DisplayName', 'Final Position',...
    'LineWidth', 1.5);
xlabel('$x$ [km]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$y$ [km]', 'Interpreter', 'latex', 'FontSize', 12);
zlabel('$z$ [km]', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
legend('Location', 'Best', 'Interpreter', 'latex');
exportgraphics(gca1, 'HW4_1.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Initialize error storage
num_steps = size(orbital_elements_history, 1);
orbital_elements_error = zeros(num_steps, 4); 

for step = 1:num_steps
    % Extract current orbital elements
    a_current = orbital_elements_history(step, 1);  
    e_current = orbital_elements_history(step, 2);  
    i_current = deg2rad(orbital_elements_history(step, 3)); 
    omega_current = deg2rad(orbital_elements_history(step, 4)); 

    % Compute errors
    a_error = a_current - a_target;         
    e_error = e_current - e_target;        
    i_error = i_current - i_target;        
    omega_error = omega_current - omega_target;
    
    % Store in error matrix
    orbital_elements_error(step, :) = [a_error, e_error, ...
        rad2deg(i_error), rad2deg(omega_error)];
end

% Plot Orbital Elements Evolution and Errors
gca2 = figure(2);
colors = lines(4);
subplot(4, 2, 1); 
plot(time_history, orbital_elements_history(:, 1), 'LineWidth',...
    2, 'Color', colors(1, :));
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$a$ [km]', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
title('Semi-Major Axis Evolution', 'Interpreter', 'latex', 'FontSize', 14);
subplot(4, 2, 3); 
plot(time_history, orbital_elements_history(:, 2), 'LineWidth',...
    2, 'Color', colors(2, :));
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$e$ [-]', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
title('Eccentricity Evolution', 'Interpreter', 'latex', 'FontSize', 14);
subplot(4, 2, 5);
plot(time_history, orbital_elements_history(:, 3), 'LineWidth', 2, ...
    'Color', colors(3, :));
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$i$ [deg]', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
title('Inclination Evolution', 'Interpreter', 'latex', 'FontSize', 14);
subplot(4, 2, 7); 
plot(time_history, orbital_elements_history(:, 4), 'LineWidth', 2,...
    'Color', colors(4, :));
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\omega$ [deg]', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
title('Argument of Perigee Evolution', 'Interpreter', 'latex', 'FontSize', 14);

% Right Side: Error Compared to Target Values
subplot(4, 2, 2);
semilogy(time_history, abs(orbital_elements_error(:, 1)), 'LineWidth',...
    2, 'Color', colors(1, :));
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\Delta a$ [km]', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
title('Semi-Major Axis Error', 'Interpreter', 'latex', 'FontSize', 14);
subplot(4, 2, 4); 
semilogy(time_history, abs(orbital_elements_error(:, 2)), 'LineWidth', ...
    2, 'Color', colors(2, :));
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\Delta e$ [-]', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
title('Eccentricity Error', 'Interpreter', 'latex', 'FontSize', 14);
subplot(4, 2, 6); 
semilogy(time_history, abs(orbital_elements_error(:, 3)), 'LineWidth', 2,...
    'Color', colors(3, :));
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\Delta i$ [deg]', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
title('Inclination Error', 'Interpreter', 'latex', 'FontSize', 14);
subplot(4, 2, 8); 
semilogy(time_history, abs(orbital_elements_error(:, 4)), 'LineWidth', 2,...
    'Color', colors(4, :));
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\Delta \omega$ [deg]', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
title('Argument of Perigee Error', 'Interpreter', 'latex', 'FontSize', 14);
set(gcf, 'Position', [200, 200, 800, 800]);
exportgraphics(gca2, 'HW4_2.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot Guidance Parameters Evolution
gca3 = figure(3);
subplot(2, 1, 1);
plot(time_history, guidance_history(:, 1:3), 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('A Components [-]', 'Interpreter', 'latex', 'FontSize', 12);
legend('$A_x$', '$A_y$', '$A_z$', 'Interpreter', 'latex', 'Location',...
    'Best');
grid on;
subplot(2, 1, 2);
plot(time_history, guidance_history(:, 4:6), 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('B Components [-]', 'Interpreter', 'latex', 'FontSize', 12);
legend('$B_x$', '$B_y$', '$B_z$', 'Interpreter', 'latex', 'Location', ...
    'Best');
grid on;
exportgraphics(gca3, 'HW4_3.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot Error History
gca4 = figure(4);
components = {'$|\Delta v_x|$', '$|\Delta v_y|$', '$|\Delta v_z|$', ...
              '$|\Delta r|$', '$|\Delta h|$', '$|\Delta g|$',...
              '$|\Delta \lambda|$'};
colors = lines(7); 
for i = 1:7
    subplot(7, 1, i);
    semilogy(time_history, max(1e-16, abs(error_history(i, :))),...
        'LineWidth', 2, 'Color', colors(i, :));
    xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(components{i}, 'Interpreter', 'latex', 'FontSize', 12);
    grid on;
    title(['Error History: ', components{i}], 'Interpreter',...
        'latex', 'FontSize', 14);
end
set(gcf, 'Position', [200, 200, 800, 800]);
exportgraphics(gca4, 'HW4_4.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Uncontrolled Trajectory
clear; clc; close;

% Constants
mu = 398600;  % [km^3/s^2]
t0 = 0;  
tf = 23;  
thrust_acc = 30 / 1000;  % [km/s^2]

% Target orbital elements
a_target = 10000;  
e_target = 0.0001;     
i_target = deg2rad(0); 
Omega_target = deg2rad(0); 
omega_target = deg2rad(35);
nu_target = deg2rad(0);    

% Initial state vector 
r0 = [8276; 5612; 5]; 
v0 = [-3.142; 4.672; 0];
state0 = [r0; v0];

% Initialize guidance parameters 
A = [1; 0; 0];  
B = [0; 0; 0];  

% Time span for simulation
t_span = [t0, tf];
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Simulate trajectory
[t, trajectory] = ode113(@(t, y) equations_of_motion(t, y,...
    A, B, thrust_acc, mu), t_span, state0, options);

% Compute orbital elements over time
num_steps = length(t);
orbital_elements_history = zeros(num_steps, 4); 
for step = 1:num_steps
    r = trajectory(step, 1:3)';
    v = trajectory(step, 4:6)';
    [a, e, i, omega] = cartesian_to_orbital_elements(r, v, mu);
    orbital_elements_history(step, :) = [a, e, rad2deg(i), rad2deg(omega)];
end

% Final orbital elements
final_orbital_elements = orbital_elements_history(end, :);

% Compute errors
errors = [
    final_orbital_elements(1) - a_target;  % Semi-major axis error
    final_orbital_elements(2) - e_target;  % Eccentricity error
    final_orbital_elements(3) - rad2deg(i_target);  % Inclination error
    final_orbital_elements(4) - rad2deg(omega_target)  % Argument of periapsis error
];

% Print achieved and errors with 8 decimal places
fprintf('\nFinal Orbital Elements and Errors:\n');
fprintf('-----------------------------------\n');
fprintf('Semi-major axis (a): Achieved = %.8f km, Error = %.8f km\n', ...
    final_orbital_elements(1), errors(1));
fprintf('Eccentricity (e):    Achieved = %.8f, Error = %.8f\n', ...
    final_orbital_elements(2), errors(2));
fprintf('Inclination (i):     Achieved = %.8f deg, Error = %.8f deg\n', ...
    final_orbital_elements(3), errors(3));
fprintf('Argument of periapsis (ω): Achieved = %.8f deg, Error = %.8f deg\n', ...
    final_orbital_elements(4), errors(4));


% Integrate initial and target orbits
t_span1 = [0, 23];  
[~, state_free1] = ode113(@(t, y) two_body_eom(t, y, mu), ...
    t_span1, [r0; v0], options);
t_span2 = [0, -23]; 
[~, state_free2] = ode113(@(t, y) two_body_eom(t, y, mu), t_span2,...
    [r0; v0], options);
state_free = [flipud(state_free2); state_free1];
[r_f, v_f] = orbital_elements_to_cartesian(mu, a_target, e_target,...
    i_target, Omega_target, omega_target, nu_target);
[~, state_target1] = ode113(@(t, y) two_body_eom(t, y, mu), t_span1,...
    [r_f; v_f], options);
[~, state_target2] = ode113(@(t, y) two_body_eom(t, y, mu), t_span2,...
    [r_f; v_f], options);
state_target = [flipud(state_target2); state_target1];

% Plot Trajectories
gca1 = figure(1);
plot3(state_free(:, 1), state_free(:, 2), state_free(:, 3), 'm--', ...
    'LineWidth', 1.5, 'DisplayName', 'Initial Orbit');
hold on;
plot3(state_target(:, 1), state_target(:, 2), state_target(:, 3), 'g--',...
    'LineWidth', 1.5, 'DisplayName', 'Target Orbit');
plot3(trajectory(:, 1), trajectory(:, 2), trajectory(:, 3), 'r', ...
    'LineWidth', 2, 'DisplayName', 'Simulated Trajectory');
scatter3(r0(1), r0(2), r0(3), 'bo', 'DisplayName', ...
    'Initial Position ($\mathbf{r}_0$)', 'LineWidth', 1.5);
scatter3(trajectory(end, 1), trajectory(end, 2), trajectory(end, 3),...
    'ro', 'DisplayName', 'Final Position', 'LineWidth', 1.5);
xlabel('$x$ [km]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$y$ [km]', 'Interpreter', 'latex', 'FontSize', 12);
zlabel('$z$ [km]', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
legend('Location', 'Best', 'Interpreter', 'latex');
exportgraphics(gca1, 'HW4_5.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot Orbital Elements Evolution
gca2 = figure(2);
colors = lines(4);
subplot(4, 1, 1);
plot(t, orbital_elements_history(:, 1), 'LineWidth', 2, 'Color', ...
    colors(1, :));
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$a$ [km]', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
title('Semi-Major Axis Evolution', 'Interpreter', 'latex', 'FontSize', 14);
subplot(4, 1, 2);
plot(t, orbital_elements_history(:, 2), 'LineWidth', 2, 'Color',...
    colors(2, :));
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$e$ [-]', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
title('Eccentricity Evolution', 'Interpreter', 'latex', 'FontSize', 14);
subplot(4, 1, 3);
plot(t, orbital_elements_history(:, 3), 'LineWidth', 2, 'Color', ...
    colors(3, :));
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$i$ [deg]', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
title('Inclination Evolution', 'Interpreter', 'latex', 'FontSize', 14);
subplot(4, 1, 4);
plot(t, orbital_elements_history(:, 4), 'LineWidth', 2, 'Color',...
    colors(4, :));
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\omega$ [deg]', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
title('Argument of Perigee Evolution', 'Interpreter', 'latex', 'FontSize', ...
    14);
exportgraphics(gca2, 'HW4_6.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Functions

function [A_final, B_final, tf_final, trajectory, error_final] = predictor_corrector(state0, A, B, thrust_acc, mu, a_target, e_target, t0, tf, omega_target)
    % Parameters
    tol = 1e-12;                    % Convergence tolerance
    max_iters = 20;                 % Maximum number of iterations
    delta = 1e-6;                   % Finite difference step size
    alpha = 1.0;                    % Initial line search step size
    min_step_size = 1e-12;          % Minimum allowable step size
    step_size_reduction = 0.5;      % Step size reduction factor
    max_jacobian_cond = 1e9;        % Maximum condition number Jacobian
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

    % Display header for iteration stats
    disp('Starting Predictor-Corrector Iterations');
    disp('-----------------------------------------------------------');
    disp(' Iter |   Norm of Errors   |   A   |   B   |   tf');
    disp('-----------------------------------------------------------');

    for iter = 1:max_iters
        % Simulate trajectory with current A, B, and tf 
        [~, state] = ode113(@(t, y) equations_of_motion(t, y,...
            A, B, thrust_acc, mu), [t0, tf + 1e-14], state0, options);
        r_final = state(end, 1:3)';
        v_final = state(end, 4:6)';

        % Compute errors
        errors = compute_errors(r_final, v_final, A, B, tf, mu,...
            a_target, e_target, omega_target);
        norm_errors = norm(errors);

        % Print iteration stats
        fprintf('%5d | %18.6e | [%5.2f %5.2f %5.2f] | [%5.2f %5.2f %5.2f] | %6.2f\n', ...
                iter, norm_errors, A(1), A(2), A(3), B(1), B(2), B(3), tf);

        % Check convergence
        if norm_errors < tol
            disp('Converged successfully.');
            disp('-----------------------------------------------------------');
            A_final = A;
            B_final = B;
            tf_final = tf;
            trajectory = state;
            error_final = errors;
            return;
        end

        % Compute Jacobian 
        J = compute_jacobian(@equations_of_motion, state0, A, B, tf,...
            delta, t0, thrust_acc, mu, errors, omega_target, a_target, e_target);

        % Check Jacobian condition number
        J_cond = cond(J);
        if J_cond > max_jacobian_cond
            warning('Jacobian condition number too high (%e). Skipping correction.', J_cond);
            continue;
        end

        % Gradient descent with line search
        correction = -J \ errors;  
        step_size = alpha;

        while step_size > min_step_size
            % Test updated parameters
            A_test = A + step_size * correction(1:3);
            B_test = B + step_size * correction(4:6);
            tf_test = tf + step_size * correction(7);  
            [~, state_test] = ode113(@(t, y) equations_of_motion(t,...
                y, A_test, B_test, thrust_acc, mu), [t0 tf_test], state0, options);
            r_test = state_test(end, 1:3)';
            v_test = state_test(end, 4:6)';
            errors_test = compute_errors(r_test, v_test, A_test,...
                B_test, tf_test, mu, a_target, e_target, omega_target);

            % If error decreases, accept step
            if norm(errors_test) < norm_errors
                A = A_test;
                B = B_test;
                tf = tf_test;
                break;
            else
                % Reduce step size
                step_size = step_size * step_size_reduction;  
            end
        end

        % If step size reaches the minimum threshold give warn 
        if step_size <= min_step_size
            warning('Step size reduced to minimum threshold. Guidance may be stuck.');
        end
    end

    % If maximum iterations are reached, return the best solution found
    disp('Failed to converge within the maximum number of iterations.');
    disp('Returning the best solution achieved.');
    disp('-----------------------------------------------------------');
    A_final = A;
    B_final = B;
    tf_final = tf;
    trajectory = state;
    error_final = errors;
end

function J = compute_jacobian(eom, state0, A, B, tf, delta, t0, thrust_acc, mu, errors, omega_target, a_target, e_target)
    % Compute Jacobian using finite differences
    J = zeros(length(errors), 7);  
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

    % Perturb A components
    for i = 1:3
        A_pert = A;
        A_pert(i) = A_pert(i) + delta;
        [~, state_pert] = ode113(@(t, y) eom(t, y, A_pert, B,...
            thrust_acc, mu), [t0 tf], state0, options);
        errors_pert = compute_errors(state_pert(end, 1:3)',...
            state_pert(end, 4:6)', A_pert, B, tf, mu, a_target,...
            e_target, omega_target);
        J(:, i) = (errors_pert - errors) / delta;
    end

    % Perturb B components
    for i = 1:3
        B_pert = B;
        B_pert(i) = B_pert(i) + delta;
        [~, state_pert] = ode113(@(t, y) eom(t, y, A, B_pert,...
            thrust_acc, mu), [t0 tf], state0, options);
        errors_pert = compute_errors(state_pert(end, 1:3)',...
            state_pert(end, 4:6)', A, B_pert, tf, mu, a_target,...
            e_target, omega_target);
        J(:, i + 3) = (errors_pert - errors) / delta;
    end

    % Perturb tf
    tf_pert = tf + delta;
    [~, state_pert] = ode113(@(t, y) eom(t, y, A, B, thrust_acc, mu), ...
        [t0 tf_pert], state0, options);
    errors_pert = compute_errors(state_pert(end, 1:3)', ...
        state_pert(end, 4:6)', A, B, tf_pert, mu, a_target, e_target,...
        omega_target);
    J(:, 7) = (errors_pert - errors) / delta;
end

function dydt = equations_of_motion(t, y, A, B, thrust_acc, mu)
    r = y(1:3);
    v = y(4:6);
    r_norm = norm(r);
    gravity = -mu / r_norm^3 * r;
    
    % LVLH Thrust Components
    u_lvlh = A + B * t; 
    
    % LVLH to Inertial Transformation
    h = cross(r, v);          
    h_hat = h / norm(h);      % Cross-track direction
    r_hat = r / r_norm;       % Radial direction
    theta_hat = cross(h_hat, r_hat);  % Along-track direction
    T_lvlh2inertial = [r_hat, theta_hat, h_hat];  % Transformation matrix

    % Convert to inertial frame
    u_inertial = T_lvlh2inertial * u_lvlh;        
    
    % Compute thrust command
    thrust = thrust_acc * u_inertial / norm(u_inertial); 
    
    % Equations of motion
    dydt = [v; gravity + thrust];
end

function errors = compute_errors(r, v, A, B, tf, mu, a_target, e_target, omega_target)
    % Compute periapsis direction
    rp_dir = [cos(omega_target); sin(omega_target); 0];
    vp_dir = [-sin(omega_target); cos(omega_target); 0]; 

    % Compute target orbital parameters
    p = a_target * (1 - e_target^2); 
    nu = atan2(dot(r, vp_dir), dot(r, rp_dir));  

    % Compute target position and velocity
    rR = p / (1 + e_target * cos(nu));
    rR_vec = rR * (cos(nu) * rp_dir + sin(nu) * vp_dir);
    vR_vec = -sqrt(mu / p) * sin(nu) * rp_dir + sqrt(mu / p) *...
        (1 + e_target * cos(nu)) * vp_dir;
    
    % Compute angular momentum direction
    hR_hat = cross(rR_vec, vR_vec) / norm(cross(rR_vec, vR_vec));
    
    % Compute gravity vector
    g_tf = -mu / norm(r)^3 * r;

    % Lambda vector
    lambda_tf = A + B * tf;

    % Characteristic scales for each error term
    velocity_scale = sqrt(mu / a_target);  
    position_scale = a_target;           
    angular_momentum_scale = 1;           
    gravity_alignment_scale = 1;          
    lambda_magnitude_scale = 1;           
    
    % Compute scaled errors
    errors = [
        (vR_vec(1) - v(1)) / velocity_scale;    
        (vR_vec(2) - v(2)) / velocity_scale;       
        (vR_vec(3) - v(3)) / velocity_scale;       
        (norm(rR_vec) - norm(r)) / position_scale; 
        dot(hR_hat, r) / angular_momentum_scale;   
        (dot(lambda_tf, g_tf) - dot(B, v)) / gravity_alignment_scale; 
        (norm(lambda_tf) - 1) / lambda_magnitude_scale; 
    ];
end

function [r, v] = orbital_elements_to_cartesian(mu, a, e, i, Omega, omega, nu)
    % Orbital elements to Cartesian 
    p = a * (1 - e^2);  % Semi-latus rectum

    % Position and velocity in perifocal coordinates
    r_perifocal = [p * cos(nu) / (1 + e * cos(nu));...
        p * sin(nu) / (1 + e * cos(nu)); 0];
    v_perifocal = [-sqrt(mu / p) * sin(nu);...
        sqrt(mu / p) * (e + cos(nu)); 0];

    % Position and velocity to inertial frame
     R_perifocal_to_inertial = perifocal_to_inertial_matrix(i, Omega, omega);
    r = R_perifocal_to_inertial * r_perifocal;
    v = R_perifocal_to_inertial * v_perifocal;
end

function R = perifocal_to_inertial_matrix(i, Omega, omega)
    % Rotation matrix from perifocal to inertial frame
    R = [ cos(Omega) * cos(omega) - sin(Omega) * sin(omega) * cos(i),...
        -cos(Omega) * sin(omega) - sin(Omega) * cos(omega) * cos(i),...
        sin(Omega) * sin(i);
          sin(Omega) * cos(omega) + cos(Omega) * sin(omega) * cos(i),...
          -sin(Omega) * sin(omega) + cos(Omega) * cos(omega) * cos(i),...
          -cos(Omega) * sin(i);
          sin(omega) * sin(i), cos(omega) * sin(i), cos(i)];
end

function [a, e, i, omega, nu] = cartesian_to_orbital_elements(r, v, mu)
    % Specific angular momentum
    h = cross(r, v);
    h_norm = norm(h);
    
    % Eccentricity vector
    e_vec = (cross(v, h) / mu) - (r / norm(r));
    e = norm(e_vec);
    
    % Semi-major axis
    energy = norm(v)^2 / 2 - mu / norm(r);  
    if energy < 0
        a = -mu / (2 * energy);  % Elliptical orbit
    else
        a = Inf;  % Parabolic or hyperbolic orbit
    end
    
    % Inclination
    i = acos(h(3) / h_norm);  

    % Argument of perigee
    n = cross([0; 0; 1], h);
    n_norm = norm(n);

    % Argument of perigee 
    omega = acos(dot(n, e_vec) / (n_norm * e));
    if e_vec(3) < 0
        omega = 2 * pi - omega;  
    end

    % True anomaly 
    if nargout > 4
        % Radial direction and eccentricity direction
        r_norm = norm(r);
        nu = acos(dot(e_vec, r) / (e * r_norm));  
        if dot(r, v) < 0
            nu = 2 * pi - nu;  
        end
    end
end


function dydt = two_body_eom(~, y, mu)
    % Two-body equations of motion
    r = y(1:3);
    v = y(4:6);
    r_norm = norm(r);
    a = -mu / r_norm^3 * r;  
    dydt = [v; a];
end

