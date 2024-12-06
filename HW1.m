%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:       Giovanni Fereoli
% DATE:         September 6, 2024
% CLASS:        ASEN 6015: Aerospace Vehicle Guidance and Control
% INSTRUCTOR:   Prof. Jay W. Mcmahon
% ASSIGNMENT:   Homework 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Ex. 1
clear; clc; close all;

% Define symbolic variables (Orbital elements and burns)
syms a e i Omega omega M vr vtheta vh mu

% Define additional symbols for anomalies and angle
syms theta E % True anomaly and eccentric anomaly

% Substitute nu with omega + theta (nu = true anomaly)
nu = omega + theta;

% Compute semi-latus rectum p
p = a * (1 - e^2);  % km

% Compute specific angular momentum h
h = sqrt(mu * p);  % km^2/s

% Compute mean motion n
n = sqrt(mu / a^3);  % rad/s

% Compute orbital radius r (using nu = omega + theta)
r = p / (1 + e * cos(theta));  % km

% Define the impulsive burn equations
Delta_a = (2 * a^2 / h) * (e * sin(nu) * vr + (p / r) * vtheta);
Delta_e = (1 / h) * (p * sin(nu) * vr + ((p + r) * cos(nu) + r * e) * vtheta);
Delta_i = (r * cos(theta) / h) * vh;
Delta_Omega = (r * sin(theta) / (h * sin(i))) * vh;
Delta_omega = (1 / (h * e)) * (-p * cos(nu) * vr + (p + r) * sin(nu) * vtheta) - (r * sin(theta) * cos(i) / (h * sin(i))) * vh;
Delta_M = n + (a * sqrt(1 - e^2) / (a * h * e)) * ((p * cos(nu) - 2 * r * e) * vr - (p + r) * sin(nu) * vtheta);

% State vector (orbital elements) - purely symbolic variables
x = [a; e; i; Omega; omega; M];

% Control vector (impulsive burns) - symbolic as well
u = [vr; vtheta; vh];

% Define the state change vector (the Delta terms)
x_dot = [Delta_a; Delta_e; Delta_i; Delta_Omega; Delta_omega; Delta_M];

% Compute the Jacobians
A = jacobian(x_dot, x);  % Jacobian with respect to the state
B = jacobian(x_dot, u);  % Jacobian with respect to the control input

% Compute controllability matrix C = [B AB A^2B A^3B A^4B A^5B]
C = [B A*B A^2*B A^3*B A^4*B A^5*B];

% Display the controllability matrix
disp('Controllability Matrix:');
disp(C);

% Compute and display the rank of the controllability matrix
rank_Controllability = rank(C);
disp('Rank of the Controllability Matrix:');
disp(rank_Controllability);

%% Ex. 2
clear; clc; close;

function dx = closed_loop_dynamics(~, x, xd, K)
    % Get the matrices A and B for the current state x
    % [~, A, B] = GVE(x, zeros(3,1));  
    
    % Get the desired A_d for the desired state xd
    % [~, Ad, ~] = GVE(xd, zeros(3,1)); 

    % Compute control input: u = -(B'B)^-1 B' [A - A_d - K(x - x_d)]
    % u = - inv(B' * B) * B' * (A - Ad + K * (x - xd));

    % Compute control input as HW requires 
    u = - K * (x - xd);
    
    % Update state using GVE with the computed control input
    dx = GVE(x, u);
end

function [x_dot, A, B] = GVE(x, u)
    % Extract the current state (orbital elements)
    a = x(1);  % Semi-major axis
    e = x(2);  % Eccentricity
    i = x(3);  % Inclination
    RAAN = x(4);  % Right Ascension of Ascending Node (RAAN)
    omega = x(5);  % Argument of Perigee
    M = x(6);  % Mean anomaly

    % Gravitational parameter of Earth (mu)
    mu = 3.986e5;  % km^3/s^2 

    % Get true anomaly
    E0 = M;  % Initial guess
    kepler_eq = @(E) E - e*sin(E) - M;
    options = optimoptions('fsolve', 'Display', 'off');  
    E = fsolve(kepler_eq, E0, options);
    nu = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));

    % Utils
    p = a * (1 - e^2);  % Semi-latus rectum
    h = sqrt(mu * p);   % Specific angular momentum
    theta = nu + omega;  % Argument of latitude
    r = p / (1 + e * cos(theta));  % Orbital radius
    n = sqrt(mu / a^3);  % Mean motion
    
    % Intrinsic dynamics 
    A = [0; 0; 0; 0; 0; n];  
    
    % Control influence matrix 
    B = zeros(6, 3);  
    B(1, :) = [2 * a^2 / h * e * sin(nu), 2 * a^2 / h * (p / r), 0];  
    B(2, :) = [p / h * sin(nu), (1 / h) * ((p + r) * cos(nu) + r * e), 0];  
    B(3, :) = [0, 0, (r * cos(theta) / h)]; 
    B(4, :) = [0, 0, (r * sin(theta) / (h * sin(i)))];  
    B(5, :) = [-p * cos(nu) / (h * e), (p + r) * sin(nu) / (h * e),...
        -(r * sin(theta) * cos(i) / (h * sin(i)))];  
    B(6, :) = [(a * sqrt(1 - e^2) / (a * h * e)) * (p * cos(nu) - 2 * r * e),...
        -(a * sqrt(1 - e^2) / (a * h * e)) * (p + r) * sin(nu), 0];  

    % Compute the state derivative (x_dot = A + B*u)
    x_dot = A + B * u;
end

% Simulation parameters
t_end = 6000;    % Final time [s]
dt = 1;       % Time step [s]

% Initial and desired state 
x0 = [7000; 0.1; deg2rad(10); deg2rad(10); deg2rad(10); deg2rad(10)];  
xd = [7001; 0.2; deg2rad(15); deg2rad(15); deg2rad(15); deg2rad(15)];

% Gain matrix 
% K = 0.0001 * diag([1, 1, 1, 1, 1, 1]);  
K = 0.001 * [1, 1, 0, 0, 1, 0; 
               1, 1, 0, 0, 1, 1; 
                  0, 0, 1, 1, 1, 0];

% Simulation
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[time, x_store] = ode45(@(t, x) closed_loop_dynamics(t, x, xd, K),...
    [0 t_end], x0, options);

% Compute deltas for each state variable: delta = x - xd
delta_store = x_store - repmat(xd', length(time), 1);

% Plot delta results 
gca = figure(1);
lineWidth = 1.5;
colors = lines(6); 
subplot(3,2,1);
semilogy(time, abs(delta_store(:,1)), 'Color', colors(1,:), 'LineWidth', lineWidth);
xlabel('Time [s]', 'FontSize', 12);
ylabel('|\Delta a| [km]', 'FontSize', 12, 'Interpreter', 'tex');
grid on;
subplot(3,2,2);
semilogy(time, abs(delta_store(:,2)), 'Color', colors(2,:), 'LineWidth', lineWidth);
xlabel('Time [s]', 'FontSize', 12);
ylabel('|\Delta e| [-]', 'FontSize', 12, 'Interpreter', 'tex');
grid on;
subplot(3,2,3);
semilogy(time, abs(delta_store(:,3)), 'Color', colors(3,:), 'LineWidth', lineWidth);
xlabel('Time [s]', 'FontSize', 12);
ylabel('|\Delta i| [rad]', 'FontSize', 12, 'Interpreter', 'tex');
grid on;
subplot(3,2,4);
semilogy(time, abs(delta_store(:,4)), 'Color', colors(4,:), 'LineWidth', lineWidth);
xlabel('Time [s]', 'FontSize', 12);
ylabel('|\Delta \Omega| [rad]', 'FontSize', 12, 'Interpreter', 'tex');
grid on;
subplot(3,2,5);
semilogy(time, abs(delta_store(:,5)), 'Color', colors(5,:), 'LineWidth', lineWidth);
xlabel('Time [s]', 'FontSize', 12);
ylabel('|\Delta \omega| [rad]', 'FontSize', 12, 'Interpreter', 'tex');
grid on;
subplot(3,2,6);
semilogy(time, abs(delta_store(:,6)), 'Color', colors(6,:), 'LineWidth', lineWidth);
xlabel('Time [s]', 'FontSize', 12);
ylabel('|\Delta M| [rad]', 'FontSize', 12, 'Interpreter', 'tex');
grid on;
set(gcf, 'Position', [100, 100, 800, 600]);  
exportgraphics(gca, 'HW1_a.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Ex. 3a
clear; close; clc;

% Constants
mu = 3.986e5;   % Gravitational parameter of Earth [km^3/s^2]
T = 86400;      % Orbital period GEO (24 hours in seconds)
Re = 6378;      % Radius of the Earth [km]
r_initial = Re + 200;       % Radius of parking orbit [km]
r_geo = (mu * T^2 / (4 * pi^2))^(1/3);  % GEO orbit radius
r_desired = r_geo;  

% Misc
T_initial = 2 * pi * sqrt(r_initial^3 / mu);
T_desired = 2 * pi * sqrt(r_desired^3 / mu);
v_initial = sqrt(mu / r_initial);  
v_desired = sqrt(mu / r_desired);  

% Velocities for the transfer orbit
a_t = (r_initial + r_desired) / 2;  % Semi-major axis of transfer orbit
e_t = (r_desired - r_initial) / (r_initial + r_desired);  % Eccentricity of transfer orbit
T_t = pi * sqrt(a_t^3 / mu);
v_transfer_perigee = sqrt(mu * (1 + e_t) / ((1 - e_t) * a_t));  
v_transfer_apogee = sqrt(mu * (1 - e_t) / ((1 + e_t) * a_t)); 

% Delta-v's
Delta_v1 = v_transfer_perigee - v_initial; 
Delta_v2 = v_desired - v_transfer_apogee; 
fprintf('Delta-v1: %.3f km/s\n', Delta_v1);
fprintf('Delta-v2: %.3f km/s\n', Delta_v2);

% Two-body dynamics function
function dx = two_body_dynamics(~, x, mu)
    r = x(1:3); 
    v = x(4:6);
    r_norm = norm(r);
    a = - mu / r_norm^3 * r;
    dx = [v; a];
end

% Simulate Initial Parking Orbit (before the transfer)
r0_parking = [r_initial; 0; 0];  
v0_parking = [0; v_initial; 0];  
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[~, X_parking] = ode45(@(t,x) two_body_dynamics(t, x, mu), [0, T_initial],...
    [r0_parking; v0_parking], options);

% Simulate Hohmann Transfer (after first impulse)
r0_transfer = r0_parking;  % Final position from parking orbit
v1_hat = v0_parking / norm(v0_parking);
v0_transfer = v0_parking + Delta_v1 * v1_hat;  % Apply first impulse
[~, X_transfer] = ode45(@(t,x) two_body_dynamics(t, x, mu), [0, T_t],...
    [r0_transfer; v0_transfer], options);

% Simulate Hohmann Transfer (after second impulse)
r0_final = X_transfer(end, 1:3)';  % Final position from transfer orbit
v2_hat = - [0; v_desired; 0] / norm([0; v_desired; 0]);
v0_final = X_transfer(end, 4:6)' + Delta_v2 * v2_hat;  % Apply second impulse
[~, X_final] = ode45(@(t,x) two_body_dynamics(t, x, mu), [0, T_desired],...
    [r0_final; v0_final], options);

% Simulate Final GEO Orbit (after second impulse)
r0_desired = [r_desired; 0; 0];  
v0_desired = [0; v_desired; 0];
[~, X_desired] = ode45(@(t,x) two_body_dynamics(t, x, mu), [0, T_desired],...
    [r0_desired; v0_desired], options);

% Plot the results
gca = figure(2);
hold on;
plot3(X_parking(:,1), X_parking(:,2), X_parking(:,3), 'k--', 'LineWidth', 3);  
plot3(X_transfer(:,1), X_transfer(:,2), X_transfer(:,3), 'r', 'LineWidth', 2);  
plot3(X_final(:,1), X_final(:,2), X_final(:,3), 'b', 'LineWidth', 2);  
plot3(X_desired(:,1), X_desired(:,2), X_desired(:,3), 'k--', 'LineWidth', 3);  
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
plot_earth();
legend('Parking Orbit', 'Transfer Orbit', 'Final Orbit', 'Desired GEO Orbit');
grid on;
axis equal;
hold off;
exportgraphics(gca, 'HW1_b.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Ex. 3b 
clear; close; clc;

% Constants
mu = 3.986e5; 
T = 86400;      
Re = 6378;     
r_initial = Re + 200;      
r_geo = (mu * T^2 / (4 * pi^2))^(1/3); 
r_desired = r_geo;  

% Misc
T_initial = 2 * pi * sqrt(r_initial^3 / mu);
T_desired = 2 * pi * sqrt(r_desired^3 / mu);
v_initial = sqrt(mu / r_initial); 
v_desired = sqrt(mu / r_desired);  

% Velocities for the transfer orbit
a_t = (r_initial + r_desired) / 2;  % Semi-major axis of transfer orbit
e_t = (r_desired - r_initial) / (r_initial + r_desired);  % Eccentricity of transfer orbit
T_t = pi * sqrt(a_t^3 / mu);
v_transfer_perigee = sqrt(mu * (1 + e_t) / ((1 - e_t) * a_t));  
v_transfer_apogee = sqrt(mu * (1 - e_t) / ((1 + e_t) * a_t)); 

% Delta-v's
Delta_v1 = v_transfer_perigee - v_initial;  % First impulse
Delta_v2 = v_desired - v_transfer_apogee;  % Second impulse
fprintf('Delta-v1: %.3f km/s\n', Delta_v1);
fprintf('Delta-v2: %.3f km/s\n', Delta_v2);

% Simulate Initial Parking Orbit (before the transfer)
r0_parking = [r_initial; 0; 0];  
v0_parking = [0; v_initial; 0];  
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[~, X_parking] = ode45(@(t,x) two_body_dynamics(t, x, mu),...
    [0, T_initial], [r0_parking; v0_parking], options);

% Simulate Hohmann Transfer (after first impulse)
angle_error = deg2rad(1);  % 1 degree error
rotation_matrix = [cos(angle_error), -sin(angle_error), 0;...
    sin(angle_error), cos(angle_error), 0; 0, 0, 1];
r0_transfer = r0_parking;  % Final position from parking orbit
v1_hat = v0_parking / norm(v0_parking);
v0_transfer = v0_parking + rotation_matrix * Delta_v1 * v1_hat;  % Apply first impulse
[~, X_transfer] = ode45(@(t,x) two_body_dynamics(t, x, mu), ...
    [0, T_t], [r0_transfer; v0_transfer], options);

% Simulate Final Orbit (after second impulse)
r0_final = X_transfer(end, 1:3)';  % Final position from transfer orbit
v2_hat = - [0; v_desired; 0] / norm([0; v_desired; 0]);
v0_final = X_transfer(end, 4:6)' + Delta_v2 * v2_hat;  % Apply second impulse
[~, X_final] = ode45(@(t,x) two_body_dynamics(t, x, mu),...
    [0, T_desired], [r0_final; v0_final], options);

% Nominal final orbit (without pointing error)
r0_desired = [r_desired; 0; 0];  
v0_desired = [0; v_desired; 0];
[~, X_desired] = ode45(@(t,x) two_body_dynamics(t, x, mu),...
    [0, T_desired], [r0_desired; v0_desired], options);

% Compute Orbital Elements
function [a, e] = rv2oe(r, v, mu)
    % Compute semi-major axis and eccentricity from position and velocity vectors
    r_norm = norm(r);
    v_norm = norm(v);
    
    % Specific energy
    energy = v_norm^2 / 2 - mu / r_norm;
    
    % Semi-major axis
    a = -mu / (2 * energy);
    
    % Eccentricity vector
    e_vec = (1 / mu) * ((v_norm^2 - mu / r_norm) * r - dot(r, v) * v);
    e = norm(e_vec);
end

% Nominal Transfer Orbit
v1_hat = v0_parking / norm(v0_parking);
v0_transfer_ideal = v0_parking + Delta_v1 * v1_hat;  % Apply first impulse
[~, X_transfer_ideal] = ode45(@(t,x) two_body_dynamics(t, x, mu), ...
    [0, T_t], [r0_transfer; v0_transfer_ideal], options);
[a_nominal, e_nominal] = rv2oe(X_transfer_ideal(1, 1:3)', X_transfer_ideal(1, 4:6)', mu);

% Transfer Orbit with Pointing Error
[a_transfer, e_transfer] = rv2oe(X_transfer(end, 1:3)', X_transfer(end, 4:6)', mu);

% Final Orbit with Pointing Error
[a_final, e_final] = rv2oe(X_final(end, 1:3)', X_final(end, 4:6)', mu);

% Desired Final Orbit (Nominal)
[a_desired, e_desired] = rv2oe(X_desired(end, 1:3)', X_desired(end, 4:6)', mu);

% Display Orbital Elements and Errors
fprintf('\nNominal Transfer Orbit: a = %.10f km, e = %.10f\n', a_nominal, e_nominal);
fprintf('Transfer Orbit (With Error): a = %.10f km, e = %.10f\n', a_transfer, e_transfer);
fprintf('Final Orbit (With Error): a = %.10f km, e = %.10f\n', a_final, e_final);
fprintf('Desired Final Orbit: a = %.10f km, e = %.10f\n', a_desired, e_desired);

% Errors in transfer orbit and final orbit
fprintf('\nErrors in Transfer Orbit:\n');
fprintf('Semi-major axis error (transfer): %.10f km\n', abs(a_transfer - a_nominal));
fprintf('Eccentricity error (transfer): %.10f\n', abs(e_transfer - e_nominal));

fprintf('\nErrors in Final Orbit:\n');
fprintf('Semi-major axis error (final): %.10f km\n', abs(a_final - a_desired));
fprintf('Eccentricity error (final): %.10f\n', abs(e_final - e_desired));

% Plot the results
gca = figure(3);
hold on;
plot3(X_parking(:,1), X_parking(:,2), X_parking(:,3), 'k--', 'LineWidth', 3); 
plot3(X_transfer(:,1), X_transfer(:,2), X_transfer(:,3), 'r', 'LineWidth', 2);  
plot3(X_final(:,1), X_final(:,2), X_final(:,3), 'b', 'LineWidth', 2);  
plot3(X_desired(:,1), X_desired(:,2), X_desired(:,3), 'k--', 'LineWidth', 3);  
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
plot_earth();
legend('Parking Orbit', 'Transfer Orbit', 'Final Orbit', 'Desired GEO Orbit');
grid on;
axis equal;
hold off;
exportgraphics(gca, 'HW1_c.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Ex. 3c 
clear; close; clc;

% Constants
mu = 3.986e5;   
T = 86400;      
Re = 6378;      
r_initial = Re + 200;       
r_geo = (mu * T^2 / (4 * pi^2))^(1/3);  
r_desired = r_geo; 

% Compute Hohmann Transfer and Delta-V's
% T_initial = 2 * pi * sqrt(r_initial^3 / mu);
% T_desired = 2 * pi * sqrt(r_desired^3 / mu);
v_initial = sqrt(mu / r_initial); 
v_desired = sqrt(mu / r_desired);  

% Velocities for the transfer orbit
a_t = (r_initial + r_desired) / 2;  % Semi-major axis of transfer orbit
e_t = (r_desired - r_initial) / (r_initial + r_desired);  % Eccentricity of transfer orbit
% T_t = pi * sqrt(a_t^3 / mu);
v_transfer_perigee = sqrt(mu * (1 + e_t) / ((1 - e_t) * a_t));  
v_transfer_apogee = sqrt(mu * (1 - e_t) / ((1 + e_t) * a_t));  

% Delta-v's
Delta_v1 = v_transfer_perigee - v_initial;  % First impulse
Delta_v2 = v_desired - v_transfer_apogee;  % Second impulse
fprintf('Delta-v1: %.3f km/s\n', Delta_v1);
fprintf('Delta-v2: %.3f km/s\n', Delta_v2);

% Define pointing angle perturbations (in radians)
perturb_angle = deg2rad(1);  

% Simulate the nominal trajectory (no pointing error)
[~, a_nominal, e_nominal] = simulate_hohmann(mu, r_initial, r_desired,...
    Delta_v1, Delta_v2, 0);

% Simulate with small positive perturbation
[~, a_pos_perturb, e_pos_perturb] = simulate_hohmann(mu, r_initial, r_desired,...
    Delta_v1, Delta_v2, perturb_angle);

% Simulate with small negative perturbation
[~, a_neg_perturb, e_neg_perturb] = simulate_hohmann(mu, r_initial, r_desired,...
    Delta_v1, Delta_v2, -perturb_angle);

% Compute Sensitivities (finite difference method)
sensitivity_a = (a_pos_perturb - a_nominal) / (perturb_angle);  % da/dtheta
sensitivity_e = (e_pos_perturb - e_nominal) / (perturb_angle);  % de/dtheta

fprintf('Sensitivity of semi-major axis: %.10f km/rad\n', sensitivity_a);
fprintf('Sensitivity of eccentricity: %.10f /rad\n', sensitivity_e);

% Calculate Expected Errors for 0.1deg and -0.3deg pointing errors
error_1 = deg2rad(0.1);  % 0.1 degree
error_2 = deg2rad(-0.3);  % -0.3 degree

% Expected errors in the final state
expected_error_a_1 = sensitivity_a * error_1;  % Da for 0.1deg
expected_error_e_1 = sensitivity_e * error_1;  % Da for 0.1deg
expected_error_a_2 = sensitivity_a * error_2;  % De for -0.3deg
expected_error_e_2 = sensitivity_e * error_2;  % De for -0.3deg

fprintf('\nExpected errors for 0.1° pointing error:\n');
fprintf('∆a = %.10f km\n', expected_error_a_1);
fprintf('∆e = %.10f\n', expected_error_e_1);

fprintf('\nExpected errors for -0.3° pointing error:\n');
fprintf('∆a = %.10f km\n', expected_error_a_2);
fprintf('∆e = %.10f\n', expected_error_e_2);

% Compare Expected Errors with Simulation Results
% Simulate with 0.1° pointing error
[~, a_sim_1, e_sim_1] = simulate_hohmann(mu, r_initial, r_desired, ...
    Delta_v1, Delta_v2, error_1);

% Simulate with -0.3° pointing error
[~, a_sim_2, e_sim_2] = simulate_hohmann(mu, r_initial, r_desired, ...
    Delta_v1, Delta_v2, error_2);

% Compute simulation errors
sim_error_a_1 = a_sim_1 - a_nominal;  % da from simulation for 0.1°
sim_error_e_1 = e_sim_1 - e_nominal;  % de from simulation for 0.1°
sim_error_a_2 = a_sim_2 - a_nominal;  % da from simulation for -0.3°
sim_error_e_2 = e_sim_2 - e_nominal;  % de from simulation for -0.3°

fprintf('\nSimulation errors for 0.1° pointing error:\n');
fprintf('∆a = %.10f km\n', sim_error_a_1);
fprintf('∆e = %.10f\n', sim_error_e_1);

fprintf('\nSimulation errors for -0.3° pointing error:\n');
fprintf('∆a = %.10f km\n', sim_error_a_2);
fprintf('∆e = %.10f\n', sim_error_e_2);

% Hohmann transfer function (with pointing error)
function [X_final, a_final, e_final] = simulate_hohmann(mu, r_initial,...
    r_desired, Delta_v1, Delta_v2, angle_error)
    % Initial parking orbit
    v_initial = sqrt(mu / r_initial);
    r0_parking = [r_initial; 0; 0];  
    v0_parking = [0; v_initial; 0];  

    % Apply first impulse with pointing error
    v1_hat = [0; 1; 0];
    rotation_matrix = [cos(angle_error), -sin(angle_error), 0; ...
        sin(angle_error), cos(angle_error), 0; 0, 0, 1];
    v0_transfer = v0_parking + rotation_matrix * Delta_v1 * v1_hat;  

    % Simulate transfer orbit
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    [~, X_transfer] = ode45(@(t,x) two_body_dynamics(t, x, mu), ...
        [0, pi * sqrt(((r_initial + r_desired) / 2)^3 / mu)],...
        [r0_parking; v0_transfer], options);

    % Apply second impulse (no error here)
    r0_final = X_transfer(end, 1:3)';  % Final position from transfer orbit
    v2_hat = [0; -1; 0];
    v0_final = X_transfer(end, 4:6)' + Delta_v2 * v2_hat;  % Apply second impulse

    % Simulate final orbit
    [~, X_final] = ode45(@(t,x) two_body_dynamics(t, x, mu), ...
        [0, 3600 * 24], [r0_final; v0_final], options);

    % Compute orbital elements of the final state
    [a_final, e_final] = rv2oe(X_final(end, 1:3)', X_final(end, 4:6)', mu);
end

%% Ex. 3d
clear; close; clc;

% Constants
mu = 3.986e5;   
T = 86400;      
Re = 6378;     
a_acceleration = 30e-3;  
r_initial = Re + 200;       
r_geo = (mu * T^2 / (4 * pi^2))^(1/3);  
r_desired = r_geo; 

% Compute Hohmann Transfer and Delta-V's
T_initial = 2 * pi * sqrt(r_initial^3 / mu);
T_desired = 2 * pi * sqrt(r_desired^3 / mu);
v_initial = sqrt(mu / r_initial);  
v_desired = sqrt(mu / r_desired);  

% Velocities for the transfer orbit
a_t = (r_initial + r_desired) / 2;  % Semi-major axis of transfer orbit
e_t = (r_desired - r_initial) / (r_initial + r_desired);  % Eccentricity of transfer orbit
T_t = pi * sqrt(a_t^3 / mu);
v_transfer_perigee = sqrt(mu * (1 + e_t) / ((1 - e_t) * a_t)); 
v_transfer_apogee = sqrt(mu * (1 - e_t) / ((1 + e_t) * a_t));  

% Delta-v's
Delta_v1 = v_transfer_perigee - v_initial;  % First impulse
Delta_v2 = v_desired - v_transfer_apogee;  % Second impulse
fprintf('Delta-v1: %.3f km/s\n', Delta_v1);
fprintf('Delta-v2: %.3f km/s\n', Delta_v2);

% Burn Time Calculation (finite burn model)
t_burn1 = Delta_v1 / a_acceleration;  % Burn time for first burn [s]
t_burn2 = Delta_v2 / a_acceleration;  % Burn time for second burn [s]
fprintf('\nBurn time for first burn: %.6f s\n', t_burn1);
fprintf('Burn time for second burn: %.6f s\n', t_burn2);

% Simulate the finite burn for the Hohmann transfer
% Initial parking orbit
v_initial = sqrt(mu / r_initial);  % Initial velocity in parking orbit
r0_parking = [r_initial; 0; 0];  
v0_parking = [0; v_initial; 0];  

% Apply first burn over the duration of t_burn1
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[t1, X1] = ode45(@(t,x) finite_burn_dynamics(t, x, mu, a_acceleration,...
    t_burn1), [0, t_burn1], [r0_parking; v0_parking], options);

% Transfer orbit after first burn
r_transfer = X1(end, 1:3)';  % Final position after first burn
v_transfer = X1(end, 4:6)';  % Final velocity after first burn

% Simulate coast phase in transfer orbit (no thrust)
T_transfer = pi * sqrt(((r_initial + r_geo) / 2)^3 / mu); 
[~, X_transfer] = ode45(@(t,x) two_body_dynamics(t, x, mu),...
    [0, T_transfer], [r_transfer; v_transfer], options);

% Apply second burn over the duration of t_burn2 at apogee
r_final = X_transfer(end, 1:3)';  % Position at apogee
v_apogee = X_transfer(end, 4:6)';  % Velocity at apogee

% Simulate second burn
[t2, X2] = ode45(@(t,x) finite_burn_dynamics(t, x, mu, a_acceleration,...
    t_burn2), [0, t_burn2], [r_final; v_apogee], options);
r_final = X2(end, 1:3)';
v_final = X2(end, 4:6)';
[~, X_final] = ode45(@(t,x) two_body_dynamics(t, x, mu),...
    [0, T_desired], [r_final; v_final], options);

% Compare Final Orbit with Desired GEO Orbit
r_desired = [r_geo; 0; 0];
v_desired = [0; sqrt(mu / r_geo); 0];

% Final Orbit with Pointing Error
[a_final, e_final] = rv2oe(r_final', v_final', mu);

% Desired Final Orbit (Nominal)
[a_desired, e_desired] = rv2oe(r_desired', v_desired', mu);

% Display Orbital Elements and Errors
fprintf('\nFinite burn vs desired GEO:\n');
fprintf('Final Orbit (With Finite Burn): a = %.10f km, e = %.10f\n', a_final, e_final);
fprintf('Desired Final Orbit: a = %.10f km, e = %.10f\n', a_desired, e_desired);

% Errors in transfer orbit and final orbit
fprintf('\nErrors in Final Orbit:\n');
fprintf('Semi-major axis error (final): %.10f km\n', abs(a_final - a_desired));
fprintf('Eccentricity error (final): %.10f\n', abs(e_final - e_desired));


% Plots
gca = figure(4);
hold on;
theta = linspace(0, 2*pi, 100);
parking_orbit_x = r_initial * cos(theta);
parking_orbit_y = r_initial * sin(theta);
plot(parking_orbit_x, parking_orbit_y, 'k--', 'LineWidth', 2);  
plot(X_transfer(:,1), X_transfer(:,2), 'r', 'LineWidth', 1.5);  
geo_orbit_x = r_geo * cos(theta);
geo_orbit_y = r_geo * sin(theta);
plot(geo_orbit_x, geo_orbit_y, 'k--', 'LineWidth', 2);  
plot(X_final(:,1), X_final(:,2), 'b', 'LineWidth', 1.5); 
plot_earth();

% Plot the finite burns as quiver arrows during the burn times
plot_burn_quivers(t1, X1, a_acceleration, 'First Burn', 'm');
plot_burn_quivers(t2, X2, a_acceleration, 'Second Burn', 'm');

% Labels and formatting
xlabel('x [km]');
ylabel('y [km]');
legend('Parking Orbit', 'Transfer Orbit', 'Desired GEO Orbit',...
    'Final Reached Orbit');
grid on;
axis equal;
hold off;
exportgraphics(gca, 'HW1_d.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Function to plot quiver arrows during finite burns
function plot_burn_quivers(t, X, a_acceleration, ~, color)
    % Determine time intervals for quiver plotting
    burn_time_interval = linspace(min(t), max(t), 10);  % 10 arrows along the burn
    for i = 1:length(burn_time_interval)
        % Find the corresponding state at the burn time
        idx = find(t >= burn_time_interval(i), 1);
        if isempty(idx), continue; end
        
        r_burn = X(idx, 1:3);  % Position
        v_burn = X(idx, 4:6);  % Velocity
        norm_v_burn = norm(v_burn);  % Magnitude of velocity
        
        % Plot the quiver arrow representing acceleration
        quiver(r_burn(1), r_burn(2), ...
            200000 * a_acceleration * v_burn(1) / norm_v_burn, ...
            200000 * a_acceleration * v_burn(2) / norm_v_burn, ...
               0, 'MaxHeadSize', 1.5, 'LineWidth', 1.5, 'Color', color);
    end
end

% Finite burn dynamics (constant thrust model)
function dx = finite_burn_dynamics(t, x, mu, a_acceleration, t_burn)
    r = x(1:3); v = x(4:6);
    r_norm = norm(r);
    
    % Gravity acceleration
    a_gravity = -mu / r_norm^3 * r;
    
    % Constant thrust acceleration (during burn)
    if t <= t_burn
        a_thrust = a_acceleration * v / norm(v);  % Thrust in the direction of velocity
    else
        a_thrust = [0; 0; 0];  % No thrust after burn
    end
    
    % Total acceleration
    a_total = a_gravity + a_thrust;
    
    % State derivatives
    dx = [v; a_total];
end

%% Utils

function plot_earth()
    % Earth's radius in kilometers
    earth_radius = 6378;  % km

    % Create a sphere to represent the Earth
    [x, y, z] = sphere(50);          % Create a sphere with a 50x50 grid
    
    % Scale the sphere to Earth's radius
    x = x * earth_radius;
    y = y * earth_radius;
    z = z * earth_radius;

    % Load the topographic data provided by MATLAB (you can use your own Earth texture if needed)
    load topo;                       % Load MATLAB's topographic data

    % Plot the spherical surface
    s = surface(x, y, z);
    
    % Apply texture mapping and set topographic data as the color data
    s.FaceColor = 'texturemap';      % Use texture mapping
    s.CData = topo;                  % Set color data to topographic data
    s.EdgeColor = 'none';            % Remove edges
    s.FaceLighting = 'gouraud';      % Set Gouraud lighting for smooth curved surface lighting
    s.SpecularStrength = 0.4;        % Adjust strength of the specular reflection

    % Add a light source
    light('Position', [-1, 0, 1]);   % Position the light
end

