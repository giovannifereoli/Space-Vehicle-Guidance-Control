%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:       Giovanni Fereoli
% DATE:         October 19, 2024
% CLASS:        ASEN 6015: Aerospace Vehicle Guidance and Control
% INSTRUCTOR:   Prof. Jay W. Mcmahon
% ASSIGNMENT:   Homework 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Ex. 1
clear; clc; close all;

% Define system matrices
A = [5  0  3  0  1;
     3  0  0 -2  0;
     0 -2  4  1  0;
     1  3 -4  1  3;
     0  2  2  0 -1];

B = [0 1;
     0 2;
     0 0;
     1 3;
     1 1];

% Define initial state
x0 = [1; -1; 2; 1; -0.5];  

% Time span for simulation
tspan = [0 15]; 

% Define cases and parameters
cases = struct('gamma', {1, 1, 1}, ...
               'a', {1, 0, 1}, ...
               'b', {1, 0, 0}, ...
               'label', {'Case (a)', ...
                         'Case (b)', ...
                         'Case (c)'});

% Preallocate storage for results
results = struct('t', [], 'x', [], 'cost', [], 'label', []);

% LQR simulation and cost function computation 
function [t, x, J] = simulate_lqr_and_cost_trapz(A, B, Q, R, x0, tspan)
    try
        % Solve the Riccati equation and compute LQR gain
        [K, ~, ~] = lqr(A, B, Q, R);

        % Define closed-loop dynamics
        dynamics = @(t, x) (A - B * K) * x;

        % Integrate dynamics 
        options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
        [t, x] = ode45(dynamics, tspan, x0, options);

        % Compute control input u(t) and instantaneous cost
        cost = zeros(length(t), 1);  % Initialize cost array
        for i = 1:length(t)
            u = -K * x(i, :)';  % Compute control input
            cost(i) = x(i, :) * Q * x(i, :)' + u' * R * u;  % Instantaneous cost
        end

        % Compute total cost using trapz for numerical integration
        J = trapz(t, cost);

    catch ME
        % If an error occurs, display the error message and return NaNs
        disp('Error using lqr:');
        disp(ME.message);
        t = NaN;
        x = NaN;
        J = NaN;
    end
end

% Loop through each case
for i = 1:length(cases)
    % Print the current case number
    fprintf('Case %d\n:', i);

    % Define cost function weights for each case
    Q = diag([0 0 cases(i).gamma 0 0]);
    R = [cases(i).a 0; 0 cases(i).b];

    % Simulate LQR and compute the cost function using trapz
    [results(i).t, results(i).x, results(i).cost] = simulate_lqr_and_cost_trapz(A,...
        B, Q, R, x0, tspan);
    results(i).label = cases(i).label;

    % Display the computed cost function
    fprintf('Cost for %s: %.10f\n', cases(i).label, results(i).cost);
    fprintf('\n');
end

% Plots
for i = 1:length(cases)
    if ~isnan(results(i).t)
        figure(i); 
        plot(results(i).t, results(i).x, 'LineWidth', 1.5);
        legend({'$x_1$', '$x_2$', '$x_3$', '$x_4$', '$x_5$'},...
            'Interpreter', 'latex', 'Location', 'northeast');
        grid on;
        xlabel('Time [s]', 'Interpreter', 'latex');
        ylabel('States [-]', 'Interpreter', 'latex');
    end
end
exportgraphics(gca, 'HW3_1a.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot control action, Case 1
[K_1, S_1, E_1] = lqr(A, B, diag([0 0 cases(1).gamma 0 0]),...
    [cases(1).a 0; 0 cases(1).b]);
figure(i+1); 
plot(results(1).t, - K_1 * (results(1).x)', 'LineWidth', 1.5);
legend({'$u_1$', '$u_2$'},...
    'Interpreter', 'latex', 'Location', 'northeast');
grid on;
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Control [-]', 'Interpreter', 'latex');
exportgraphics(gca, 'HW3_1b.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Ex. 2a
clear; clc; close all;

% Define system matrices (from Problem 1.a)
A = [5  0  3  0  1;
     3  0  0 -2  0;
     0 -2  4  1  0;
     1  3 -4  1  3;
     0  2  2  0 -1];

B = [0 1;
     0 2;
     0 0;
     1 3;
     1 1];

C = [0 0 1 0 0];

% Compute matrix K (from Problem 1.a)
a = 1;
b = 1;
gamma = 1;
Q = diag([0 0 gamma 0 0]);
R = [a 0; 0 b];
[K_a, ~, ~] = lqr(A, B, Q, R);

% Define matrix H = A - B * K
H = A - B * K_a;

% Define L
L = eye(5);

% Power spectral density for white noise
R_w = diag([0.1, 0.2, 0.1, 0.3, 0.2]);

% Define time span 
t0 = 0;
tf = 3;  

% Define the adjoint system 
adjoint_system = @(t, p) -(H') * p;

% Set final condition for p(tf)
p_final = C';

% Solve adjoint system backward in time
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[t_adj, p_sol] = ode45(adjoint_system, [tf t0], p_final, options);

% Plot the adjoint vector p over time
figure;
plot(t_adj, p_sol,'LineWidth', 1.5);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Adjoints [-]', 'Interpreter', 'latex');
grid on;
legend({'$p_1$', '$p_2$', '$p_3$', '$p_4$', '$p_5$'},...
            'Interpreter', 'latex', 'Location', 'northwest');
exportgraphics(gca, 'HW3_2a.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Function to compute sigma^2 for each time step
function sigma_squared_dtau = compute_sigma_squared_dtau(t_adj, p_sol, L, R_w)
    % Initialize 
    sigma_squared_dtau = zeros(length(t_adj), 5); 
    % Loop through each time step
    for i = 1:length(t_adj)
        % Get adjoint vector at current time
        p_t = p_sol(i, :)'; 
        % Compute sigma^2 for each channel
        for k = 1:5
            sigma_squared_dtau(i, k) = (p_t(k)^2) * L(k,k) * R_w(k,k) * L(k,k);
        end
    end
end

% Compute sigma_squared_dtau for each time step
sigma_squared_dtau_sol = compute_sigma_squared_dtau(t_adj, p_sol, L, R_w);

% Use trapezoidal rule to integrate for each channel
integral_sigma_squared = zeros(1, 5);
for k = 1:5
    integral_sigma_squared(k) = -trapz(t_adj, sigma_squared_dtau_sol(:, k));  
end

% Compute the total effect 
total_effect = sum(integral_sigma_squared);

% Display the final result for each channel and the total effect
disp('Final integral values of sigma^2 for each channel:');
disp(integral_sigma_squared);
disp('Final total effect:');
disp(total_effect);

% Create the figure 
figure;
hold on;
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; ...
          0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330];
for i = 1:5
    bar(i, integral_sigma_squared(i), 'FaceColor', colors(i, :));
end
bar(6, total_effect, 'FaceColor', colors(6, :));
ylabel('$\sigma^2$ [-]', 'Interpreter', 'latex');
xticks(1:6);  
xticklabels({'Channel 1', 'Channel 2', 'Channel 3', ...
    'Channel 4', 'Channel 5', 'Total'}); 
set(gca, 'TickLabelInterpreter', 'latex');
grid on;
hold off;
exportgraphics(gca, 'HW3_2b.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Ex. 2b
clear; close all; clc;

% Define system matrices
A = [5  0  3  0  1;
     3  0  0 -2  0;
     0 -2  4  1  0;
     1  3 -4  1  3;
     0  2  2  0 -1];

B = [0 1;
     0 2;
     0 0;
     1 3;
     1 1];

% Define Q and R matrices for LQR
a = 1;
b = 1;
gamma = 1;
Q = diag([0 0 gamma 0 0]);
R = [a 0; 0 b];

% Define a initial state
x0 = [1; -1; 2; 1; -0.5];  

% Time span for simulation
tspan = [0 15]; 

% Final condition for the Riccati equation
Pf = zeros(5); 

% Solve the Riccati equation backward in time 
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[t_backward, P_sol] = ode45(@(t, P_vec) riccati_ode(t, P_vec, A, B, Q, R),...
    tspan(end:-1:1), Pf(:), options);

% Reverse the time vector 
t_backward = flipud(t_backward);
P_sol = flipud(P_sol);

% Reshape the Riccati solution into matrices
P_sol = reshape(P_sol.', size(A,1), size(A,1), []);

% Preallocate memory for control gain matrix K(t)
K_sol = zeros(length(t_backward), size(B, 2), size(A, 1));

% Compute the time-varying gain matrix K(t) at each time step
for i = 1:length(t_backward)
    K_sol(i, :, :) = R \ (B' * P_sol(:, :, i));   
end

% Integrate state equations forward in time 
[t_forward, x_sol] = ode45(@(t, x) state_dynamics(t, x, t_backward, K_sol, A, B),...
    tspan, x0, options);

% Compute the control input and the cost function J
cost = zeros(length(t_forward), 1);  
for i = 1:length(t_forward)
    % Interpolate the gain matrix K(t) for the current time
    K_t = interp1(t_backward, K_sol, t_forward(i), 'linear');
    
    % Compute the control input
    u = -squeeze(K_t) * x_sol(i, :)';
    
    % Compute the instantaneous cost
    cost(i) = x_sol(i, :) * Q * x_sol(i, :)' + u' * R * u;
end

% Compute the total cost J using trapz to integrate over time
J = trapz(t_forward, cost);

% Print the computed cost function J
fprintf('The computed cost function J is: %.10f\n', J);

% Plot the state variables over time
figure;
plot(t_forward, x_sol, 'LineWidth', 1.5);
legend({'$x_1$', '$x_2$', '$x_3$', '$x_4$', '$x_5$'},...
    'Interpreter', 'latex', 'Location', 'best');
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('States [-]', 'Interpreter', 'latex');
grid on;
exportgraphics(gca, 'HW3_3a.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot all gains in K_sol
figure;
plot(t_backward, squeeze(K_sol(:, 1, :)), 'LineWidth', 1.5);
hold on;
plot(t_backward, squeeze(K_sol(:, 2, :)), 'LineWidth', 1.5);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Gain values [-]', 'Interpreter', 'latex');
legend({'$K_{1,1}$', '$K_{1,2}$', '$K_{1,3}$', '$K_{1,4}$', '$K_{1,5}$',...
    '$K_{2,1}$', '$K_{2,2}$', '$K_{2,3}$', '$K_{2,4}$', '$K_{2,5}$'}, ...
    'Interpreter', 'latex', 'Location', 'northeast');
grid on;
exportgraphics(gca, 'HW3_3b.pdf', 'ContentType','image',...
    'Resolution', 1000);


% Riccati ODE function (Backward integration)
function dP_vec = riccati_ode(~, P_vec, A, B, Q, R)
    % Reshape P_vec back to a matrix form
    P = reshape(P_vec, size(A));
    
    % Compute the derivative of P using Riccati Equation
    dP = -(A' * P + P * A - P * B * (R \ B') * P + Q);
    
    % Convert dP back to vector form
    dP_vec = dP(:);
end

% State dynamics function (Forward integration)
function dx = state_dynamics(t, x, t_backward, K_sol, A, B)
    % Interpolate K(t) at the current time step
    K_t = interp1(t_backward, K_sol, t, 'linear');
    
    % Compute control input
    u = - squeeze(K_t) * x;
    
    % Compute state dynamics
    dx = A * x + B * u;
end


%% Ex. 2d
clear; clc; close all;

% Define system matrices
A = [5  0  3  0  1;
     3  0  0 -2  0;
     0 -2  4  1  0;
     1  3 -4  1  3;
     0  2  2  0 -1];

B = [0 1;
     0 2;
     0 0;
     1 3;
     1 1];

% Define the reference trajectory
r = @(t) [0; 0; sin(pi/5 * t); 0; 0];

% Define the cost function parameters
a = 100;
b = 100;
Q = eye(5);
R = [a 0; 0 b];
S = eye(5); 

% Initial and terminal conditions
x0 = [0; 0; 0; 0; 0];  
tf = 500;  % Final time
r_tf = r(tf);  % Final reference

% Final boundary conditions for s and K
K_tf = S;  % Terminal condition for K
s_tf = - S * r_tf;  % Terminal condition for s

% Solve for K(t) backward in time
[t_K, K_sol] = ode45(@(t, K_vec) dKdt(t, K_vec, A, B, Q, R), [tf 0], K_tf(:));
K_sol = reshape(K_sol.', size(A,1), size(A,1), []);
K_sol = permute(K_sol, [3, 1, 2]);
K_sol = flipud(K_sol);  
t_K = flipud(t_K);  

% Solve for s(t) backward in time
[t_s, s_sol] = ode45(@(t, s) dsdt(t, s, A, B, K_sol, Q, r, R, t_K), [tf 0], s_tf);
s_sol = flipud(s_sol); 
t_s = flipud(t_s);

% Solve the state equation forward in time 
[t_x, x_sol] = ode45(@(t, x) state_dynamics2(t, x, A, B, K_sol, s_sol,...
    R, t_K, t_s), [0 tf], x0);

% Control inputs
u_sol = zeros(length(t_x), 2);  
for i = 1:length(t_x)
    K_t = interp1(t_K, K_sol, t_x(i)); 
    s_t = interp1(t_s, s_sol, t_x(i)); 
    u_sol(i, :) = -inv(R) * (B' * squeeze(K_t) * x_sol(i, :)' - B' * s_t');
end

% Compute cost J
cost = zeros(length(t_x), 1);  
error_sol = zeros(length(t_x), 5);  
for i = 1:length(t_x)
    r_t = r(t_x(i));  % Reference at time t
    state_error = x_sol(i, :)' - r_t;  % State deviation from reference
    error_sol(i, :) = state_error';  
    control_input = u_sol(i, :)'; 
    
    % Instantaneous cost
    cost(i) = state_error' * Q * state_error + control_input' * R * control_input;
end

% Integrate cost
J_integration = 0.5 * trapz(t_x, cost);

% Compute the terminal cost
x_tf = x_sol(end, :)';  % State at final time
terminal_cost = 0.5 * (x_tf - r_tf)' * S * (x_tf - r_tf);

% Total cost
J_total = J_integration + terminal_cost;

% Print the total cost function
fprintf('The total cost function J is: %.10f\n', J_total);

% Plot the state variables over time
figure;
plot(t_x, x_sol, 'LineWidth', 1.5);
legend({'$x_1$', '$x_2$', '$x_3$', '$x_4$', '$x_5$'}, ...
    'Interpreter', 'latex', 'Location', 'northeast');
grid on;
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('States [-]', 'Interpreter', 'latex');
exportgraphics(gca, 'HW3_4a.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot the control inputs over time
figure;
plot(t_x, u_sol, 'LineWidth', 1.5);
legend({'$u_1$', '$u_2$'}, 'Interpreter', 'latex', 'Location', 'northeast');
grid on;
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Control Inputs [-]', 'Interpreter', 'latex');
exportgraphics(gca, 'HW3_4b.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot the error between the system and the reference over time
figure;
plot(t_x, error_sol, 'LineWidth', 1.5);
legend({'$e_1$', '$e_2$', '$e_3$', '$e_4$', '$e_5$'}, ...
    'Interpreter', 'latex', 'Location', 'northeast');
grid on;
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Error [-]', 'Interpreter', 'latex');
exportgraphics(gca, 'HW3_4c.pdf', 'ContentType','image',...
    'Resolution', 1000);


function dK = dKdt(~, K_vec, A, B, Q, R)
    % Reshape K_vec back to a matrix form
    K = reshape(K_vec, 5, 5);
    
    % Compute the derivative of K 
    dK = -(K * A + A' * K - K * B * (R \ B') * K + Q);
    
    % Convert dK back to vector form
    dK = dK(:);
end

function ds = dsdt(t, s, A, B, K_sol, Q, r, R, t_K)
    % Interpolate K(t) at the current time step
    K_t = interp1(t_K, K_sol, t, 'linear');
    
    % Compute the reference trajectory
    r_t = r(t);
    
    % Compute the adjoint vector dynamics
    ds = -((A' - squeeze(K_t) * B * (R \ B')) * s - Q * r_t);
end

function dx = state_dynamics2(t, x, A, B, K_sol, s_sol, R, t_K, t_s)
    % Interpolate K(t) and s(t) at the current time step
    K_t = interp1(t_K, K_sol, t);
    s_t = interp1(t_s, s_sol, t)';
    
    % Compute the control input
    u = -inv(R) * (B' * squeeze(K_t) * x - B' * s_t);
    
    % Compute the state dynamics
    dx = A * x + B * u;
end

%% Check Controllability
clear; clc; close;

% Define your matrices A and B
A = [5  0  3  0  1;
     3  0  0 -2  0;
     0 -2  4  1  0;
     1  3 -4  1  3;
     0  2  2  0 -1];

B = [0 1;
     0 2;
     0 0;
     1 3;
     1 1];

% Calculate the controllability matrix
Co = ctrb(A, B);

% Check if the system is fully controllable by comparing the rank
n = size(A, 1);  
rank_of_Co = rank(Co);

if rank_of_Co == n
    disp('The system is fully controllable.');
else
    disp('The system is NOT fully controllable.');
end

%% Check zeros system
clear; clc; close;

% Define system matrices A and B
A = [5  0  3  0  1;
     3  0  0 -2  0;
     0 -2  4  1  0;
     1  3 -4  1  3;
     0  2  2  0 -1];

B = [0 1;
     0 2;
     0 0;
     1 3;
     1 1];

C = [0 0 1 0 0];

D = [0 0];

% Create the state-space system
sys = ss(A, B, C, D);

% Get the transfer function from the state-space representation
sys_tf = tf(sys);

% Get the zeros of the system
z = tzero(sys_tf);

% For continuous-time systems (check if zeros are in the right half-plane)
if isct(sys)
    if any(real(z) > 0)
        disp('The system exhibits non-minimum phase behavior (has zeros in the right half-plane).');
    else
        disp('The system is minimum phase (all zeros are in the left half-plane).');
    end
else  % For discrete-time systems (check if zeros are outside the unit circle)
    if any(abs(z) > 1)
        disp('The system exhibits non-minimum phase behavior (has zeros outside the unit circle).');
    else
        disp('The system is minimum phase (all zeros are inside the unit circle).');
    end
end


