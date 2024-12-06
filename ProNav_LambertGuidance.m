%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:       Giovanni Fereoli
% DATE:         October 10, 2024
% CLASS:        ASEN 6015: Aerospace Vehicle Guidance and Control
% INSTRUCTOR:   Prof. Jay W. Mcmahon
% ASSIGNMENT:   Homework 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Ex. 1a
clear; clc; close all;

% Constants
VM = 3000; % Missile velocity (ft/s)
VT = 1000; % Target velocity (ft/s)
Vc = VM + VT; % Closing velocity (ft/s)
R0 = 40000; % Initial range (ft)
g = 32.2; % Gravity (ft/s^2)
tgo = R0 / Vc; % Time to go (s)

% Simulation parameters
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-7);
disc = 1000;
tspan = linspace(0, tgo, disc);

% Initial conditions
y0 = 0; % Initial lateral displacement (ft)

% Guidance constants
N_values = [3, 4, 5]; 
theta_HE_values = [30, 0, 30]; % Heading error angles (degrees)
nT_values = [0, 3 * g, 3 * g]; % Target accelerations (ft/s^2)

% ProNav commanded acceleration
ProNav = @(N, lambda_dot, Vc) N * Vc * lambda_dot;

% Augmented ProNav commanded acceleration
AugmentedProNav = @(N, lambda_dot, Vc, nT) N * (Vc * lambda_dot + 0.5 * nT);

% Loop through each case
for case_num = 1:3
    theta_HE = theta_HE_values(case_num); 
    nT = nT_values(case_num);
    
    % Initial conditions for state vector [y, y_dot, nT]
    y_dot_0 = VM * deg2rad(theta_HE); % Initial lateral velocity (ft/s)
    x0 = [y0; y_dot_0; nT]; 

    % Allocate arrays for results
    y_proNav = zeros(length(tspan), length(N_values));
    y_augProNav = zeros(length(tspan), length(N_values));
    nc_proNav = zeros(length(tspan), length(N_values));
    nc_augProNav = zeros(length(tspan), length(N_values));

    % Loop 
    for i = 1:length(N_values)
        N = N_values(i);

        % ProNav dynamics
        ode_proNav = @(t, x) [x(2);  
            x(3)-ProNav(N, (x(1) + x(2) * (tgo - t))/(Vc * (tgo - t)^2),...
            Vc); 
            0]; 
        
        % Augmented ProNav dynamics
        ode_augProNav = @(t, x) [x(2);  
            x(3)-AugmentedProNav(N, ...
            (x(1) + x(2) * (tgo - t))/(Vc * (tgo - t)^2), Vc, x(3)); 
            0];  

        % Solve for ProNav
        [T_proNav, X_proNav] = ode45(ode_proNav, tspan, x0, options);

        % Solve for Augmented ProNav
        [T_augProNav, X_augProNav] = ode45(ode_augProNav, tspan, x0, options);

        % Store results for ProNav
        y_proNav(:, i) = X_proNav(:, 1);

        % Store results for Augmented ProNav
        y_augProNav(:, i) = X_augProNav(:, 1);

        % Store nc for plotting
        lambda_dot_ProNav = (X_proNav(:, 1) + X_proNav(:, 2) .* ...
            (tgo - T_proNav)) ./ (Vc .* (tgo - T_proNav).^2);
        nc_proNav(:, i) = ProNav(N, lambda_dot_ProNav, Vc);
        lambda_dot_augProNav = (X_augProNav(:, 1) + X_augProNav(:, 2) .* ...
            (tgo - T_augProNav)) ./ (Vc .* (tgo - T_augProNav).^2);
        nc_augProNav(:, i) = AugmentedProNav(N, lambda_dot_augProNav,...
            Vc, X_augProNav(:, 3));
    end

    % Plot results for ProNav
    gca1 = figure(1);
    subplot(3, 1, case_num);
    plot(T_proNav, y_proNav(:, 1), 'r', T_proNav, y_proNav(:, 2), ...
        'g', T_proNav, y_proNav(:, 3), 'b', 'LineWidth', 1.5);
    title(['\textbf{ProNav:} $y$ \textbf{vs time, Case} ',...
        num2str(case_num)], 'Interpreter', 'latex');
    xlabel('\textbf{Time [s]}', 'Interpreter', 'latex');
    ylabel('$y$ [ft]', 'Interpreter', 'latex');
    legend('$N=3$', '$N=4$', '$N=5$', 'Interpreter', 'latex');
    grid on;
    
    % Plot results for Augmented ProNav
    gca2 = figure(2);
    subplot(3, 1, case_num);
    plot(T_augProNav, y_augProNav(:, 1), 'r', T_augProNav,...
        y_augProNav(:, 2), 'g', T_augProNav, y_augProNav(:, 3),...
        'b', 'LineWidth', 1.5);
    title(['\textbf{Augmented ProNav:} $y$ \textbf{vs time, Case} ',...
        num2str(case_num)], 'Interpreter', 'latex');
    xlabel('\textbf{Time [s]}', 'Interpreter', 'latex');
    ylabel('$y$ [ft]', 'Interpreter', 'latex');
    legend('$N=3$', '$N=4$', '$N=5$', 'Interpreter', 'latex');
    grid on;
    
    % Plot commanded accelerations 
    gca3 = figure(3);
    subplot(3, 1, case_num);
    plot(T_proNav, nc_proNav(:, 1), 'r', T_proNav, nc_proNav(:, 2),...
        'g', T_proNav, nc_proNav(:, 3), 'b', 'LineWidth', 1.5);
    title(['\textbf{ProNav:} $n_c$ \textbf{vs time, Case} ',...
        num2str(case_num)], 'Interpreter', 'latex');
    xlabel('\textbf{Time [s]}', 'Interpreter', 'latex');
    ylabel('$n_c$ [ft/s$^2$]', 'Interpreter', 'latex');
    legend('$N=3$', '$N=4$', '$N=5$', 'Interpreter', 'latex');
    grid on;
    
    % Plot augmented commanded accelerations 
    gca4 = figure(4);
    subplot(3, 1, case_num);
    plot(T_augProNav, nc_augProNav(:, 1), 'r', T_augProNav, ...
        nc_augProNav(:, 2), 'g', T_augProNav, nc_augProNav(:, 3), 'b', 'LineWidth', 1.5);
    title(['\textbf{Augmented ProNav:} $n_c$ \textbf{vs time, Case} ',...
        num2str(case_num)], 'Interpreter', 'latex');
    xlabel('\textbf{Time [s]}', 'Interpreter', 'latex');
    ylabel('$n_c$ [ft/s$^2$]', 'Interpreter', 'latex');
    legend('$N=3$', '$N=4$', '$N=5$', 'Interpreter', 'latex');
    grid on;

end

% Export graphics
exportgraphics(gca1, 'HW2_1a_a.pdf', 'ContentType','image',...
    'Resolution', 1000);
exportgraphics(gca2, 'HW2_1a_b.pdf', 'ContentType','image',...
    'Resolution', 1000);
exportgraphics(gca3, 'HW2_1a_c.pdf', 'ContentType','image',...
    'Resolution', 1000);
exportgraphics(gca4, 'HW2_1a_d.pdf', 'ContentType','image',...
    'Resolution', 1000);


%% Ex. 1b
close all; clc; clear;

% Constants
VM = 3000; 
VT = 1000; 
Vc = VM + VT; 
R0 = 40000; 
g = 32.2; 

% Simulation parameters
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-7); 
disc = 1000; % Number of time steps
MC_runs = 1000; % Number of Monte Carlo runs
sigma_tgo = sqrt(0.3); % Standard deviation for tgo error (s)

% Case details 
theta_HE = 30;
nT = 3 * g; 

% Initial conditions for state vector 
y0 = 0; 
y_dot_0 = VM * deg2rad(theta_HE); 
x0 = [y0; y_dot_0; nT];

% Time to go 
tgo_nominal = R0 / Vc; 
tspan = linspace(0, tgo_nominal, disc); 

% Augmented ProNav parameters
N_augProNav = 4;

% Augmented ProNav commanded acceleration
AugmentedProNav = @(N, lambda_dot, Vc, nT) N * (Vc * lambda_dot + 0.5 * nT);

% Preallocate for Monte Carlo results
y_MC = zeros(disc, MC_runs);
nc_MC = zeros(disc, MC_runs);
y_final = zeros(1, MC_runs); 

% Monte Carlo Simulation
for k = 1:MC_runs
    % Introduce random error in tgo
    % NOTE: ZEM should be computed considering t_max = tgo_nominal
    tgo_perturbed = normrnd(tgo_nominal, sigma_tgo); % Perturbed tgo
    tspan_perturbed = linspace(0, min(10, tgo_perturbed), disc); % Perturbed time span
    
    % Augmented ProNav with perturbed tgo
    ode_augProNav = @(t, x) [x(2);  
        x(3)-AugmentedProNav(N_augProNav, ...
        (x(1) + x(2) * (tgo_perturbed  - t))/(Vc * (tgo_perturbed - t)^2),...
        Vc, x(3)); 
        0];  
    
    % Solve 
    [T_augProNav, X_augProNav] = ode45(ode_augProNav, ...
        tspan_perturbed, x0, options);
    
    % Store results for each Monte Carlo run
    % OSS: Interpolate to match nominal time span
    y_MC(:, k) = interp1(T_augProNav, X_augProNav(:, 1), tspan); 
    y_MC(isnan(y_MC(:, k)), k) = 0; % Replace NaNs with 0 in y_MC
    nc_MC(:, k) = interp1(T_augProNav,...
        AugmentedProNav(N_augProNav,...
        (X_augProNav(:, 1) + X_augProNav(:, 2) .* ...
        (tgo_perturbed * ones(disc,1)  - T_augProNav)) ./ (Vc .* ...
        (tgo_perturbed * ones(disc,1) - T_augProNav).^2), Vc,...
         X_augProNav(:, 3)), tspan);
    nc_MC(isnan(nc_MC(:, k)), k) = 0; % Replace NaNs with 0 in nc_MC
    y_final(k) = X_augProNav(end, 1); % Store the final value of y
end

% Compute mean and standard deviation of results
y_mean = mean(y_MC, 2);
y_std = std(y_MC, 0, 2);
nc_mean = mean(nc_MC, 2);
nc_std = std(nc_MC, 0, 2);

% Plotting results of Monte Carlo simulation
y_std_upper_positive = y_mean + 3 * y_std;
y_std_lower_positive = y_mean - 3 * y_std;
y_std_lower_positive(y_std_lower_positive  < 0) = NaN;
gca1 = figure(1);
subplot(2, 1, 1);
plot(tspan, y_mean, 'r', 'LineWidth', 1.5); hold on;
plot(tspan, y_std_upper_positive, 'b--', 'LineWidth', 1);
plot(tspan, y_std_lower_positive, 'b--', 'LineWidth', 1);
title('\textbf{Monte Carlo Simulation:} $y$ \textbf{vs time, Augmented ProNav,} $N = 4$', ...
    'Interpreter', 'latex');
xlabel('\textbf{Time [s]}', 'Interpreter', 'latex');
ylabel('$y$ [ft]', 'Interpreter', 'latex');
legend('Mean $y$', '$\pm 3 \sigma$', 'Interpreter', 'latex');
grid on;
subplot(2, 1, 2);
plot(tspan, nc_mean, 'r', 'LineWidth', 1.5); hold on;
plot(tspan, nc_mean + 3 * nc_std, 'b--', 'LineWidth', 1);
plot(tspan, nc_mean - 3 * nc_std, 'b--', 'LineWidth', 1);
title('\textbf{Monte Carlo Simulation:} $n_c$ \textbf{vs time, Augmented ProNav,} $N = 4$',...
    'Interpreter', 'latex');
xlabel('\textbf{Time [s]}', 'Interpreter', 'latex');
ylabel('$n_c$ [ft/s$^2$]', 'Interpreter', 'latex');
legend('Mean $n_c$', '$\pm 3 \sigma$', 'Interpreter', 'latex');
grid on;

% Monte Carlo: Final y distribution and Gaussian fit
gca2 = figure(2);
histogram(y_final, 50, 'Normalization', 'pdf'); % Histogram of final y
hold on;
pd = fitdist(y_final', 'Normal'); % Fit a Gaussian distribution
x_values = linspace(min(y_final), max(y_final), 1000);
y_gauss = pdf(pd, x_values');
plot(x_values, y_gauss, 'r', 'LineWidth', 2); % Plot Gaussian fit
xlabel('$y$ [ft]', 'Interpreter', 'latex');
ylabel('Probability Density [-]', 'Interpreter', 'latex');
legend('Final $y$', 'Gaussian Fit', 'Interpreter', 'latex', 'Location', 'best');
grid on;

% Display the results
disp('Gaussian Fit, ZEM Mean [ft]:');
disp(pd.mu);
disp('Gaussian Fit, ZEM Standard Deviation [ft]:');
disp(pd.sigma);

% Export graphics
exportgraphics(gca1, 'HW2_1b_a.pdf', 'ContentType','image',...
    'Resolution', 1000);
exportgraphics(gca2, 'HW2_1b_b.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Ex. 2
clear; clc; close all;

% Define problem parameters
mu_earth = 398600.4418;  % [km^3/s^2]
r0 = [6578; 0; 0];  % Initial position [km]
v0 = [0; sqrt(mu_earth / norm(r0)); 0];  % [km/s] 
rf = [-42164 * cosd(5); 42164 * sind(5); 0];  % Final position [km]
dt = 37864;  % Transfer time [seconds]
DM = 1;  % Assuming prograde 

% Call the lambert solver function
[v1, v2, errorl, count] = lambert(r0, rf, dt, DM);

% Check if the solution converged
if strcmp(errorl, 'ok')
    disp('Lambert solution found successfully:');
    disp(['Number of iterations: ', num2str(count)]);
    disp('Initial velocity (v1) [km/s]:');
    disp(v1);
    disp('Final velocity (v2) [km/s] (velocity at rf):');
    disp(v2);
else
    disp('Error in solving Lambert problem:');
    disp(errorl);
    return;
end

% Two-body dynamics function
function dx = two_body_dynamics(~, x, mu)
    r = x(1:3); 
    v = x(4:6);
    r_norm = norm(r);
    a = - mu / r_norm^3 * r;
    dx = [v; a];
end

% Define conditions and integrate
initial_conditions = [r0; v0];  
transfer_conditions = [r0; v1]; 
final_conditions = [rf; v2];    
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[t_circ, state_circ] = ode45(@(t, y) two_body_dynamics(t, y, mu_earth),...
    [0 dt], initial_conditions, options);
[t_transfer, state_transfer] = ode45(@(t, y) two_body_dynamics(t, y,...
    mu_earth), [0 dt], transfer_conditions, options);
[t_final, state_final] = ode45(@(t, y) two_body_dynamics(t, y, mu_earth),...
    [0 2 * dt], final_conditions, options);

% Plot the orbits and the transfer
gca = figure;
hold on;
plot3(state_circ(:,1), state_circ(:,2), state_circ(:,3), 'b', ...
    'LineWidth', 1.5); % Initial orbit
plot3(state_final(:,1), state_final(:,2), state_final(:,3), 'r', ...
    'LineWidth', 1.5); % Final orbit
plot3(state_transfer(:,1), state_transfer(:,2), state_transfer(:,3), ...
    'k--', 'LineWidth', 2); % Transfer path
plot3(r0(1), r0(2), r0(3), 'kd', 'MarkerSize', 10, 'MarkerFaceColor', 'g'); 
plot3(rf(1), rf(2), rf(3), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  
plot_earth();
grid on;
axis equal;
xlabel('x $\cdot 10^{4}$ [km]', 'Interpreter', 'latex');
ylabel('y $\cdot 10^{4}$ [km]', 'Interpreter', 'latex');
zlabel('z $\cdot 10^{4}$ [km]', 'Interpreter', 'latex');
legend('Initial Orbit', 'Final Orbit', 'Transfer Path', '$\mathbf{r}_0$',...
    '$\mathbf{r}_f$', 'Earth', 'Interpreter', 'latex');
hold off;
exportgraphics(gca, 'HW2_2a.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Lambert solver function
function [v1, v2, errorl, count] = lambert(r1, r2, dtsec, DM)
    % Initialization
    mu = 398600.4418;  % Gravitational parameter [km^3/s^2]
    tol = 1e-5;        % Tolerance for convergence
    num_iter = 10000;  % Max number of iterations
    errorl = 'ok';     % Error flag

    % Magnitudes of r1 and r2
    r1_norm = norm(r1);
    r2_norm = norm(r2);
    
    % Calculate cos of the true anomaly change
    cos_delta_theta = dot(r1, r2) / (r1_norm * r2_norm);
    A = DM * sqrt(r1_norm * r2_norm * (1 + cos_delta_theta));

    % Check for impossible cases
    if cos_delta_theta == 0 || A == 0
        errorl = 'Error: Cannot calculate with given values.';
        v1 = []; v2 = [];
        return;
    end
    
    % Set initial bounds (assuming single revolution)
    lower = -4.0 * pi^2;  
    upper = 4.0 * pi^2;

    % Initial guess for psi
    psi = 0;
    count = 0;
    dtnew = 0; % Initializing dtnew to start the loop

    % Initial values for universal variable functions
    c2 = 1/2;
    c3 = 1/6;

    % Main iteration loop
    while abs(dtnew - dtsec) > tol && count < num_iter
        % Compute y
        y = compute_y(c2, c3, A, psi, r1_norm, r2_norm);
        
        % Handle y < 0 
        if A > 0 && y < 0
            [psi, y] = handle_negative_y(psi, A, r1_norm, r2_norm, c2, c3);
        end

        % Compute Chi
        if abs(c2) > tol
            Chi = sqrt(y / c2);
        else
            Chi = 0.0;
        end
        
        % Compute time
        dtnew = (Chi^3 * c3 + A * sqrt(y)) / sqrt(mu);
        
        % Update psi using Bisection method
        if dtnew <= dtsec
            lower = psi;
        else
            upper = psi;
        end 
        psi = (upper + lower) * 0.5;

        % Update c2 and c3 using the universal variable auxiliary functions
        [c2, c3] = find_c2_c3(psi);

        % Update iteration count
        count = count + 1;
    end

    % Check if the iteration converged
    if abs(dtnew - dtsec) > tol
        errorl = 'Error: Not converged within tolerance.';
        v1 = []; v2 = [];
        return;
    end
    
    % Calculate velocity vectors
    f = 1 - y / r1_norm;
    gdot = 1 - y / r2_norm;
    g = 1 / (A * sqrt(y / mu));
    
    % Compute velocities
    v1 = (r2 - f * r1) * g;
    v2 = (gdot * r2 - r1) * g;
end

% Computes y for Lambert problem using universal variable functions.
function y = compute_y(c2, c3, A, psi, r1_norm, r2_norm)
    if abs(c2) > 1e-5
        y = r1_norm + r2_norm - (A * (1.0 - psi * c3) / sqrt(c2));
    else
        y = r1_norm + r2_norm;
    end
end

% Handles negative y by adjusting psi
function [lower, y] = handle_negative_y(A, r1_norm, r2_norm, c2, c3)
    count = 1;
    while y < 0 && count < 100
        lower = 0.8 * (1.0 / c3) * (1.0 - (r1_norm + r2_norm) * sqrt(c2) / A);
        [c2, c3] = find_c2_c3(lower);
        y = compute_y(c2, c3, A, lower, r1_norm, r2_norm);
        count = count + 1;
    end
end

% Computes c2 and c3 based on psi using universal variable functions.
function [c2, c3] = find_c2_c3(psi)
    if psi > 1e-5
        sqrtpsi = sqrt(psi);
        c2 = (1.0 - cos(sqrtpsi)) / psi;
        c3 = (sqrtpsi - sin(sqrtpsi)) / (sqrtpsi^3);
    elseif psi < -1e-5
        sqrtnegpsi = sqrt(-psi);
        c2 = (1.0 - cosh(sqrtnegpsi)) / psi;
        c3 = (sinh(sqrtnegpsi) - sqrtnegpsi) / (-sqrtnegpsi^3);
    else
        c2 = 0.5;
        c3 = 1.0 / 6.0;
    end
end

%% Ex. 2b
close all; clc; clear;

% Define problem parameters
mu_earth = 398600.4418; 
r0 = [6578; 0; 0];  
v0 = [0; sqrt(mu_earth / norm(r0)); 0]; 
rf = [-42164 * cosd(5); 42164 * sind(5); 0]; 
dt = 37864; 
DM = 1;  

% Thrust acceleration 
a_T = 0.03;  % [km/s^2]

% Tolerance for achieving desired velocity
tol = 1e-12;  % [km/s]

% Time step for guidance updates
dt_guidance = 1;  % [seconds]

% Initialize state and time
state = [r0; v0];  
current_time = 0;  

% Set up ODE solver options
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% History storage for plotting
position_history = [];
velocity_history = [];
thrust_history = []; 
time_history = [];  

% Run simulation with closed-loop Lambert guidance
while current_time < dt
    % Call Lambert solver to get desired velocity at this point in time
    [v_desired, ~, errorl, ~] = lambert(state(1:3), rf, dt - current_time, DM);
    
    if strcmp(errorl, 'ok')
        % Compute thrust direction based on desired velocity
        thrust_direction = (v_desired - state(4:6)) / norm(v_desired - state(4:6));
    else
        disp('Error in Lambert solver:');
        disp(errorl);
        return;
    end
    
    % Define the dynamics function for 2BP with thrust
    dynamics = @(t, y) two_body_with_thrust(t, y, mu_earth, a_T,...
        thrust_direction, v_desired, tol);
    
    % Integrate the dynamics for the current guidance time step
    [t, state_new] = ode113(dynamics,...
        [current_time current_time + dt_guidance], state, options);
    
    % Update the state to the last value after integration
    state = state_new(end, :)';
    
    % Update time
    current_time = current_time + dt_guidance;
    
    % Store history for plotting
    position_history = [position_history; state(1:3)'];
    thrust_history = [thrust_history; a_T * thrust_direction'];
    time_history = [time_history; current_time];
    
    % Check if the burn is complete 
    if norm(state(4:6) - v_desired) < tol
        disp('Burn complete: Target velocity achieved.');
        break;
    end
end

% Plot the trajectory results
gca1 = figure(1);
plot3(position_history(:,1), position_history(:,2),...
    position_history(:,3), 'b', 'LineWidth', 2);  % Trajectory
hold on;
plot3(r0(1), r0(2), r0(3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g'); 
plot3(rf(1), rf(2), rf(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); 
plot3(position_history(end,1), position_history(end,2),...
    position_history(end,3), 'kx', 'MarkerSize', 12, 'LineWidth', 2);  
plot_earth();  
axis equal;
xlabel('x $\cdot 10^{4}$ [km]', 'Interpreter', 'latex');
ylabel('y $\cdot 10^{4}$ [km]', 'Interpreter', 'latex');
zlabel('z $\cdot 10^{4}$ [km]', 'Interpreter', 'latex');
legend('Trajectory', '$\mathbf{r}_0$', '$\mathbf{r}_f$',...
    '$\mathbf{r}_{final}$', 'Earth', 'Interpreter', 'latex');
grid on;
view(2);
hold off;

% Display the difference between final and desired position
position_diff = position_history(end,:)' - rf;
disp(['Difference between final and desired position: [', num2str(position_diff(1)), ', ',...
    num2str(position_diff(2)), ', ', num2str(position_diff(3)), '] km']);
disp(['Magnitude of the difference: ', num2str(norm(position_diff)), ' km']);

% Plot the control (thrust) history
gca2 = figure(2);
plot(time_history / 3600, thrust_history(:,1), 'r', 'LineWidth', 1);  
hold on;
plot(time_history / 3600, thrust_history(:,2), 'b', 'LineWidth', 1);  
plot(time_history / 3600, thrust_history(:,3), 'g', 'LineWidth', 1);
xlim([time_history(1) / 3600, time_history(end) / 3600]);
xlabel('Time [hr]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Thrust [$\mathrm{km/s^2}$]', 'Interpreter', 'latex', 'FontSize', 12);
legend('$u_x$', '$u_y$', '$u_z$', 'Location', 'northeast', 'Interpreter', 'latex');
grid on;
grid minor;
hold off;

% Export graphics
exportgraphics(gca1, 'HW2_2b_a.pdf', 'ContentType','image',...
    'Resolution', 1000);
exportgraphics(gca2, 'HW2_2b_b.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Dynamics function for two-body problem with thrust
function dydt = two_body_with_thrust(~, y, mu, a_T, thrust_dir, v_desired, tol)
    % Extract from y
    r = y(1:3);  
    v = y(4:6);  
    
    % Two-body problem 
    r_norm = norm(r);
    a_gravity = -mu / r_norm^3 * r;
    
    % Thrust acceleration
    if norm(v_desired - v) > tol
        a_thrust = a_T * thrust_dir;
    else
        a_thrust = [0; 0; 0];  % No thrust once the velocity is close to desired
    end
    
    % Total acceleration
    a_total = a_gravity + a_thrust;
    
    % Return the derivatives 
    dydt = [v; a_total];
end

%% Ex. 2c
clear; clc; close;

% Define problem parameters
mu_earth = 398600.4418;  
r0 = [6578; 0; 0];  
v0 = [0; sqrt(mu_earth / norm(r0)); 0];  
rf = [-42164 * cosd(5); 42164 * sind(5); 0]; 
dt = 37864;  
DM = 1;  
t_final = dt;  

% Thrust acceleration  
a_T = 0.03;  % [km/s^2]

% Tolerance for achieving desired velocity
tol = 1e-12;  

% Time step for guidance updates
dt_guidance = 1; 

% Initialize state and time
state = [r0; v0]; 
current_time = 0; 

% Set up ODE solver options
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% History storage for plotting
position_history = [];
velocity_history = [];
thrust_history = [];
time_history = [];

% Nominal Lambert transfer for gravity gradient calculation
[v1_nominal, v2_nominal, errorl, count] = lambert(r0, rf, dt, DM);
if strcmp(errorl, 'ok')
    disp('Lambert nominal solution found successfully.');
else
    disp('Error in Lambert solver:');
    disp(errorl);
    return;
end

% Gravity-gradient matrix G computation function
function G = gravity_gradient_matrix(r, mu)
    r_norm = norm(r);
    G = mu / r_norm^5 * (3 * (r * r') - r_norm^2 * eye(3));
end

% Cross-product steering with Q matrix
while current_time < dt
    % Call Lambert solver to get desired velocity at this point in time
    [v_desired, ~, errorl, ~] = lambert(state(1:3), rf, t_final - current_time, DM);
    
    if strcmp(errorl, 'ok')
        % Compute gravity-gradient matrix G along the trajectory
        G = gravity_gradient_matrix(state(1:3), mu_earth);
        
        % Compute the Q matrix
        Q = -1 / (dt_guidance) * (eye(3) + (dt_guidance)^2 / 2 * G);
        
        % Compute thrust direction based on Q and desired velocity
        v_diff = v_desired - state(4:6);      
        v_diff_dir = v_diff / norm(v_diff);  
        
        b_vec = - Q * v_diff;                
        b_vec_dot = dot(b_vec, v_diff_dir);  
        
        % Compute q using the quadratic formula approach
        q = -b_vec_dot + sqrt(b_vec_dot^2 - norm(b_vec)^2 + a_T^2);
        
        % Calculate and normalize the thrust direction
        thrust_command = b_vec + q * v_diff_dir;
    else
        disp('Error in Lambert solver:');
        disp(errorl);
        return;
    end
    
    % Define the dynamics function for 2BP with cross-product steering using Q
    dynamics = @(t, y) two_body_with_thrust2(t, y, mu_earth,...
        thrust_command, v_desired, tol);
    
    % Integrate the dynamics for the current guidance time step
    [t, state_new] = ode45(dynamics, ...
        [current_time current_time + dt_guidance], state, options);
    
    % Update the state to the last value after integration
    state = state_new(end, :)';
    
    % Update time
    current_time = current_time + dt_guidance;
    
    % Store the position and velocity history for plotting
    position_history = [position_history; state(1:3)'];
    velocity_history = [velocity_history; state(4:6)'];

    % Store the thrust magnitude
    thrust_history = [thrust_history; thrust_command'];
    
    % Store the time
    time_history = [time_history; current_time];
    
    % Check if the burn is complete 
    if norm(state(4:6) - v_desired) < tol
        disp('Burn complete: Target velocity achieved.');
        break;
    end
end

% Plot the trajectory results
gca1 = figure(1);
plot3(position_history(:,1), position_history(:,2),...
    position_history(:,3), 'b', 'LineWidth', 2);  
hold on;
plot3(r0(1), r0(2), r0(3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');  
plot3(rf(1), rf(2), rf(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); 
plot3(position_history(end,1), position_history(end,2),...
    position_history(end,3), 'kx', 'MarkerSize', 12, 'LineWidth', 2);  
plot_earth();  
axis equal;
xlabel('x $\cdot 10^{4}$ [km]', 'Interpreter', 'latex');
ylabel('y $\cdot 10^{4}$ [km]', 'Interpreter', 'latex');
zlabel('z $\cdot 10^{4}$ [km]', 'Interpreter', 'latex');
legend('Trajectory', '$\mathbf{r}_0$', '$\mathbf{r}_f$',...
    '$\mathbf{r}_{final}$', 'Earth', 'Interpreter', 'latex');
grid on;
view(2);
hold off;

% Display the difference between final and desired position
position_diff = position_history(end,:)' - rf;
disp(['Difference between final and desired position: [', num2str(position_diff(1)), ', ',...
    num2str(position_diff(2)), ', ', num2str(position_diff(3)), '] km']);
disp(['Magnitude of the difference: ', num2str(norm(position_diff)), ' km']);

% Plot the control (thrust) history
gca2 = figure(2);
plot(time_history / 3600, thrust_history(:,1), 'r', 'LineWidth', 1);  
hold on;
plot(time_history / 3600, thrust_history(:,2), 'b', 'LineWidth', 1);  
plot(time_history / 3600, thrust_history(:,3), 'g', 'LineWidth', 1);
xlim([time_history(1) / 3600, time_history(end) / 3600])
xlabel('Time [hr]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Thrust [$\mathrm{km/s^2}$]', 'Interpreter', 'latex', 'FontSize', 12);
legend('$u_x$', '$u_y$', '$u_z$', 'Location', 'northeast', 'Interpreter', 'latex');
grid on;
grid minor;
hold off;

% Export graphics
exportgraphics(gca1, 'HW2_2c_a.pdf', 'ContentType','image',...
    'Resolution', 1000);
exportgraphics(gca2, 'HW2_2c_b.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Dynamics function for two-body problem with thrust
function dydt = two_body_with_thrust2(~, y, mu, thrust_command, v_desired, tol)
    % Extract y
    r = y(1:3); 
    v = y(4:6);  
    
    % Two-body problem 
    r_norm = norm(r);
    a_gravity = -mu / r_norm^3 * r;
    
    % Thrust acceleration 
    if norm(v_desired - v) > tol
        a_thrust = thrust_command;
    else
        a_thrust = [0; 0; 0]; 
    end
    
    % Total acceleration
    a_total = a_gravity + a_thrust;
    
    % Return the derivatives
    dydt = [v; a_total];
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

