%% 2-DOF Planar Robot Arm Control via Feedback Linearization
% System parameters
L1 = 1.0;   % [m] length of link 1
L2 = 0.8;   % [m] length of link 2
m1 = 2.0;   % [kg] mass of link 1 (concentrated at end)
m2 = 1.5;   % [kg] mass of link 2 (concentrated at end)
g = 9.81;   % [m/s^2] gravity

% PD control gains
Kp = diag([100, 100]);  % Position gain
Kv = diag([20, 20]);    % Velocity gain

% Simulation time
tspan = [0 10];  % [s]
dt = 0.01;       % [s] time step
t = tspan(1):dt:tspan(2);

%% Desired trajectory (sinusoidal)
theta_d = [sin(t); 0.5*cos(t)];              % Desired positions
theta_dot_d = [cos(t); -0.5*sin(t)];         % Desired velocities
theta_ddot_d = [-sin(t); -0.5*cos(t)];       % Desired accelerations

%% Initial conditions
theta0 = [0; 0];         % Initial position
theta_dot0 = [0; 0];     % Initial velocity

%% Simulate the system
state = [theta0; theta_dot0];  % Initial state [theta1; theta2; theta1_dot; theta2_dot]

% Preallocate storage
theta = zeros(2, length(t));
theta_dot = zeros(2, length(t));
tau = zeros(2, length(t));

for i = 1:length(t)
    % Current state
    theta(:,i) = state(1:2);
    theta_dot(:,i) = state(3:4);
    
    % Compute tracking errors
    e = theta_d(:,i) - theta(:,i);
    e_dot = theta_dot_d(:,i) - theta_dot(:,i);
    
    % Feedback linearization control law
    u = theta_ddot_d(:,i) + Kv*e_dot + Kp*e;
    
    % Compute dynamics matrices
    [M, C, G] = computeDynamics(theta(:,i), theta_dot(:,i), L1, L2, m1, m2, g);
    
    % Compute required torque
    tau(:,i) = M*u + C*theta_dot(:,i) + G;
    
    % Simulate forward (using Euler integration for simplicity)
    theta_ddot = M\(tau(:,i) - C*theta_dot(:,i) - G);
    state = state + [theta_dot(:,i); theta_ddot]*dt;
end

%% Plot results
figure('Position', [100, 100, 1200, 800]);

% Joint positions
subplot(3,2,1);
plot(t, theta(1,:), 'b', t, theta_d(1,:), 'r--');
title('Joint 1 Position');
xlabel('Time [s]'); ylabel('\theta_1 [rad]');
legend('Actual', 'Desired');
grid on;

subplot(3,2,2);
plot(t, theta(2,:), 'b', t, theta_d(2,:), 'r--');
title('Joint 2 Position');
xlabel('Time [s]'); ylabel('\theta_2 [rad]');
legend('Actual', 'Desired');
grid on;

% Joint velocities
subplot(3,2,3);
plot(t, theta_dot(1,:), 'b', t, theta_dot_d(1,:), 'r--');
title('Joint 1 Velocity');
xlabel('Time [s]'); ylabel('d\theta_1/dt [rad/s]');
legend('Actual', 'Desired');
grid on;

subplot(3,2,4);
plot(t, theta_dot(2,:), 'b', t, theta_dot_d(2,:), 'r--');
title('Joint 2 Velocity');
xlabel('Time [s]'); ylabel('d\theta_2/dt [rad/s]');
legend('Actual', 'Desired');
grid on;

% Control torques
subplot(3,2,5);
plot(t, tau(1,:), 'b');
title('Joint 1 Torque');
xlabel('Time [s]'); ylabel('\tau_1 [Nm]');
grid on;

subplot(3,2,6);
plot(t, tau(2,:), 'b');
title('Joint 2 Torque');
xlabel('Time [s]'); ylabel('\tau_2 [Nm]');
grid on;

%% Robot animation
figure;
for i = 1:20:length(t)  % Animate every 20 steps
    % Forward kinematics
    x1 = L1*cos(theta(1,i));
    y1 = L1*sin(theta(1,i));
    x2 = x1 + L2*cos(theta(1,i)+theta(2,i));
    y2 = y1 + L2*sin(theta(1,i)+theta(2,i));
    
    % Plot
    plot([0 x1 x2], [0 y1 y2], 'b-o', 'LineWidth', 2);
    hold on;
    plot(x2, y2, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot(theta_d(1,1:i), theta_d(2,1:i), 'g--');
    hold off;
    
    axis equal;
    xlim([-2 2]); ylim([-2 2]);
    title(sprintf('Robot Arm Animation (t=%.2f s)', t(i)));
    xlabel('X [m]'); ylabel('Y [m]');
    grid on;
    drawnow;
end

%% Dynamics computation function
function [M, C, G] = computeDynamics(theta, theta_dot, L1, L2, m1, m2, g)
    theta1 = theta(1);
    theta2 = theta(2);
    theta1_dot = theta_dot(1);
    theta2_dot = theta_dot(2);
    
    % Inertia matrix
    M11 = (m1 + m2)*L1^2 + m2*L2^2 + 2*m2*L1*L2*cos(theta2);
    M12 = m2*L2^2 + m2*L1*L2*cos(theta2);
    M22 = m2*L2^2;
    M = [M11 M12; M12 M22];
    
    % Coriolis matrix
    C11 = -2*m2*L1*L2*sin(theta2)*theta2_dot;
    C12 = -m2*L1*L2*sin(theta2)*theta2_dot;
    C21 = m2*L1*L2*sin(theta2)*theta1_dot;
    C22 = 0;
    C = [C11 C12; C21 C22];
    
    % Gravity vector
    G1 = (m1 + m2)*g*L1*cos(theta1) + m2*g*L2*cos(theta1 + theta2);
    G2 = m2*g*L2*cos(theta1 + theta2);
    G = [G1; G2];
end