clc; clear; close all

addpath('../quadrotor_sim/src/')

g = 9.81;  % The gravitational acceleration [m/s^2]
l = 0.2;  % Distance from the center of mass to each rotor [m]
m = 0.5;  % Total mass of the quadrotor [kg]
I = [1.24, 1.24, 2.48];  % Mass moment of inertia [kg m^2]
mu = 3.0;  % Maximum thrust of each rotor [N]
sigma = 0.01;  % The proportionality constant relating thrust to torque [m]

quad = quadrotor(g, l, m, diag(I), mu, sigma);

%% Simple Go-To Point Controller
speed = 1/2;
desired_x = @(t) [cos(t*speed) sin(t*speed) 2+sin(t*speed)];
zd = @(t) [desired_x(t) zeros([1,9])]';

% State and Control Values to Linearize Around
zs = zeros(12, 1);
zs(3) = 1;
us = m*g/4 * ones([4,1]); % enough to counteract gravity

% Linearize the system
[A, B] = linearize(zs, us);
C = ctrb(A,B);

% Check if the system is controllable
if rank(C) < min(size(C))
    fprintf("This linearized system is not controllable!")
    return
else
    fprintf("This linearized system is controllable :)\n")
end

% Cost values for LQR
x1_gain = 10;
x2_gain = 10;
x3_gain = 20;
a1_gain = 100;
a2_gain = 100;
a3_gain = 200;
v1_gain = 1;
v2_gain = 1;
v3_gain = 2;
w1_gain = 10;
w2_gain = 10;
w3_gain = 20;



v_gain = 1/10;
Q = diag([x1_gain x2_gain x3_gain a1_gain a2_gain a3_gain ...
          v1_gain v2_gain v3_gain w1_gain w2_gain w3_gain]);
u_gain = 1/mu^2;
R = eye(4,4) * u_gain;
[K,S,P] = lqr(A,B,Q,R);

u = @(t, z) us+K*(zd(t)-z);

%% Run Simulation
z0 = zeros(12,1);
tspan = [0, 10];

[t,z] = quad.solve(tspan, z0, u);

% Results
quad.animate(t,z)
quad.plot(t,z)
final = z(end, 1:3)
error = desired_x(max(tspan)) - final
