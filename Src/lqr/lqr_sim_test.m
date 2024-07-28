
clc; clear; close all

addpath('../quadrotor_sim/src');

% QUADROTOR

g = 9.81;  % The gravitational acceleration [m/s^2]
l = 0.2;  % Distance from the center of mass to each rotor [m]
m = 0.5;  % Total mass of the quadrotor [kg]
I = [1.24, 1.24, 2.48];  % Mass moment of inertia [kg m^2]
mu = 3.0;  % Maximum thrust of each rotor [N]
sigma = 0.01;  % The proportionality constant relating thrust to torque [m]

quad = quadrotor(g, l, m, diag(I), mu, sigma);
dt = 0.01;

% INTRUDER
speed = 0.5;
r = 3;
path = @(t) [r*cos(t*speed); r*sin(t*speed); 5];
% path = @(t) [-5 + speed*t; 2.5; 5];
% path = @(t) [5*cos(t*speed); 5*sin(t*speed); 5];
% path = @(t) [-t; t; 3];
% path = @(t) [t; 5*(t-10-5)*(t+10-5)/-100; 5];
% dist = struct("r", @(t,z)0.0*[0.5*sin(t); 0.5*sin(2*t); 0.5*sin(4*t)],...
%     "n", @(t,z)[0.1; 0.01; 0.1]);
dist = struct("r", @(t,z)0.0*[0.1; 0.01; 0.1],...
     "n", @(t,z)[0;0;0]);

intruder = uav(path, dist);

% CONTROLLER
ctrl = lqr_velocity_controller(quad, dt);

% SIMULATION

sim = simulator(quad, ctrl, intruder);
sim.simtime = [0 50];
sim.timestep = dt;
sim.epsilon = 0.1;

z0 = zeros(12,1);

fprintf("Running simulation...\n")
[t,z,u,d,y] = sim.simulate(z0);
fprintf("Simulation complete!\n")

%% ANIMATION
sim.animate(t, z, y);





