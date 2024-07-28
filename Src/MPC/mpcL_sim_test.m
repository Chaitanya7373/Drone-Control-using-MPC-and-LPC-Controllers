clc; 
clear; 
close all;

% Add path to quadrotor_sim/src if necessary
addpath('../quadrotor_sim/src');
addpath('../traj');
addpath('../lqr');

% QUADROTOR PARAMETERS
g = 9.81;  % The gravitational acceleration [m/s^2]
l = 0.2;  % Distance from the center of mass to each rotor [m]
m = 0.5;  % Total mass of the quadrotor [kg]
I = [1.24, 1.24, 2.48];  % Mass moment of inertia [kg m^2]
mu = 3.0;  % Maximum thrust of each rotor [N]
sigma = 0.01;  % The proportionality constant relating thrust to torque [m]


% Create quadrotor object
quad = quadrotor(g, l, m, diag(I), mu, sigma);

% INTRUDER PARAMETERS
path = @(t) [-cos(t); sin(t); 2];
path = @(t) [0;0;2];
dist = struct("r", @(t)0.1*[sin(t); sin(2*t); sin(4*t)],...
    "n", @(t)[0.1; 0.01; 0.1]);

intruder = uav(path, dist);

% Create MPC_L controller

controller = mpc_controller(quad);

% SIMULATION PARAMETERS
simtime = [0 10]; % Simulation time [s]
timestep = 0.01; % Timestep [s]
z0 = zeros(12, 1); % Initial state

% SIMULATION
sim = simulator(quad, controller, intruder);
sim.simtime = simtime;
sim.timestep = timestep;

fprintf("Running simulation...\n")
[t, z, u, d, y] = sim.simulate(z0);
fprintf("Simulation complete!\n")

%% ANIMATION (if available)
sim.animate(t, z, y);
quad.plot(t, z);

