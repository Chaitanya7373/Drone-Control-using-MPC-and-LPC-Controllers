% lets see if we can get direct collocation to work...
clear;
clc;
%% Example path
speed = 0.5;
path = @(t) [speed*t; 0; 0.1];

quadInit = zeros([16,1]);

%% Make vector of decision variables
T = 1;
N = 10;
dt = T/N;
ts = linspace(0, T, N);

pathSamples = cell2mat(arrayfun(path, ts, UniformOutput=false))';
uavEnd = pathSamples(end, :);

xs = linspace(quadInit(1),uavEnd(1), N);
ys = linspace(quadInit(2),uavEnd(2), N);
zs = linspace(quadInit(3),uavEnd(3), N);


% 16*N decision variables
x0 = [xs, ys, zs, zeros([1,N]), zeros([1,N]), zeros([1,N]), ...
      zeros([1,N]), zeros([1,N]), zeros([1,N]), ...
      zeros([1,N]), zeros([1,N]), zeros([1,N]), ...
      zeros([1,N]), zeros([1,N]), zeros([1,N]), zeros([1,N])]; % inputs

%% Constraints
thrust_upper = 3;
thrust_lower = 0;
A = [];
b = [];
Aeq = [];
beq = [];
lb = [-Inf*ones([1,N]), -Inf*ones([1,N]), zeros([1,N]), ...
      -Inf*ones([1,N]), -Inf*ones([1,N]), -Inf*ones([1,N]), ...
      -Inf*ones([1,N]), -Inf*ones([1,N]), -Inf*ones([1,N]), ...
      -Inf*ones([1,N]), -Inf*ones([1,N]), -Inf*ones([1,N]), ...
      thrust_lower*ones([1,N]), thrust_lower*ones([1,N]), ...
      thrust_lower*ones([1,N]), thrust_lower*ones([1,N])];
ub = [Inf*ones([1,N]), Inf*ones([1,N]), Inf*ones([1,N]), ...
      Inf*ones([1,N]), Inf*ones([1,N]), Inf*ones([1,N]), ...
      Inf*ones([1,N]), Inf*ones([1,N]), Inf*ones([1,N]), ...
      Inf*ones([1,N]), Inf*ones([1,N]), Inf*ones([1,N]), ...
      thrust_upper*ones([1,N]), thrust_upper*ones([1,N]), ...
      thrust_upper*ones([1,N]), thrust_upper*ones([1,N])];

%% Running fmincon
opts = optimoptions('fmincon');
opts.MaxFunctionEvaluations = 3000;
fprintf("Running fmincon...\n");
[x,fval] = fmincon(@(x)cost(x,N,pathSamples),x0,A,b,Aeq,beq,lb,ub,...
    @(x)constraints(x,N,ts,quadInit), opts);

plot3(pathSamples(:,1), pathSamples(:,2), pathSamples(:,3), ...
      x(1:N), x(N+1:2*N), x(2*N+1:3*N));




