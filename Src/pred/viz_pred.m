

%% Example Path
speed = 0.5;
path = @(t) [cos(t*speed); sin(t*speed); 2];
% path = @(t) [cos(t*speed); (t*speed); 1+sin(t*speed)];
% path = @(t) [-1+ speed*t; 0; 1];

%% Prediction Options
tp = 1; % time prediction occurs at
Nb = 0.1; % number of seconds looking back
Nf = 4; % number of seconds to show looking forwards
dt = 0.01; % timestep for discrete samples

%% Make Prediction
[samples, ts] = genSamples(path, tp, dt, Nb/dt);
polyFun = polyPred(samples, ts);

%% Vizualize
[futureSamples, ~] = genSamples(path, tp+Nf, dt, Nf/dt);
xReal = samples(1,:);
yReal = samples(2,:);
zReal = samples(3,:);
xFuture = futureSamples(1,:);
yFuture = futureSamples(2,:);
zFuture = futureSamples(3,:);


predT = linspace(tp, tp+Nf, Nf/dt);
predictions = polyFun(predT);
xPred = predictions(1,:);
yPred = predictions(2,:);
zPred = predictions(3,:);

plot3(xReal, yReal, zReal, xFuture, yFuture, zFuture, xPred, yPred, zPred, '--');
axis([-2 2 -2 2 0 4]);



%% Functions
function [samples, ts] = genSamples(path, tp, dt, Ns)
% Ns is number of SAMPLES
    samples = zeros([3, Ns]);
    ts = linspace(tp - Ns*dt, tp, Ns);
    i = 1;
    for t = ts
        samples(:, i) = path(t);
        i = i+1;
    end
end


