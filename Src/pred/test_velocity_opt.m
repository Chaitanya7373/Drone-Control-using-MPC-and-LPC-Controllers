%% Problem Setup
dt = 0.01;

y = [-1; 1; 0];
x = [-0; 0; 0];

vhat = [-1; 0; 0];
s = 2;

%% Gradient Descent

% start with vd = vhat as initial condition
vd = vhat;

% do updates
n = 1000;
lr = 1;
disp("Running gradient descent for velocity optimization...")
for i = 1:n
    e = e_next(x, y, vhat, vd, s, dt);
    J = -dt*e/norm(e); % +vd for l2 regularization
    vd = vd - lr*J;
end
disp("...done")
disp("Result: vd = ")

disp(vd);
disp(norm(vd));
disp(vd/norm(vd) * s)

%% Function for error
function en = e_next(x, y, vhat, vd, s, dt)
    e = y - x;
    eta = norm(e) / s
    eta = 3
%     eta = dt;
    en = (y + vhat*eta) - (x + vd*eta);
end

