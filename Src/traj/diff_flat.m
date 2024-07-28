%% Example Path
start = [0;0;0];
final = [1;1;1];
time = 1;

g = 9.81;
% parameterized straight line path with constant yaw
syms t
pos_path = [(final-start)/time * t; 0];
vel_path = simplify(jacobian()


