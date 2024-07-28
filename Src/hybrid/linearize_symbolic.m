function [A, B] = linearize_symbolic()
    
    % To take jacobion as function of these
    syms x1 x2 x3 a1 a2 a3 v1 v2 v3 w1 w2 w3
    syms u1 u2 u3 u4

    x = [x1 x2 x3]';
    a = [a1 a2 a3]';
    v = [v1 v2 v3]';
    w = [w1 w2 w3]';

    z = [x; a; v; w;];
    u = [u1 u2 u3 u4]';

    % Known Quantities
    g = 9.81;  % The gravitational acceleration [m/s^2]
    l = 0.2;  % Distance from the center of mass to each rotor [m]
    m = 0.5;  % Total mass of the quadrotor [kg]
    I = diag([1.24, 1.24, 2.48]);  % Mass moment of inertia [kg m^2]
    sigma = 0.01;  % The proportionality constant relating thrust to torque [m]


    
    % Define coordinate, correct??
    e3 = [0 0 1]';

    c1 = [1 0 0]';
    c2 = [0 1 0]';
    c3 = [0 0 1]';
    
    % No external forces/moments for now
    %syms r1 r2 r3
    %r = r1*c1+r2*c2+r3*c3;
    %syms n1 n2 n3
    %n = n1*c1+n2*c2+n3*c3;
    r = zeros([3,1]);
    n = zeros([3,1]);
    
    % Euler Angle Rotation Matrices
    R3 = [cos(a3) -sin(a3) 0;
          sin(a3) cos(a3)  0;
          0 0 1];
    R2 = [cos(a2) 0 sin(a2);
          0 1 0;
          -sin(a2) 0 cos(a2)];
    R1 = [1 0 0;
          0 cos(a1) -sin(a1);
          0 sin(a1) cos(a1)];
    R_CE = simplify(R3*R2*R1);
    
    T_inv = [1 sin(a1)*tan(a2) cos(a1)*tan(a2);
             0 cos(a1) -sin(a1);
             0 sin(a1)/cos(a2) cos(a1)/cos(a2)];

    xdot = v;
    adot = T_inv * w;
    vdot = -g*e3 + 1/m*R_CE*(u1+u2+u3+u4)*c3 + (1/m)*R_CE*r;
    I_inv = I^(-1);
    wdot = (I_inv)*((u2-u4)*l*c1+(u3-u1)*l*c2 + (u1-u2+u3-u4)*sigma*c3+n-cross(w, I*w));

    f = [xdot; adot; vdot; wdot;];
    
    A = jacobian(f,z');
    B = jacobian(f,u');
end