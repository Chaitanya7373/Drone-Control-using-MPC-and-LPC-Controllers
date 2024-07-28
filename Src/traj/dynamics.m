function zdot = dynamics(z,u)
    % Known Quantities
    g = 9.81;  % The gravitational acceleration [m/s^2]
    l = 0.2;  % Distance from the center of mass to each rotor [m]
    m = 0.5;  % Total mass of the quadrotor [kg]
    I = diag([1.24, 1.24, 2.48]);  % Mass moment of inertia [kg m^2]
    sigma = 0.01;  % The proportionality constant relating thrust to torque [m]

    v = z(4:6);
    a = z(7:9);
    w = z(10:12);

    a1 = a(1);
    a2 = a(2);
    a3 = a(3);

    u1 = u(1);
    u2 = u(2);
    u3 = u(3);
    u4 = u(4);

    % Define coordinate
    e3 = [0 0 1]';

    c1 = [1 0 0]';
    c2 = [0 1 0]';
    c3 = [0 0 1]';
    
    % No external forces/moments for now
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
    R_CE = R3*R2*R1;
    
    T_inv = [1 sin(a1)*tan(a2) cos(a1)*tan(a2);
             0 cos(a1) -sin(a1);
             0 sin(a1)/cos(a2) cos(a1)/cos(a2)];

    xdot = v;
    adot = T_inv * w;
    vdot = -g*e3 + 1/m*R_CE*(u1+u2+u3+u4)*c3 + (1/m)*R_CE*r;
    I_inv = I^(-1);
    wdot = (I_inv)*((u2-u4)*l*c1+(u3-u1)*l*c2 + (u1-u2+u3-u4)*sigma*c3+n-cross(w, I*w));

    zdot = [xdot; adot; vdot; wdot;];
end