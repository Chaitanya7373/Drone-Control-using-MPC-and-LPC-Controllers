function [nEst, rEst] = estimateDisturbance(zActual, zPred, quad, dt)
    z = zActual;
    % need rotation to calculate r
    a1 = z(4);
    a2 = z(5);
    a3 = z(6);
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

    % actually look at the error 
    predErr = z-zPred;
    errWdot = predErr(10:12) / dt;
    errVdot = predErr(7:9) / dt;
    nEst = quad.I*errWdot;
    rEst = R_CE'*quad.m*errVdot;
end