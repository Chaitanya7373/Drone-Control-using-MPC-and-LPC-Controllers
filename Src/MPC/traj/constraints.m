function [c,ceq] = constraints(x, N, ts, posInit)
    c = [];
    % Start with inital condition constraints
    ceq = zeros([12*N,1]);

    x = reshape(x, [16, N]);

    z = x(1:12, :);
    u = x(13:16, :);

    % Enforce z0 constraints
    ceq(1:12) = posInit(1:12) - z(:,1);

    for i = 1:N-1
        z0 = z(:, i);
        z1 = z(:, i+1);
        u0 = u(:, i);
        u1 = u(:, i+1);
%         tc = 0.5*(ts(i) + ts(i+1));
        h = ts(i+1) - ts(i);

        zdot0 = dynamics(z0, u0);
        zdot1 = dynamics(z1, u1);

        uc = 0.5*(u0+u1);
        zc = 0.5*(z0+z1) + h/8*(zdot0 - zdot1);
        zdotc = -3/(2*h)*(z0-z1) - 0.25*(zdot0 + zdot1);

        ceq(13+(i-1)*N:(i-1)*N+24) = zdotc - dynamics(zc, uc);
    end
end