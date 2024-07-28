function c = cost(x, N, target)
    x_pos = x(1:N);
    y_pos = x(N+1:2*N);
    z_pos = x(2*N+1:3*N);

    u = x(12*N+1:end);

    x_targ = target(:, 1)';
    y_targ = target(:, 2)';
    z_targ = target(:, 3)';

    x_err = x_pos - x_targ;
    y_err = y_pos - y_targ;
    z_err = z_pos - z_targ;

    c = x_err*x_err' + y_err*y_err' + z_err*z_err' + u*u';
end