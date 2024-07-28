function [A, B] = linearize(zs, us)
    [A, B] = linearize_symbolic();
    [A, B] = eval_linear(A, B, zs, us);
end