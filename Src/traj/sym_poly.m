function sym_fun = sym_poly(p)
syms t;
sym_fun = p(1) + p(2)*t + p(3)*t^2 + p(4)*t^3 + p(5)*t^4 + p(6)*t^5;
end