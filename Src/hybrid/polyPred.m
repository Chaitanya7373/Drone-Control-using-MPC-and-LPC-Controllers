function traj_pred = polyPred(samples, ts, n)
    xP = polyfit(ts, samples(1, :), n);
    yP = polyfit(ts, samples(2, :), n);
    zP = polyfit(ts, samples(3, :), n);
    traj_pred = @(t) [polyval(xP, t); polyval(yP, t); polyval(zP, t)];
end