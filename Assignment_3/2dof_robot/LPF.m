function y = LPF(x, pre_y, ts, tau)
    y = (pre_y.*tau + x.*ts) /(tau + ts);
end