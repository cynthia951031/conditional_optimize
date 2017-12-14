function s = target(n, x)
    s = 0;
    for i = 1: 29
        s = s + r(i, n, x)^2;
    end
    s = s + (x(1))^2 + (x(2) - x(1)^2 - 1)^2;
end