function s = r(i, n, x)
    t = i / 29;
    tmp = 0;
    s = 0;
    for j = 2: n
        s = s + (j - 1) * x(j) * t ^ (j - 2);
    end
    for k = 1: n
       tmp = tmp + x(k) * t ^ (k - 1);
    end
    s = s - tmp ^ 2 - 1;
end


