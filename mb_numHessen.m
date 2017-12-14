function [hp, count] = mb_numHessen(deltaFunc, x, p, h)
    n = numel(p);
    if nargin < 4
        h  = 1e-8;
    end
    oldDiffValue = subs(deltaFunc, x, p);
    hp = zeros(n, n);
    for i = 1:n
       for j = 1:n
          p_new = p;
          p_new(j) = p_new(j) + h;
          
          newDiffValue = subs(deltaFunc, x, p_new);
          hp(i, j) = (newDiffValue(i) - oldDiffValue(i)) / h;
          
       end
    end
    count = n * n;
end