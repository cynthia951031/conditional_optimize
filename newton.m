function [lambda, count] = newton(f, x, init_val)
    iteration = 100;
    threshold = 1e-6;
    x_new = init_val;
    
    count = 0;
    delta_f = diff(f, x);
    delta2_f = diff(delta_f, x);
    for i = 1: iteration
       g_f = subs(delta_f, x, x_new);
       if abs(g_f) < threshold
          break; 
       end
       g2_f = subs(delta2_f, x, x_new);
       count = count + 2;
       x_new = x_new - g_f / g2_f;
    end
    lambda = x_new;
end