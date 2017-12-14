function [dp, count] = mb_numDiff(func,p,h)

count = 0;
if nargin < 3
    h  = 1e-8;
end

dp = NaN*p;

oldObjFuncValue = func(p);

for i = 1:numel(p)
    
    p_new    = p;
    p_new(i) = p_new(i) + h;
    
    newObjFuncValue = func(p_new);
    
    dp(i) = (newObjFuncValue-oldObjFuncValue)/h;
    count = count + 1;

end

end