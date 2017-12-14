function [alpha, count] = mb_wolfeLineSearch(objFunc,objFuncValue,x,dx,dir)
alphaMax     = 1;
alpha        = alphaMax;
fac          = 1/2; 
c_1          = 1e-1;
c_2          = 2e-1;
count = 0;

while objFunc(x+alpha*dir) > objFuncValue + c_1*alpha*dir'*dx && mb_numDiff(objFunc, x+alpha*dir)'*dir < c_2*dx'*dir
    count = count + 1;
    alpha = fac*alpha;
    
    if alpha < 10*eps
        error('Error in Line search - alpha close to working precision');
    end
    
end
end