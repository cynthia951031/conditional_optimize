function mb_quasiNewton(searchMode, mode, n)
fprintf(['N = ' num2str(n) '\n' ...
           'Using ' mode ' method update\n'...
           'Using ' searchMode ' method update\n']);

x       = zeros(n, 1); 
H       = eye(n);  

objFunc         = @(x) target(n, x);
objFuncValue    = objFunc(x);
callCounter = 1;
oldObjFuncValue = objFuncValue + 1;
[dx, tempCounter] = mb_numDiff(objFunc,x);
diffCounter = tempCounter;

iter      = 0;
numOfIter = 100;
prec      = 1e-6;

while iter < numOfIter && abs((oldObjFuncValue-objFuncValue)/objFuncValue)>prec && norm(dx)>prec
    iter = iter + 1;
    oldObjFuncValue = objFuncValue;
    dir = -(H\dx); 
    if strcmp(searchMode, 'ARMIJO')
        [alpha, tempCounter] = mb_armijoLineSearch(objFunc,objFuncValue,x,dx,dir);
    elseif strcmp(searchMode, 'WOLFE')
        [alpha, tempCounter] = mb_wolfeLineSearch(objFunc,objFuncValue,x,dx,dir);
    end
    callCounter = callCounter + tempCounter;
    p = alpha*dir;
    x = x + p;
    objFuncValue = objFunc(x);
    callCounter = callCounter + 1;
    dx_old = dx;
    [dx, tempCounter] = mb_numDiff(objFunc,x);
    diffCounter = diffCounter + tempCounter;
    q = dx-dx_old;

    if strcmp(mode,'DFP')
        H = H + (q*q')/(q'*p) - ((H*p)*(H*p)')/(p'*H*p);
    elseif strcmp(mode,'SR1')
        H = H + ((q-H*p)*(q-H*p)')/((q-H*p)'*p);
    elseif strcmp(mode, 'BFGS')
        H = H + (1+(p'*H*p)/(p'*q))*(q*q')/(p'*q) - (q*p'*H+H*p*q')/(p'*q);
    end
    
    fprintf(1,'Iteration %d: alpha_min=%f, OF=%f\n',iter,alpha,objFuncValue);
    
end
fprintf(['\n' num2str(iter) ' iteration(s) performed to converge\n'])
fprintf(['diffCounter: ' num2str(diffCounter) '\n callCounter: ' num2str(callCounter) '\n'])
fprintf(1,'Final solution: \n');
display(x);
display(objFunc(x));
end