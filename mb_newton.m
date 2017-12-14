function mb_newton(searchMode, mode, n)

fprintf(['\n N = ' num2str(n) '\n' ...
           'Using ' mode ' method update\n'...
           'Using ' searchMode ' method update\n']);
epslon = 20;

x = zeros(n, 1);
X = sym('x', [n, 1]);
func = target(n, X);
dFunc = jacobian(func, X);
objFunc = @(x) target(n, x);
objFuncValue = objFunc(x);
callCounter = 1;
oldObjFuncValue = objFuncValue + 1;
[dx, tempCounter] = mb_numDiff(objFunc, x);
diffCounter = tempCounter;
[G, tempCounter] = mb_numHessen(dFunc, X, x);
diffCounter = diffCounter + tempCounter;

iter = 0;
numOfIter = 100;
prec = 1e-6;

while iter < numOfIter && norm(dx) > prec && abs((oldObjFuncValue-objFuncValue)/objFuncValue)>prec
   iter = iter + 1;
   oldObjFuncValue = objFuncValue;
   if strcmp(mode, 'DUMP')
       dir = -(G\dx);
   elseif strcmp(mode, 'MODIFIED')
       dir = -((G + epslon * eye(n)) \ dx);
   end
   
   if strcmp(searchMode, 'ARMIJO')
        [lambda, tempCounter] = mb_armijoLineSearch(objFunc,objFuncValue,x,dx,dir);
   elseif strcmp(searchMode, 'WOLFE')
        [lambda, tempCounter] = mb_wolfeLineSearch(objFunc,objFuncValue,x,dx,dir);
   end
   callCounter = callCounter + tempCounter;
   p = lambda * dir;
   x = x + p;
   objFuncValue = objFunc(x);
   callCounter = callCounter + 1;
   [dx, tempCounter] = mb_numDiff(objFunc, x);
   diffCounter = diffCounter + tempCounter;
   [G, tempCounter] = mb_numHessen(dFunc, X, x);
   diffCounter = diffCounter + tempCounter;
   
   fprintf(1,'Iteration %d: alpha_min=%f, OF=%f\n',iter,lambda,objFuncValue);
   
end
fprintf(['\n' num2str(iter) ' iteration(s) performed to converge\n'])
fprintf(['diffCounter: ' num2str(diffCounter) '\ncallCounter: ' num2str(callCounter) '\n'])
fprintf('Final solution: \n');
display(x);
display(objFunc(x))
end