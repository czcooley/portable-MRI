function os = itermsg(itermeth,tol,~,i,flag,iter,relres)
%ITERMSG   Displays the final message for iterative methods.
%   ITERMSG(ITERMETH,TOL,MAXIT,I,FLAG,ITER,RELRES)
%
%   See also BICG, BICGSTAB, BICGSTABL, CGS, GMRES, LSQR, MINRES, PCG, QMR,
%   SYMMLQ, TFQMR.

%   Copyright 1984-2013 The MathWorks, Inc. 

if flag == 0
    if iter == 0
        if isnan(relres)
            os = getString(message('MATLAB:itermsg:RHSVectorAllZero', itermeth));
        else
            os = getString(message('MATLAB:itermsg:InitialGuessHasRelativeWithinTol', ...
                sprintf('%0.2g',relres), sprintf('%0.2g',tol), itermeth));
        end
    else
        os = getString(message('MATLAB:itermsg:ConvergedWithRelativeResidual', itermeth, getIterationInfo(iter, true), ...
            sprintf('%0.2g',relres)));
    end
else
    switch flag
        case 1,
            ncnv =getString(message('MATLAB:itermsg:StoppedMaxIterations', ...
                itermeth, getIterationInfo(i, true), sprintf('%0.2g',tol)));
        case 2,
            ncnv =getString(message('MATLAB:itermsg:StoppedPreconditionerSystemIllCond', ...
                itermeth, getIterationInfo(i, true), sprintf('%0.2g',tol)));
        case 3,
            ncnv =getString(message('MATLAB:itermsg:StoppedMethodStagnated', ...
                itermeth, getIterationInfo(i, true), sprintf('%0.2g',tol)));
        case 4,
            ncnv =getString(message('MATLAB:itermsg:StoppedScalarTooSmallOrLarge', ...
                itermeth, getIterationInfo(i, true), sprintf('%0.2g',tol)));
        case 5,
            ncnv = getString(message('MATLAB:itermsg:StoppedPreconditionerNotSymPosDef', ...
                itermeth, getIterationInfo(i, true), sprintf('%0.2g',tol)));
    end
    retStr = getString(message('MATLAB:itermsg:IterateReturnedHasRelativeResidual', ...
        getIterationInfo(iter, false), sprintf('%0.2g',relres)));
    os = sprintf('%s\n%s', ncnv, retStr);
end
disp(os)

function itstr = getIterationInfo(it, verbose)
if length(it) == 2 % gmres
    if verbose
        itstr = getString(message('MATLAB:itermsg:OuterIteration', ...
            it(1), it(2)));
    else
        itstr = getString(message('MATLAB:itermsg:NumberIntParen', it(1), it(2)));
    end
elseif fix(it) ~= it % bicgstab
    if verbose
        itstr = getString(message('MATLAB:itermsg:IterationFloat', sprintf('%.1f',it)));
    else
        itstr = getString(message('MATLAB:itermsg:NumberFloat', sprintf('%.1f',it)));
    end
else
    if verbose
        itstr = getString(message('MATLAB:itermsg:IterationInteger', it));
    else
        itstr = getString(message('MATLAB:itermsg:NumberInteger', it));
    end
end

