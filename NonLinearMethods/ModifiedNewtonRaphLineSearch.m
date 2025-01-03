function [di, iCon, resF] = ModifiedNewtonRaphLineSearch(Fi, d0, x)
%solve 2d problem using modified newton-raphson with line search.

Fi = [Fi ; 0];
d0 = [d0 ; d0];

tol = 10^(-4); %tolerance
maxiter = 15; %maximum number of iterations

Ktilde = tangentMatrix(d0, x);
dd = Ktilde\ residual(Fi, d0, x); %initial solve
s = lineSearch(dd, d0, x);
diOld = d0 + s * dd;

for i=2:maxiter
dd = Ktilde \ residual(Fi, diOld, x) ;
s = lineSearch(dd, diOld, x);
diNew = diOld + s * dd;

if  norm(residual(Fi, diNew, x))  < tol *norm(residual(Fi, d0,x))
    di = diNew(1);
    iCon = i;
    %disp(['Converged at iteration ', num2str(i)])
    break

else
    diOld = diNew;
end

end

resF = residual(Fi, diNew, x);
resF = resF(1);

if i == maxiter
disp(['Did not converge. Last value is %f ', num2str(diOld(1))])
di = real(diNew(1));
iCon = i;
end




    function FminN = residual(Fterm, d,xt)
        if nargin ==2
            xt=x;
        end
        N = [(xt* d(1))/(10 - d(1)) - 0.5 * d(2)^2; d(2) - d(1)];

        FminN = Fterm - N; 
        
    end

    function DN = tangentMatrix(d, xt)
        if nargin ==1
            xt = x;
        end

        DN = [10 *xt /(10 - d(1))^2 , -d(2); -1 , 1];

    end

    function si = lineSearch(deltad, Old_d, xt)
        if nargin==2
           xt=x; 
        end

        G = @(s) deltad.' * residual(Fi, Old_d + s*deltad, xt);
        dGds = @(s) deltad.' * tangentMatrix(Old_d + s*deltad, xt) * deltad;

        lineTol = 0.5*G(0); %line search tolerance
        lineMaxiter = 5;  %line search max iterations

        ds = dGds(0) \ G(0);
        sOld = ds;

        for j = 2:lineMaxiter
            ds = dGds(sOld) \ G(sOld);
            sNew = sOld + ds;

            if norm(G(sNew)) < lineTol
                si = sNew;
                %disp(['LineSearch converged at iteration ', num2str(j)])
                break
            else
                sOld = sNew;
            end
        end

        if j == lineMaxiter
        %disp(['Line search did not converge. Last value is %f ', ...
            %num2str(sOld)])
        si = sNew;
        end

    end




end