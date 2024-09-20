function [di, iCon, resF] = NewtonRaphNoLineSearch(Fi, d0, x)
%solve 2d problem using newton-raphson without line search.

Fi = [Fi ; 0];
d0 = [d0 ; d0];

tol = 10^(-4); %tolerance
maxiter = 15; %maximum number of iterations

dd = tangentMatrix(d0, x)\ residual(Fi, d0, x); %initial solve
diOld = d0 + dd;

for i=2:maxiter
dd = tangentMatrix(diOld, x) \ residual(Fi, diOld, x) ;
diNew = diOld + dd;

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
di = diNew(1);
end




    function FminN = residual(Fterm, d,xt)
        if nargin ==1
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






end