function [di, iCon, resF] = MNRwithBFGSlineSearchable(Fi, d0, x, LSflag)
%Modified Newton with BFGS updates with/without line search

Fi = [Fi ; 0];
d0 = [d0 ; d0];

tol = 10^(-4); %tolerance
maxiter = 15; %maximum number of iterations

Ktilde = tangentMatrix(d0, x);

if LSflag == false %no line search

    dd = Ktilde \ residual(Fi, d0, x); %initial solve
    diOld = d0;
    si = 1;
    diNew = d0 + si * dd;

    Rnew =  residual( Fi, diNew, x);
    Rold = residual(Fi, diOld, x);

    for i = 2:maxiter
        
        Gi = @(s) dd.' * residual(Fi, diOld + s *dd, x);
        
        vk = dd / (Gi(si) - Gi(0));
        alpha_k = sqrt( - si* (Gi(si) - Gi(0))/ Gi(0)  );
        wk = -(Rnew - Rold) + alpha_k * Rnew;

        a = Rnew + (vk.' * Rnew)*wk ;
        b = Ktilde \ a;
        
        dd = b + (wk.' * b)*vk ;

        diOld = diNew;
        si = 1;
        diNew = diOld + si * dd;
        
        Rnew = residual(Fi, diNew, x);
        Rold = residual(Fi, diOld, x);

        if  norm(Rnew)  < tol *norm(residual(Fi, d0,x))
            di = diNew(1);
            iCon = i;
            %disp(['Converged at iteration ', num2str(i)])
            break
        end
    end

    if i == maxiter
        %disp(['Did not converge. Last value is  ', num2str(diOld(1))])
        di = 0;
        iCon = i;
    end
    resF = residual(Fi, diNew, x);
    resF = resF(1);

else %with line search

    dd = Ktilde \ residual(Fi, d0, x); %initial solve
    diOld = d0;
    si = lineSearch(dd, diOld, x);
    diNew = d0 + si * dd;

    Rnew =  residual( Fi, diNew, x);
    Rold = residual(Fi, diOld, x);

    for i = 2:maxiter
        
        Gi = @(s) dd.' * residual(Fi, diOld + s *dd, x);
        
        vk = dd / (Gi(si) - Gi(0));
        alpha_k = sqrt( - si* (Gi(si) - Gi(0))/ Gi(0)  );
        wk = -(Rnew - Rold) + alpha_k * Rnew;

        a = Rnew + (vk.' * Rnew)*wk ;
        b = Ktilde \ a;
        
        dd = b + (wk.' * b) * vk ;

        diOld = diNew;
        si = lineSearch(dd, diOld, x);
        diNew = diOld + si * dd;
        
        Rnew = residual(Fi, diNew, x);
        Rold = residual(Fi, diOld, x);

        if  norm(Rnew)  < tol * norm(residual(Fi, d0,x))
            di = diNew(1);
            iCon = i;
            %disp(['Converged at iteration ', num2str(i)])
            break
        end
    end

    if i == maxiter
        disp(['Did not converge. Last value is  ', num2str(diOld(1))])
        di = diNew(1);
    end
    resF = residual(Fi, diNew, x);
    resF = resF(1);




end








    function FminN = residual(Fterm, d, xt)
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
    end


end