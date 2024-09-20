function di = ExactN(x, Fi)
%Exact solution to nonlinear problem

p = [0.5, -5, (x+Fi), -10*Fi];

root = roots(p);

for i =1:length(root)
    if imag(root(i))== 0 %cubic will have exactly 1 real root
        di = (root(i)); 
    end
end



end
