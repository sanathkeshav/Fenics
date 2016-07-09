function I = gauss(f,a,b,n)
%Gaussian quadrature of order 1,2,3
% n = order of gaussian quadrature
clear I;
xi = [0 -1/sqrt(3) -sqrt(3/5) ; 0 1/sqrt(3) 0 ; 0 0 sqrt(3/5)];
w = [2 1 5/9 ; 0 1 8/9 ; 0 0 5/9];
p = (b-a)/2;
q = (b+a)/2;
I = 0;
for i=1:3
    I = I + p*w(i,n)*f(p*xi(i,n)+q);
end
end

