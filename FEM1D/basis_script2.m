clear
clc
syms x h;

% Order of polynomial 
N = 3;

a = 0;  b = 1;
h1 = (b-a)/N;

D = 0:h1:1;

D = h*D;
A = sym(zeros(N+1,1));
for i=1:size(D,2)
    poly = 1;
    for j=1:size(D,2)
        if(i ~= j)
            poly = poly*(x-D(j))/(D(i)-D(j));
        end
    end
    
    A(i) = simplify(poly);

    %ezplot(poly,[a,b]);
    %axis([a,h*b,a,h*b+0.1]);
end

K = sym(zeros(N+1,N+1));
for i = 1:N+1
    for j = 1:N+1
        K(i,j) = int(diff(A(i))*diff(A(j)),x,0,h) + int(A(i)*A(j),x,0,h);
    end
end
K