%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%       -u" + u = x^2 + x - 2   0 < x <= 1               %%%%
%%%%            u(0) = 0, u(1) + u'(1) = 5                  %%%%
%%%%         Exact solution u = x*(x+1)                     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
clear all;
close all;

format long;

N = 20;

x0 = 0;
xN = 1;

h = 1/N;

for j = 1:N+1
X(j) = x0 + (j-1)*h;
end

f = @(x)(x^2 + x - 2);

A = zeros(N+1,N+1);

F = zeros(N+1,1);

a = [1/h -1/h ; -1/h 1/h] + [h/3 h/6; h/6 h/3]; %%%% Local stiffness matrix %%%%


for i=1:N
    phi1 = @(x)(X(i+1)-x)/h;            %%%% linear basis function %%%%
    phi2 = @(x)(x-X(i))/h;              %%%% linear basis function %%%%
    f1 = @(x)f(x)*phi1(x);              %%%% integrand for load vector %%%%
    f2 = @(x)f(x)*phi2(x);              %%%% integrand for load vector %%%%
    v(1,i) = gauss(f1,X(i),X(i+1),3);   %%%% element wise values of
    v(2,i) = gauss(f2,X(i),X(i+1),3);   %%%% load vector
end

%%%% Assembling %%%%

for i=1:N
    A([i i+1],[i i+1]) = A([i i+1],[i i+1]) + a;
    F([i i+1],1) = F([i i+1],1) + v([1 2],i);
end

%%%% Mixed Boundary condition %%%%

A(N+1,N+1) = A(N+1,N+1) + 1;
F(N+1,1) = F(N+1,1) + 5;

fullnodes = [1:N+1];

%%%%% Dirichlet boundary condition %%%%%

freenodes=setdiff(fullnodes,[1]);

Uh = zeros(N+1,1);

%%%% Approximate solution %%%%

Uh(freenodes)=A(freenodes,freenodes)\F(freenodes,1);

%%%% Exact solution %%%%

U = zeros(N+1,1);

for i =1:N+1
    U(i) = X(i)*(1+X(i));
end

error = max(abs(U-Uh));


[U Uh]

plot(X,U,X,Uh,'o')
