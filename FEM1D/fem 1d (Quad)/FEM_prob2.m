%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%         -u" = x*(3+x)*e^x on 0 < x < 1                 %%%%
%%%%            u(0) = 0   u(1) = 0                         %%%%
%%%%         Exact solution u = x(1-x)*e^x                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




clear all;
close all;

format long;
    
N = 10;
x0 = 0;
xN = 1;
    


for p=1:4
    h(p) = (xN-x0)/N;
    for j = 1:2*N+1
        X(j) = x0 + (j-1)*h(p)/2;
    end

    f = @(x)x*(3+x)*exp(x);

    A = zeros(2*N+1,2*N+1);

    F = zeros(2*N+1,1);

    a = 1/(3*h(p))*[7 -8 1; -8 16 -8; 1 -8 7];               %%%% Local stiffness matrix %%%%


    for i=1:2:2*N-1
        phi1 = @(x)(x-X(i+1))*(x-X(i+2))/(h(p)^2/2);    %%%% quadratic basis function %%%%
        phi2 = @(x)(x-X(i))*(x-X(i+2))/(-h(p)^2/4);     %%%% quadratic basis function %%%%
        phi3 = @(x)(x-X(i))*(x-X(i+1))/(h(p)^2/2);      %%%% quadratic basis function %%%%
        f1 = @(x)f(x)*phi1(x);                    %%%% integrand for load vector %%%%
        f2 = @(x)f(x)*phi2(x);                    %%%% integrand for load vector %%%%
        f3 = @(x)f(x)*phi3(x);                    %%%% integrand for load vector %%%%  
        v(1,i) = gauss(f1,X(i),X(i+2),3);    %%%% element wise values of
        v(2,i) = gauss(f2,X(i),X(i+2),3);    %%%% load vector
        v(3,i) = gauss(f3,X(i),X(i+2),3);     %%%% load vector
    end

    %%%% Assembling %%%%

    for i=1:2:2*N-1
        A([i i+1 i+2],[i i+1 i+2]) = A([i i+1 i+2],[i i+1 i+2]) + a;
        F([i i+1 i+2],1) = F([i i+1 i+2],1) + v([1 2 3],i);
    end

    fullnodes = [1:2*N+1];

    %%%%% Dirichlet boundary condition %%%%%

    freenodes=setdiff(fullnodes,[1,2*N+1]);

    Uh = zeros(2*N+1,1);

    %%%% Approximate solution %%%%

    Uh(freenodes)=A(freenodes,freenodes)\F(freenodes,1);

    %%%% Exact solution %%%%

    U = zeros(2*N+1,1);

    for i =1:2*N+1
        U(i) = X(i)*(1-X(i))*exp(X(i));
    end

    error(p) = norm(U-Uh);

    N = N*2;
end

for j=1:p-1
    order(j) = log(error(j)/error(j+1))/log(2);
end
figure(1)
plot(log(h),log(error));
figure(2)
plot(X,U,X,Uh,'o');

order'