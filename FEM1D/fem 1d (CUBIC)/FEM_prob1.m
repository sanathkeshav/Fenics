%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            -u" = 1   0 < x < 1                         %%%%
%%%%            u(0) = 0   u(1) = 0                         %%%%
%%%%         Exact solution u = x*(1-x)/2                    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

format long;

N = 100;
x0 = 0;
xN = 1;
h = (xN-x0)/N;

for j = 1:3*N+1
    X(j) = x0 + (j-1)*h/3;
end

A = zeros(3*N+1,3*N+1);

F = zeros(3*N+1,1);

%a = [1/h -1/h ; -1/h 1/h];

a = [   37/(10*h), -189/(40*h),   27/(20*h),  -13/(40*h);
 -189/(40*h),    54/(5*h), -297/(40*h),   27/(20*h);
  27/(20*h), -297/(40*h),    54/(5*h), -189/(40*h);
 -13/(40*h),   27/(20*h), -189/(40*h),   37/(10*h)];        %%%% Local stiffness matrix %%%%


for i=1:3:3*N -2
    
    phi1 = @(x)(9*(X(i+3) - x)*(X(i+1) - x)*(X(i+2) - x))/(2*h^3);
                                                            %%%% cubic basis function %%%%
    phi2 = @(x)(27*(x - X(i))*(X(i+3) - x)*(X(i+2) - x))/(2*h^3);
    phi3 = @(x)-(27*(x-X(i))*(X(i+3) - x)*(X(i+1) - x))/(2*h^3);
    phi4 = @(x)(9*(x-X(i))*(X(i+1) - x)*((X(i+2) - x)))/(2*h^3);                                                %%%% cubic basis function %%%%
    f1 = @(x)phi1(x);                    %%%% integrand for load vector %%%%
    f2 = @(x)phi2(x);
    f3 = @(x)phi3(x);
    f4 = @(x)phi4(x);%%%% integrand for load vector %%%%
    v(1,i) = gauss(f1,X(i),X(i+3),3);    %%%% element wise values of
    v(2,i) = gauss(f2,X(i),X(i+3),3);    %%%% load vector
    v(3,i) = gauss(f3,X(i),X(i+3),3);
    v(4,i) = gauss(f4,X(i),X(i+3),3);
end

%%%% Assembling %%%%

for i=1:3:3*N-2
    A([i i+1 i+2 i+3],[i i+1 i+2 i+3]) = A([i i+1 i+2 i+3],[i i+1 i+2 i+3]) + a;
    F([i i+1 i+2 i+3],1) = F([i i+1 i+2 i+3],1) + v([1 2 3 4],i);
end

fullnodes = [1:3*N+1];
%%%%% Dirichlet boundary condition %%%%%

freenodes=setdiff(fullnodes,[1,3*N+1]);

Uh = zeros(3*N+1,1);

%%%% Approximate solution %%%%

Uh(freenodes)=A(freenodes,freenodes)\F(freenodes,1);

%%%% Exact solution %%%%

U = zeros(3*N+1,1);

for i =1:3*N+1
    U(i) = X(i)*(1-X(i))/2;
end
error = norm(U - Uh,Inf);

plot(X,U, X,Uh,'o')