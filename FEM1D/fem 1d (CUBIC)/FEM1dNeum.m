%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%       -u" + u = x^3 + 2x^2 - 5x - 4   0 < x <= 1       %%%%
%%%%            u(0) = 0, u'(1) = 8                         %%%%
%%%%         Exact solution u = x*(x+1)^2                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

format long;
    
N = 5;
x0 = 0;
xN = 1;
    
for p = 1:4
h(p) = 1/N;

for j = 1:3*N+1
X(j) = x0 + (j-1)*h(p)/3;
end

f = @(x)x^3 + 2*x^2 - 5*x - 4;

A = zeros(3*N+1,3*N+1);

F = zeros(3*N+1,1);

a = [   (8*h(p))/105 + 37/(10*h(p)),   (33*h(p))/560 - 189/(40*h(p)),     27/(20*h(p)) - (3*h(p))/140, (19*h(p))/1680 - 13/(40*h(p));
 (33*h(p))/560 - 189/(40*h(p)),      (27*h(p))/70 + 54/(5*h(p)), - (27*h(p))/560 - 297/(40*h(p)),   27/(20*h(p)) - (3*h(p))/140;
   27/(20*h(p)) - (3*h(p))/140, - (27*h(p))/560 - 297/(40*h(p)),      (27*h(p))/70 + 54/(5*h(p)), (33*h(p))/560 - 189/(40*h(p));
 (19*h(p))/1680 - 13/(40*h(p)),     27/(20*h(p)) - (3*h(p))/140,   (33*h(p))/560 - 189/(40*h(p)),   (8*h(p))/105 + 37/(10*h(p))];        %%%% Local stiffness matrix %%%%


for i=1:3:3*N-2
    
    phi1 = @(x)(9*(X(i+3) - x)*(X(i+1) - x)*(X(i+2) - x))/(2*h(p)^3);
                                                            %%%% cubic basis function %%%%
    phi2 = @(x)(27*(x - X(i))*(X(i+3) - x)*(X(i+2) - x))/(2*h(p)^3);
    phi3 = @(x)-(27*(x-X(i))*(X(i+3) - x)*(X(i+1) - x))/(2*h(p)^3);
    phi4 = @(x)(9*(x-X(i))*(X(i+1) - x)*((X(i+2) - x)))/(2*h(p)^3);                                                %%%% cubic basis function %%%%
    f1 = @(x)f(x)*phi1(x);                    %%%% integrand for load vector %%%%
    f2 = @(x)f(x)*phi2(x);
    f3 = @(x)f(x)*phi3(x);
    f4 = @(x)f(x)*phi4(x);%%%% integrand for load vector %%%%
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

%%%% Neumann Boundary condition %%%%

F(3*N+1,1) = F(3*N+1,1) + 8;

fullnodes = [1:3*N+1];

%%%%% Dirichlet boundary condition %%%%%

freenodes=setdiff(fullnodes,[1]);

Uh = zeros(3*N+1,1);

%%%% Approximate solution %%%%

Uh(freenodes)=A(freenodes,freenodes)\F(freenodes,1);

%%%% Exact solution %%%%

U = zeros(3*N+1,1);

for i =1:3*N+1
    U(i) = X(i)*(1+X(i))^2;
end

error(p) = max(abs(U-Uh));

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

