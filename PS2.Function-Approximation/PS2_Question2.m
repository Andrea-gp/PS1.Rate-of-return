%Question 2: unknown are capital and labor and parameters are sigma and alpha
alpha=0.5;sigma=0.25;
syms x y 
f=((1-alpha)*x.^((sigma-1)/sigma)+alpha*y.^((sigma-1)/sigma)).^(sigma/(sigma-1));
x=0:0.1:10
y=0:0.1:10
[X Y]=meshgrid(x, y)
Z=((1-alpha).*X.^((sigma-1)/sigma)+alpha.*Y.^((sigma-1)/sigma)).^(sigma/(sigma-1));
surf(X, Y, Z)
title('Evaluation of the CES function with sigma=0.25')
figure;
contourf(X, Y, Z);
title('Isoquants graph')
%Firstly, we are computing the labor share:
alpha=0.5;sigma=0.25;m=20
k=[1:10]
for i=1:10
    LS(i)=1/(((1-alpha)/alpha)*k(i)^((sigma-1)/sigma)+1)
end
plot(k,LS)
title('Labor share on the economy')
xlabel('Capital per worker')

%Steps to create Chebyshev regression algorithm
alpha=0.5;sigma=0.25;m=20 
%STEP 1: we create the Chebyshev nodes;
for k=1:20
zk(k,1)=cos((2*k-1)/(2*m)*pi);
end
%STEP 2: we adjust the nodes to the interval [a,b]=[0,10] for both
%variables
a=0;b=10;int=[a,b];
xk=(zk+1)*((int(2)-int(1))/2)+a;
%Given that both variables are defined in the same interval, we just plug
%in the adjusted nodes in the formula as follows: k = h = xk.
%STEP 3: we evaluate the CES function in the adjusted nodes:
wk=((1-alpha)*xk.^((sigma-1)/sigma)+alpha*xk.^((sigma-1)/sigma)).^(sigma/(sigma-1));
ces_eval_zk=((1-alpha)*xk.^((sigma-1)/sigma)+alpha*xk.^((sigma-1)/sigma)).^(sigma/(sigma-1));
%STEP 4: we compute the coefficient associated with Chebyshev basis. For
%this step, we have created two functions to get the basis for the
%Chebyshev polynomial and using this we have created this polynomial with
%two variables x and y. The order of the polynomial is given by n:
%Example with a polynomial of degree 11:
n=11;space_x=[0:10];space_y=[0:10]
basis_11_x=chebyshev_basis_x(n,space_x);
basis_11_y=chebyshev_basis_y(n,space_y);
basis_11new=flip(basis_11_y)
for i=1:size(basis_11_y,2)
    final_vector(i)=basis_11new(1,i)*basis_11_x(1,i)
end
syms x y 
final=zeros(size(zk,1),n)
for j=1:n
for i=1:size(zk,1)
final1(i)=subs(final_vector(1,j),[x,y],[zk(i),zk(i)])
final(i,j)=final1(i)
end
end
%We compute the Chebyshev coefficients by using OLS:
thetas=(final'*final)\(final'*wk)

%STEP 5: we approximate the CES function
for i=1:n
approx(i)=thetas(i)*subs(final_vector(1,j),[x,y],[(2*(x-a)/(b-a))-1,(2*(y-a)/(b-a))-1])
end
cheb_approx_11=sum(approx)

approx_final=zeros(n,1)
h=linspace(0,10,n);k=linspace(0,10,n);
for j=1:n
approx_final_11(j,1)=double(subs(cheb_approx_11,[x,y],[h(j),k(j)]))
end
for i=1:n
ces_eval_11(i)=((1-alpha)*k(i).^((sigma-1)/sigma)+alpha*h(i).^((sigma-1)/sigma)).^(sigma/(sigma-1));
end
e_11=abs(ces_eval_11'-approx_final_11)
figure(1)
plot(k,approx_final_11)
hold on
plot(k,ces_eval_11)
title('Chebyshev polynomial approximation of order 11')
legend('Approximation','Original function')
figure(2)
plot(k,e_11)
title('Error for Chebyshev polynomial order 11')
ylabel('Error')

%Example with a polynomial of order 15:
n=15;space_x=[0:10];space_y=[0:10]
basis_15_x=chebyshev_basis_x(n,space_x);
basis_15_y=chebyshev_basis_y(n,space_y);
basis_15new=flip(basis_15_y)
for i=1:size(basis_15_y,2)
    final_vector(i)=basis_15new(1,i)*basis_15_x(1,i)
end

syms x y 
final=zeros(size(zk,1),n)
for j=1:n
for i=1:size(zk,1)
final1(i)=subs(final_vector(1,j),[x,y],[zk(i),zk(i)])
final(i,j)=final1(i)
end
end
thetas=(final'*final)\(final'*wk)

for i=1:n
approx(i)=thetas(i)*subs(final_vector(1,j),[x,y],[(2*(x-a)/(b-a))-1,(2*(y-a)/(b-a))-1])
end
cheb_approx_15=sum(approx)
approx_final_15=zeros(n,1)
h=linspace(0,10,n);k=linspace(0,10,n);
for j=1:n
approx_final_15(j,1)=double(subs(cheb_approx_15,[x,y],[h(j),k(j)]))
end
for i=1:n
ces_eval_15(i)=((1-alpha)*k(i).^((sigma-1)/sigma)+alpha*h(i).^((sigma-1)/sigma)).^(sigma/(sigma-1));
end
e_15=abs(ces_eval_15'-approx_final_15);
figure(1)
plot(k,approx_final_15)
hold on
plot(k,ces_eval_15)
title('Chebyshev polynomial approximation of order 15')
legend('Approximation','Original function')
figure(2)
plot(k,e_15)
title('Error for Chebyshev polynomial order 15')
ylabel('Error')

%Example with a polynomial of order 3:
n=3;space_x=[0:10];space_y=[0:10]
basis_3_x=chebyshev_basis_x(n,space_x);
basis_3_y=chebyshev_basis_y(n,space_y);
basis_3new=flip(basis_3_y)
for i=1:size(basis_3_y,2)
    final_vector(i)=basis_3new(1,i)*basis_3_x(1,i)
end

syms x y 
final=zeros(size(zk,1),n)
for j=1:n
for i=1:size(zk,1)
final1(i)=subs(final_vector(1,j),[x,y],[zk(i),zk(i)])
final(i,j)=final1(i)
end
end
thetas=(final'*final)\(final'*wk)

for i=1:n
approx(i)=thetas(i)*subs(final_vector(1,j),[x,y],[(2*(x-a)/(b-a))-1,(2*(y-a)/(b-a))-1])
end
cheb_approx_3=sum(approx)
approx_final_3=zeros(n,1)
h=linspace(0,10,n);k=linspace(0,10,n);
for j=1:n
approx_final_3(j,1)=double(subs(cheb_approx_3,[x,y],[h(j),k(j)]))
end
for i=1:n
ces_eval_3(i)=((1-alpha)*k(i).^((sigma-1)/sigma)+alpha*h(i).^((sigma-1)/sigma)).^(sigma/(sigma-1));
end
e_3=abs(ces_eval_3'-approx_final_3);
figure(1)
plot(k,approx_final_3)
hold on
plot(k,ces_eval_3)
title('Chebyshev polynomial approximation of order 3')
legend('Approximation','Original function')
figure(2)
plot(k,e_3)
title('Error for Chebyshev polynomial order 3')
ylabel('Error')

%%
%Isoquants plot
clear all
alpha=0.5;sigma=0.25;m=20
%Here, we are considering the maximum amount of k and h which can be used.
ces_max=((1-alpha)*10.^((sigma-1)/sigma)+alpha*10.^((sigma-1)/sigma)).^(sigma/(sigma-1));

syms k
h1=((0.5.^3)./(2-k.^3)).^(1/3)
h2=((1.^3)./(2-k.^3)).^(1/3)
h3=((2.5.^3)./(2-k.^3)).^(1/3)
h4=((5.^3)./(2-k.^3)).^(1/3)
h5=((7.5.^3)./(2-k.^3)).^(1/3)
h6=((9.^3)./(2-k.^3)).^(1/3)
h7=((9.5.^3)./(2-k.^3)).^(1/3)
k1=[1:10];
h_eval1=zeros(1,10)
for i=1:10
h_eval1(i)=subs(h1,k,k1(i))
end
h_eval2=zeros(1,10)
for i=1:10
h_eval2(i)=subs(h2,k,k1(i))
end
h_eval3=zeros(1,10)
for i=1:10
    h_eval3(i)=subs(h3,k,k1(i))
end
h_eval4=zeros(1,10)
for i=1:10
    h_eval4(i)=subs(h4,k,k1(i))
end
h_eval5=zeros(1,10)
for i=1:10
    h_eval5(i)=subs(h5,k,k1(i))
end
h_eval6=zeros(1,10)
for i=1:10
    h_eval6(i)=subs(h6,k,k1(i))
end
h_eval7=zeros(1,10)
for i=1:10
    h_eval7(i)=subs(h7,k,k1(i))
end

 plot(h_eval1,k1)
 hold on
 plot(h_eval2,k1)
 hold on
 plot(h_eval3,k1)
hold on 
 plot(h_eval4,k1)
 hold on
 plot(h_eval5,k1)
 hold on
 plot(h_eval6,k1)
hold on
plot(h_eval7,k1)
title('Isoquant graph')
legend('Isoquant percentil 5','Isoquant percentil 10','Isoquant percentil 25','Isoquant percentil 50','Isoquant percentil 75','Isoquant percentil 90','Isoquant percentil 95')

