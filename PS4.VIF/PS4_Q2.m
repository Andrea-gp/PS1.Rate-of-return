%%Exercise 2.
%PART A:
delta=0.013;beta=0.988;theta=0.679;kappa=5.24;v=2;
%Marginal of capital:
syms h k c
% MUCK=@(h,k)beta*((1-theta)*h*k^(-theta)+1-delta)-1;
% MUCL=@(h,k,c)(1/c)*theta*h^(theta-1)*k^(1-theta)+kappa*h^(1/v);
% BC=@(h,k,c)k^(1-theta)*h^theta-delta*k-c;
%[hss,kss,css]=vpasolve([MUCK,MUCL,BC],[h,k,c]);
param=[beta,delta,theta,kappa,v];
x0=[1 1 1]';
x=fsolve('steadystate',x0,optimset,param);
kss=x(1);
hss=x(2);
css=x(3);
kmin=0.25*kss;
kmax=1.25*kss;
hmin=0.25*hss;
hmax=1.25*hss;
n=100 %number of steps in which you divide the grid.
k=linspace(kmin,kmax,n);
h=linspace(hmin,hmax,n);
M=ones(n)*(-500);
for i=1:n
capital=k(i);
 for j=1:n
    inves=k(j)-(1-delta)*capital;
    prod=capital^(1-theta)*h(j)^theta;
    cons=prod-inves;
   if cons>0
   M(j,i)=log(cons)-kappa.*(h(j).^(1+1/v))./(1+1/v);
   end
end
end

v=zeros(n,1);
v1=zeros(n,1);
tolerance=2;
kapi=zeros(n,1);
tic
i=1;
while tolerance>=0.00000001 
chi=M+beta*v.*ones(1,n)
[v1,kap]=max(chi);
v1=v1';kap=kap'
tolerance=abs(v1-v)
v=v1 
kapi=kap
i=i+1
end
time_2a=toc;

%Optimal deicision rule for capital is:
k_opti=k(kapi);
h_opti=h(kapi);
i_opti=k_opti-(1-delta)*k;
y_opti=k.^(1-theta).*h_opti.^theta;
c_opti=y_opti-i_opti;

plot(k,k_opti)
hold on
plot(k,k)
title('Capital level')
legend('Optimal level (approximation)','45º line')

plot(k,v)
title('Value function')
%%
%PART B: MONOTONICITY
v=zeros(n,1);
v1=zeros(n,1);
tolerance=2;
kapi=zeros(n,1);
k_initial=-100;
tic
i=1;
while tolerance>=0.00001
chi=M+beta*v.*ones(1,n)
[v1,kap]=max(chi);
tolerance=abs(v1-v)
v=v1 
kapi=kap
if kapi>k_initial 
        k_initial = kapi;
    else
        break
    end
i=i+1  
end
time_2b=toc;
%%
%PART C: CONCAVITY
v=zeros(n,1);
v1=zeros(n,1);
tolerance=2;
kapi=zeros(n,1);
chi0=-50;
tic
i=1;
while tolerance>=0.00001
chi=M+beta*v.*ones(1,n)
[v1,kap]=max(chi);
tolerance=abs(v1-v)
v=v1 
kapi=kap
if chi>chi0 
        chi0 = chi;
    else
        break
    end
i=i+1  
end
time_2c=toc;
%%
%PART D: Taking into account local search on the decision rule
v=zeros(n,1);
v1=zeros(n,1);
tolerance=2;
kapi=zeros(n,1);
k_initial=-50;
tic
i=1;
while tolerance>=0.00001
chi=M+beta*v.*ones(1,n)
[v1,kap]=max(chi);
tolerance=abs(v1-v)
v=v1 
kapi=kap
if kapi-k_initial<0.05 
        break
    else
        k_initial=kapi
    end
i=i+1  
end
time_2d=toc;
%%
%PART E: concavity value function + monotonicity of decision rule
v=zeros(n,1);
v1=zeros(n,1);
tolerance=2;
kapi=zeros(n,1);
k_initial=-50;
chi0=-50;
tic
i=1;
while tolerance>=0.00001 
chi=M+beta*v.*ones(1,n)
[v1,kap]=max(chi);
tolerance=abs(v1-v)
v=v1 
kapi=kap
if chi>chi0 
        chi0 = chi;
    else
        break
    end
if kapi>k_initial 
        k_initial = kapi;
    else
        break
    end
i=i+1  
end
time_2e=toc;
%%
%PART F: Howard's policy iteration
v=zeros(n,1);
v1=zeros(n,1);
tolerance=2;
kapi=zeros(n,1);
k_initial=-50;
chi0=-50;
tic
i=1;
while tolerance>=0.00001 
chi=M+beta*v.*ones(1,n)
[v1,kap]=max(chi);
tolerance=abs(v1-v)
v=v1 
kapi=kap
 for l=1:n
capital=k(l)
 for j=1:n
    inves=k(kapi(j))-(1-delta)*capital;
    prod=capital^(1-theta)*h(j)^theta;
    cons=prod-inves;
   if cons>0
   M(j,l)=log(cons);
   end
end
 end
i=i+1  
end
time_2f=toc;
%%
%PART G: 5,10,20 and 50 steps i btw policy reassessment
v5=zeros(n,1);
v1=zeros(n,1);
tolerance=2;
kapi=zeros(n,1);
k_initial=-50;
chi0=-50;
tic
iter=0;
maxiter=5;
while (tolerance>=0.0001 & iter<=maxiter)
chi=M+beta*v5.*ones(1,n)
[v1,kap]=max(chi);
tolerance=abs(v1-v5)
v5=v1 
kapi=kap
for l=1:n
capital=k(l)
for j=1:n
    inves=k(kapi(j))-(1-delta)*capital;
    prod=capital^(1-theta)*h(j)^theta;
    cons=prod-inves;
   if cons>0
   M(j,l)=log(cons);
   end
end
 end
iter=iter+1  
end
time_2g5=toc;
%%
v10=zeros(n,1);
v1=zeros(n,1);
tolerance=2;
kapi=zeros(n,1);
k_initial=-50;
chi0=-50;
tic
iter=0;
maxiter=10;
while (tolerance>=0.0001 & iter<=maxiter)
chi=M+beta*v10.*ones(1,n)
[v1,kap]=max(chi);
tolerance=abs(v1-v10)
v10=v1 
kapi=kap
for l=1:n
capital=k(l)
for j=1:n
    inves=k(kapi(j))-(1-delta)*capital;
    prod=capital^(1-theta)*h(j)^theta;
    cons=prod-inves;
   if cons>0
   M(j,l)=log(cons);
   end
end
 end
iter=iter+1  
end
time_2g10=toc

%%
v20=zeros(n,1);
v1=zeros(n,1);
tolerance=2;
kapi=zeros(n,1);
k_initial=-50;
chi0=-50;
tic
iter=0;
maxiter=20;
while (tolerance>=0.0001 & iter<=maxiter)
chi=M+beta*v20.*ones(1,n)
[v1,kap]=max(chi);
tolerance=abs(v1-v20)
v20=v1 
kapi=kap
for l=1:n
capital=k(l)
for j=1:n
    inves=k(kapi(j))-(1-delta)*capital;
    prod=capital^(1-theta)*h(j)^theta;
    cons=prod-inves;
   if cons>0
   M(j,l)=log(cons);
   end
end
 end
iter=iter+1  
end
time_2g20=toc
%%
v50=zeros(n,1);
v1=zeros(n,1);
tolerance=2;
kapi=zeros(n,1);
k_initial=-50;
chi0=-50;
tic
iter=0;
maxiter=50;
while (tolerance>=0.0001 & iter<=maxiter)
chi=M+beta*v50.*ones(1,n)
[v1,kap]=max(chi);
tolerance=abs(v1-v50)
v50=v1 
kapi=kap
for l=1:n
capital=k(l)
for j=1:n
    inves=k(kapi(j))-(1-delta)*capital;
    prod=capital^(1-theta)*h(j)^theta;
    cons=prod-inves;
   if cons>0
   M(j,l)=log(cons);
   end
end
 end
iter=iter+1  
end
time_2g50=toc

plot(k,v5)
hold on
plot(k,v10)
hold on
plot(k,v20)
hold on
plot(k,v50)
title('Howard policy function iteration')
legend('After 5 iterations','After 10 iterations','After 20 iterations','After 50 iterations')