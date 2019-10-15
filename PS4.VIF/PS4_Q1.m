%PS4. VALUE FUNCTION ITERATION. Exercise 1
%Part a
%Firstly, we are determining the grid where we are going to work on. To do
%so, we need to find the capital steady state level
delta=0.013;h=1;beta=0.988;theta=0.679;
kss=(delta/(beta*h*(1-theta)))^(-1/theta);
kmin=0.25*kss;
kmax=1.25*kss;
n=100 %number of steps in which you divide the grid.
k=linspace(kmin,kmax,n);
M=ones(n)*(-500);
%Steps 2,3 and 4.
%We fixed capital today with the first for, and then with the second one,
%we move capital for tomorrow. 
 for i=1:n
capital=k(i)
 for j=1:n
    inves=k(j)-(1-delta)*capital;
    prod=capital^(1-theta);
    cons=prod-inves;
   if cons>0
   M(j,i)=log(cons);
   end
end
 end
%Step 5:
v=zeros(n,1);
v1=zeros(n,1);
tolerance=2;
kapi=zeros(n,1);
tic
i=1;
while tolerance>=0.00001 
chi=M+beta*v.*ones(1,n)
[v1,kap]=max(chi);
v1=v1';kap=kap'
tolerance=abs(v1-v)
v=v1 
kapi=kap
i=i+1  
end
time_1a=toc;
%Optimal deicision rule for capital is:
k_opti=k(kapi);
i_opti=k_opti-(1-delta)*k;
y_opti=k.^(1-theta);
c_opti=y_opti-i_opti;

plot(k,k_opti)
hold on
plot(k,k)
hold on
plot(k,v)
legend('Optimal capital level (approx)','45º line','Value function')
%%
%PART B: MONOTONICITY
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
if kapi>k_initial 
        k_initial = kapi;
    else
        break
    end
i=i+1  
end
time_1b=toc;
plot(k,v)
%%
% PART C: CONCAVITY
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
time_1c=toc;
plot(k,v)

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
time_1d=toc

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
time_1e=toc

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
for j=1:5
    inves=k(kapi(j))-(1-delta)*capital;
    prod=capital^(1-theta);
    cons=prod-inves;
   if cons>0
   M(j,l)=log(cons);
   end
end
 end
i=i+1  
end
time_1f=toc;
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
    prod=capital^(1-theta);
    cons=prod-inves;
   if cons>0
   M(j,l)=log(cons);
   end
end
 end
iter=iter+1  
end
time_1g5=toc;
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
    prod=capital^(1-theta);
    cons=prod-inves;
   if cons>0
   M(j,l)=log(cons);
   end
end
 end
iter=iter+1  
end
time_1g10=toc

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
    prod=capital^(1-theta);
    cons=prod-inves;
   if cons>0
   M(j,l)=log(cons);
   end
end
 end
iter=iter+1  
end
time_1g20=toc
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
    prod=capital^(1-theta);
    cons=prod-inves;
   if cons>0
   M(j,l)=log(cons);
   end
end
 end
iter=iter+1  
end
time_1g50=toc

plot(k,v5)
hold on
plot(k,v10)
hold on
plot(k,v20)
hold on
plot(k,v50)
title('Howard policy function iteration')
legend('After 5 iterations','After 10 iterations','After 20 iterations','After 50 iterations')