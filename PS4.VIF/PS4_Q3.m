delta=0.013;h=1;beta=0.988;theta=0.679;
kss=(delta/(beta*h*(1-theta)))^(-1/theta);
kmin=0.25*kss;
kmax=1.25*kss;
n=100 %number of steps in which you divide the grid.
%%
tic
for j=1:n
zk(j,1)=cos((2*j-1)/(2*n)*pi); %they are the Chebyshev interpolation nodes in [-1,1]
end
a=kmin;b=kmax;int=[a,b];
xk=flip((zk+1)*((int(2)-int(1))/2)+a); %we adjust the nodes to the interval [kmin,kmax]
v=zeros(n,1);V1=zeros(n,1);kapi=zeros(n,1);
basis_k=chebyshev_approx(zk,n);%we are computing the Chebyshev polynomial basis
th_initial=basis_k\v; %we pin down our thetas by relating our basis with the value function
tolerance=1; %initial level of tolerance
maxiter=100; %maximum amount of iterations 
i=1;
while (tolerance>0.00001 & i<=maxiter)
    k_initial=xk(1);
        for j=1:n
        param=[delta beta theta n kmin kmax xk(j)];
        kapi(j)=fminunc(@v1,k_initial,[],param,th_initial); 
        k_initial=kapi(j);
        V1(j)=-v1(kapi(j),param,th_initial)
    end
    th_final=basis_k\V1;
    tolerance=abs(V1-v);
    v=V1;
    th_initial=th_final;
    i=i+1
end
time_3=toc;    
plot(xk,v)
title('Value function by Chebyshev approx')