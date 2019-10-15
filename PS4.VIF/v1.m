function residual=v1(kapi,param,thetas)
delta=param(1);
beta=param(2);
theta=param(3);
n=param(4);
kmin=param(5);
kmax=param(6);
k=param(7);
kapi=sqrt(kapi.^2);
v=value(kapi,[kmin kmax n],thetas);
c=k.^(1-theta)+(1-delta)*k-kapi;
if c<0
    util=-500
else
    util=log(c);
end
residual=-(util+beta*v);
%This function is computing the Bellman equation, where the only unknown is
%the capital for next period (called kapi in this case). 
end