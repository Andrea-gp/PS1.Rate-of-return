function val=value(k,param,thetas)
kmin=param(1);
kmax=param(2);
n=param(3);
k=2*(k-kmin)/(kmax-kmin)-1;
val=chebyshev_approx(k,n)*thetas;
%This function is the last step in the chebyshev algorithm, when we are
%evaluating our function (the value function) when our k is within [kmin,kmax]
%and multiplying by the basis.
end