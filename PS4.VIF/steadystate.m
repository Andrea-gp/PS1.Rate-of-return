function f=steadystate(x0,param)

beta=param(1);         %   Disocunt factor
delta=param(2);      	%   Depreciation rate
theta=param(3);             %   Labor elasticity
kappa=param(4); %Desiutility of working
v=param(5)  
%Initial conditions:
k=x0(1);
h=x0(2);
c=x0(3);
%We write the equilibrium conditions:
%1. Resource constraint
f(1)=(k^(1-theta))*(h^theta)-delta*k-c;
%2. Derivative wrt k+1
f(2)=beta*((1-theta)*h*(k^(-theta))+1-delta)-1;
%3. Derivative wrt labor
f(3)=(1/c)*theta*(h^(theta-1))*k^(1-theta)-kappa*h^(1/v);
f=f';