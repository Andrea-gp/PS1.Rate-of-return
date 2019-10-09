function f=steadystate(x0,param)

beta=param(1);         %   Disocunt factor
delta=param(2);      	%   Depreciation rate
theta=param(3);             %   Labor elasticity

z=param(4); %Shock
h=param(5)  %Exogenos labor
%Initial conditions:
k1=x0(1);
k=x0(2);
c1=x0(3);
c=x0(4);
%We write the equilibrium conditions:
%1. Resource constraint
f(1)=(k^(1-theta))*((z*h)^theta)-c-k1+(1-delta)*k
%2. Derivative wrt k+1
f(2)=(-1/c)+beta*(1/c1)*((1-theta)*(k1^(-theta))*((z*h)^theta)+1-delta);
%3. In the steady state k1=k and c1=c, so k1-k=0 and c1-c=0
f(3)=k1-k;
f(4)=c1-c;
f=f';