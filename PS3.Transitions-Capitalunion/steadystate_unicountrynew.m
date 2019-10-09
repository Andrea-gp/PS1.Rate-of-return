function f=steadystate_unicountrynew(x0,param)

kappa(1)=param(1);         %   Desutility factor
sigma(1)=param(2);      	%   Elasticity of consumption
v(1)=param(3);             %   Elasticity of work
nh(1)=param(4);
nl(1)=param(5);
phi(1)=param(6);
zeta(1)=param(7);
kd(1)=param(8)
kls(1)=param(9)
khs(1)=param(10)
theta(1)=param(11)
lambda(1)=param(12)

%Initial conditions:
cl=x0(1);
ch=x0(2);
hl=x0(3);
hh=x0(4);
w=x0(5);
r=x0(6);
%We write the equilibrium conditions:
f(1)=(1-phi(1))*cl^(-sigma(1))*lambda(1)*(w*hl*nl(1))^(-phi(1))*w*nl(1)-kappa(1)*hl^(1/v(1));
f(2)=(1-phi(1))*ch^(-sigma(1))*lambda(1)*(w*hh*nh(1))^(-phi(1))*w*nh(1)-kappa(1)*hh^(1/v(1));
f(3)=(1-theta(1))*zeta(1)*kd(1)^(-theta(1))*(nl(1)*hl+nh(1)*hh)^(theta(1))-r
f(4)=theta(1)*zeta(1)*kd(1)^(1-theta(1))*(nl(1)*hl+nh(1)*hh)^(theta(1)-1)-w
f(5)=lambda(1)*(w*hl*nl(1))^(1-phi(1))+r*kls(1)-cl;
f(6)=lambda(1)*(w*hh*nh(1))^(1-phi(1))+r*khs(1)-ch;
f=f';