function f=steadystate_unicountrynewB(x0,param)

kappa(2)=param(1);         %   Desutility factor
sigma(2)=param(2);      	%   Elasticity of consumption
v(2)=param(3);             %   Elasticity of work
nh(2)=param(4);
nl(2)=param(5);
phi(2)=param(6);
zeta(2)=param(7);
kd(2)=param(8)
kls(2)=param(9)
khs(2)=param(10)
theta(2)=param(11)
lambda(2)=param(12)

%Initial conditions:
cl=x0(1);
ch=x0(2);
hl=x0(3);
hh=x0(4);
w=x0(5);
r=x0(6);
%We write the equilibrium conditions:
f(1)=(1-phi(2))*cl^(-sigma(2))*lambda(2)*(w*hl*nl(2))^(-phi(2))*w*nl(2)-kappa(2)*hl^(1/v(2));
f(2)=(1-phi(2))*ch^(-sigma(2))*lambda(2)*(w*hh*nh(2))^(-phi(2))*w*nh(2)-kappa(2)*hh^(1/v(2));
f(3)=(1-theta(2))*zeta(2)*kd(2)^(-theta(2))*(nl(2)*hl+nh(2)*hh)^(theta(2))-r
f(4)=theta(2)*zeta(2)*kd(2)^(1-theta(2))*(nl(2)*hl+nh(2)*hh)^(theta(2)-1)-w
f(5)=lambda(2)*(w*hl*nl(2))^(1-phi(2))+r*kls(2)-cl;
f(6)=lambda(2)*(w*hh*nh(2))^(1-phi(2))+r*khs(2)-ch;
f=f';