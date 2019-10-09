function f=steadystate_countryb(x0,param)

kappa=param(1);         %   Desutility factor
sigma=param(2);      	%   Elasticity of consumption
v=param(3);             %   Elasticity of work
etah(2)=param(4);
etal(2)=param(5);
prob_h=param(6);
prob_l=param(7);
ksb=param(8);
theta=param(9);
z=param(10);
%Initial conditions:
hsb=x0(1);
wb=x0(2);
rb=x0(3);
%We write the equilibrium conditions:
f(1)=hsb.*-5.0+(wb.*1.0./(hsb.*wb.*(5.0./2.0)).^(1.0./5.0).*(3.99e+2./2.5e+2)+wb.*1.0./(hsb.*wb.*(7.0./2.0)).^(1.0./5.0).*1.176e-1).*1.0./(rb.*2.0+(hsb.*wb.*(5.0./2.0)).^(4.0./5.0).*(3.99e+2./5.0e+2)+(hsb.*wb.*(7.0./2.0)).^(4.0./5.0).*(2.1e+1./5.0e+2)).^(4.0./5.0);
f(2)=theta*z*2^(1-theta)*(prob_h*hsb*etah(2)+prob_l*hsb*etal(2))^(theta-1)-wb;
f(3)=(1-theta)*z*2^(-theta)*(prob_h*hsb*etah(2)+prob_l*hsb*etal(2))^(theta)-rb;
f=f';