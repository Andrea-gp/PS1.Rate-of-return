function f=steadystate_countrya(x0,param)

kappa=param(1);         %   Desutility factor
sigma=param(2);      	%   Elasticity of consumption
v=param(3);             %   Elasticity of work
etah(1)=param(4);            %   High-skill workers productivity
etal(1)=param(5);            %   Low-skill workers productivity
prob_h=param(6);        %   Distribution of high skill workers
prob_l=param(7);        %   Distribution of low skill workers
ksa=param(8);
theta=param(9);
z=param(10);
%Initial conditions:
hsa=x0(1);
wa=x0(2);
ra=x0(3);
%We write the equilibrium conditions:
f(1)=hsa.*-5.0+(wa.*1.0./((hsa.*wa)./2.0).^(1.0./5.0).*(3.61e+2./1.0e+3)+wa.*1.0./(hsa.*wa.*(1.1e+1./2.0)).^(1.0./5.0).*(2.09e+2./1.0e+3)).*1.0./(ksa.*ra+((hsa.*wa)./2.0).^(4.0./5.0).*(3.61e+2./4.0e+2)+(hsa.*wa.*(1.1e+1./2.0)).^(4.0./5.0).*(1.9e+1./4.0e+2)).^(4.0./5.0);
f(2)=theta.*z.*2.^(1-theta).*(prob_h.*hsa.*etah(1)+prob_l.*hsa.*etal(1)).^(theta-1)-wa;
f(3)=(1-theta).*z.*2.^(-theta).*(prob_h.*hsa.*etah(1)+prob_l.*hsa.*etal(1)).^(theta)-ra
%f(5)=2-kda;
%f(6)=hda-prob_h*hsa*etah(1)-prob_l*hsa*etal(1);
%f(7)=ksa-kda;
%f(7)=prob_h*(lambda(1)*(wa*hsa*etah(1))^(1-phi))+prob_l*(lambda(1)*(wa*hsa*etal(1))^(1-phi))+ra*ksa-z*((kda^(1-theta))*(hda^theta));
f=f';