function f=steadystate_bicountry(x0,param)

kappa=param(1);         %   Desutility factor
sigma=param(2);      	%   Elasticity of consumption
v=param(3);             %   Elasticity of work
etah(1)=param(4);           %   High-skill workers country a productivity
etal(1)=param(5);           %   Low-skill workers country a productivity
etah(2)=param(6);           %   High-skill workers country b productivity
etal(2)=param(7);           %   Low-skill workers country b productivity
prob_l=param(8);        %   Distribution of low skill workers
prob_h=param(9);        %   Distribution of high skill workers
theta=param(10);
z=param(11);
%Initial conditions:
hsa=x0(1);
hsb=x0(2);
ksa=x0(3);
ksb=x0(4);
wa=x0(5);
wb=x0(6);
ra=x0(7);
rb=x0(8);
%We write the equilibrium conditions:
f(1)=(ra-rb).*1.0./(ksa.*ra-rb.*(ksa-2.0)+((hsa.*wa)./2.0).^(4.0./5.0).*(3.61e+2./4.0e+2)+(hsa.*wa.*(1.1e+1./2.0)).^(4.0./5.0).*(1.9e+1./4.0e+2)).^(4.0./5.0);
f(2)=-(ra-rb).*1.0./(ksb.*rb-ra.*(ksb-2.0)+(hsb.*wb.*(5.0./2.0)).^(4.0./5.0).*(3.99e+2./5.0e+2)+(hsb.*wb.*(7.0./2.0)).^(4.0./5.0).*(2.1e+1./5.0e+2)).^(4.0./5.0);
f(3)=hsa.*-5.0+(wa.*1.0./((hsa.*wa)./2.0).^(1.0./5.0).*(3.61e+2./1.0e+3)+wa.*1.0./(hsa.*wa.*(1.1e+1./2.0)).^(1.0./5.0).*(2.09e+2./1.0e+3)).*1.0./(ksa.*ra-rb.*(ksa-2.0)+((hsa.*wa)./2.0).^(4.0./5.0).*(3.61e+2./4.0e+2)+(hsa.*wa.*(1.1e+1./2.0)).^(4.0./5.0).*(1.9e+1./4.0e+2)).^(4.0./5.0);
f(4)=hsb.*-5.0+(wb.*1.0./(hsb.*wb.*(5.0./2.0)).^(1.0./5.0).*(3.99e+2./2.5e+2)+wb.*1.0./(hsb.*wb.*(7.0./2.0)).^(1.0./5.0).*1.176e-1).*1.0./(ksb.*rb-ra.*(ksb-2.0)+(hsb.*wb.*(5.0./2.0)).^(4.0./5.0).*(3.99e+2./5.0e+2)+(hsb.*wb.*(7.0./2.0)).^(4.0./5.0).*(2.1e+1./5.0e+2)).^(4.0./5.0);
f(5)=ra-((1-theta).*z.*(ksa+ksb-2).^(-theta).*(prob_h.*hsa.*etah(1)+prob_l.*hsa.*etal(1)).^(theta));
f(6)=wa-((theta.*z.*((ksa+ksb-2).^(1-theta)).*(prob_h.*hsa.*etah(1)+prob_l.*hsa.*etal(1)).^(theta-1)));
f(7)=rb-((1-theta).*z.*(ksa+ksb-2).^(-theta).*(prob_h.*hsb.*etah(2)+prob_l.*hsb.*etal(2)).^(theta)); 
f(8)=wb-(theta.*z.*((ksa+ksb-2).^(1-theta)).*(prob_h.*hsb.*etah(2)+prob_l.*hsb.*etal(2)).^(theta));
f=f';