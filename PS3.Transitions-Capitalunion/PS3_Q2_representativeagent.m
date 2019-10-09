%QUESTION 2:
%Part 2.1:
%To solve this question, we have provide a random distribution for low and
%high skill workers:
%Country A
sigma=0.8;kappa=5;v=1;phi=0.2;lambda=[0.95,0.84];etal=[5.5,3.5];etah=[0.5,2.5];z=1;theta=0.6;prob_h=0.95;prob_l=1-prob_h;
ksa=2;
syms kda hda hsa wa ra
%For households:
ca=prob_h*(lambda(1)*(wa*hsa*etah(1))^(1-phi))+prob_l*(lambda(1)*(wa*hsa*etal(1))^(1-phi))+ra*ksa;
consa=(((ca)^(1-sigma))/(1-sigma))-kappa*(hsa^(1+(1/v)))/(1+(1/v));
MUCa=matlabFunction(diff(consa,hsa));
%MUC1a=matlabFunction(diff(consa,ksa));
%For firms:
profa=(z*(kda^(1-theta))*(prob_h*hsa*etah(1)-prob_l*hsa*etal(1))^theta)-wa*(prob_h*hsa*etah(1)-prob_l*hsa*etal(1))-ra*kda;
MPKa=matlabFunction(diff(profa,kda));
MPHa=matlabFunction(diff(profa,hda));
param=[kappa sigma v etah(1) etal(1) prob_h prob_l ksa theta z];
x0=[1 1 1 ]';
x=fsolve('steadystate_countrya',x0,optimset,param)
hs_eqa=x(1);
w_eqa=x(2);
r_eqa=x(3);
%%
%Country B
%sigma=0.8;kappa=5;v=1;phi=0.2;lambda=0.95;nl=5.5;nh=0.5;z=1;theta=0.6;prob_h=0.95;prob_l=1-prob_h;
syms kdb hdb hsb wb rb
ksb=2;
%For households:
cb=prob_h*(lambda(2)*(wb*hsb*etah(2))^(1-phi))+prob_l*(lambda(2)*(wb*hsb*etal(2))^(1-phi))+rb*ksb;
consb=(((cb)^(1-sigma))/(1-sigma))-kappa*(hsb^(1+(1/v)))/(1+1/v);
MUCb=matlabFunction(diff(consb,hsb));
%MUC1b=matlabFunction(diff(consb,ksb));
%For firms:
profb=(z*(kdb^(1-theta))*(prob_h*hsb*etah(2)-prob_l*hsb*etal(2))^theta)-wb*(prob_h*hsb*etah(2)-prob_l*hsb*etal(2))-rb*kdb;
MPKb=matlabFunction(diff(profb,kdb));
MPHb=matlabFunction(diff(profb,hdb));
param=[kappa sigma v etah(2) etal(2) prob_h prob_l ksb theta z];
x0=[1 1 1]';
x=fsolve('steadystate_countryb',x0,optimset,param)
hs_eqb=x(1);
w_eqb=x(2);
r_eqb=x(3);
%%
%Part 2.2:
sigma=0.8;kappa=5;v=1;phi=0.2;lambda=[0.95,0.84];etal=[5.5,3.5];etah=[0.5,2.5];z=1;theta=0.6;prob_h=0.95;prob_l=1-prob_h;
kda=2;kdb=2;khat=2;
syms ksa ksb hsa hsb wa wb ra rb
%For households:
c22a=prob_h*(lambda(1)*(wa*hsa*etah(1))^(1-phi))+prob_l*(lambda(1)*(wa*hsa*etal(1))^(1-phi))+ra*ksa+rb*(khat-ksa);
cons22a=(((c22a)^(1-sigma))/(1-sigma))-kappa*(hsa^(1+(1/v)))/(1+1/v);
MUCa=matlabFunction(diff(cons22a,hsa));
MUC1a=matlabFunction(diff(cons22a,ksa));
c22b=prob_h*(lambda(2)*(wb*hsb*etah(2))^(1-phi))+prob_l*(lambda(2)*(wb*hsb*etal(2))^(1-phi))+rb*ksb+ra*(khat-ksb);
cons22b=(((c22b)^(1-sigma))/(1-sigma))-kappa*(hsb^(1+(1/v)))/(1+1/v);
MUCb=matlabFunction(diff(cons22b,hsb));
MUC1b=matlabFunction(diff(cons22b,ksb));
param=[kappa sigma v etah(1) etal(1) etah(2) etal(2) prob_l prob_h theta z];
x0=[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]';
x=fsolve('steadystate_bicountry',x0,optimset,param);
hsa_eq=x(1);
hsb_eq=x(2);
ksa_eq=x(3);
ksb_eq=x(4);
wa_eq=x(5);
wb_eq=x(6);
ra_eq=x(7);
rb_eq=x(8);

%% 
%Part 2.3:
syms t ha hb ka kb
ra=MPKa;rb=MPKb;wa=MPHa;wb=MPHb;g=0.3846;
c222a=prob_h*(lam_a*((1-t)*wa*ha*nha)^(1-phi))+prob_l*(lam_a*((1-t)*wa*ha*nla)^(1-phi))+ra*ka+rb*(khat-ka);
cons222a=(((c222a)^(1-sigma))/(1-sigma))-kappa*(hsa^(1+(1/v)))/(1+1/v);
c222b=prob_h*(lam_b*((1-t)*wb*hb*nhb)^(1-phi))+prob_l*(lam_b*((1-t)*wb*hb*nlb)^(1-phi))+rb*kb+ra*(khat-kb);
cons222b=(((c222b)^(1-sigma))/(1-sigma))-kappa*(hsb^(1+(1/v)))/(1+1/v);
SWF=cons222a+cons222b+g;
MSWF_ka=matlabFunction(diff(SWF,ka));
MSWF_kb=matlabFunction(diff(SWF,kb));
MSWF_ha=matlabFunction(diff(SWF,ha));
MSWF_hb=matlabFunction(diff(SWF,hb));
param=[kappa sigma v nha nla nhb nlb prob_h prob_l phi theta lam_a lam_b];
x0=[1 1 1 1 1]';
x=fsolve('steadystate_social',x0,optimset,param);
hs_eq=x(1);
hd_eq=x(2);
ks_eq=x(3);
kd_eq=x(4);
t_eq=x(5);