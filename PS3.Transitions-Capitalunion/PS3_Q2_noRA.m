%Exercise 2: NO REPRESENTATIVE AGENT ECONOMY
%Part 2.1 Economy A
sigma=[0.8 0.8];kappa=[5 5];v=[1 1];phi=[0.2 0.2];lambda=[0.95 0.84];nl=[5.5 3.5];zeta=[1 1]; nh=[0.5 2.5]; theta=[0.6 0.6]; kd=[2 2]; kls=[1 1]; khs=[1 1];
syms cl ch hl hh w r
param=[kappa(1) sigma(1) v(1) nh(1) nl(1) phi(1) zeta(1) kd(1) kls(1) khs(1) theta(1) lambda(1)];
    x0=[1 1 1 1 1 1]';
    x=fsolve('steadystate_unicountrynew',x0,optimset,param);
    cl_eq=x(1);
    ch_eq=x(2);
    hl_eq=x(3);
    hh_eq=x(4);
    w_eq=x(5);
    r_eq=x(6);
    
  %Economy B
 syms cl ch hl hh w r
param=[kappa(2) sigma(2) v(2) nh(2) nl(2) phi(2) zeta(2) kd(2) kls(2) khs(2) theta(2) lambda(2)];
    x0=[1 1 1 1 1 1]';
    x=fsolve('steadystate_unicountrynewB',x0,optimset,param);
    cl_eq2=x(1);
    ch_eq2=x(2);
    hl_eq2=x(3);
    hh_eq2=x(4);
    w_eq2=x(5);
    r_eq2=x(6); 
%%    
%Part 2.2: Capital union
kappa=5;
v=1;
sigma=0.8;
nl=5.5;
nh=0.5;
zeta=1;
theta=0.6;
KAd=2;
KLAmax=1;
KHAmax=1;
lambda=0.95;
phi=0.2; 
kappaB=5;
vB=1;
sigmaB=0.8;
nlb=3.5;
nhb=2.5;
thetaB=0.6;
KBd=2;
KLBmax=1;
KHBmax=1;
lambB=0.84;
phiB=0.2; 
    syms cla cha hla hha wa ra kla kha clb chb hlb hhb wb rb klb khb
   param=[kappa v sigma nl nh zeta theta KAd KLAmax KHAmax lambda phi kappaB vB sigmaB nlb nhb thetaB KBd KLBmax KHBmax lambB phiB]
   x0=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';
    x=fsolve('steadystate_union',x0,optimset,param);
    claa=x(1);
    chaa=x(2);
    hlaa=x(3);
    hhaa=x(4);
    waa=x(5);
    raa=x(6);
    klaa=x(7);
    khaa=x(8);
    clbb=x(9);
    chbb=x(10);
    hlbb=x(11);
    hhbb=x(12);
    wbb=x(13);
    rbb=x(14);
    klaa=x(15);
    khaa=x(16);