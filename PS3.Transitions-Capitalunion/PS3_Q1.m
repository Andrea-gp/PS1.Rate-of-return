% In order to compute the level of steady state capital, we need to derive
% from the function for consumption. So, we assume that beta=0.99 and the
% depreciation rate is computed as follows: at the steady state kt+1=kt=k
%so from the investment function we got it=sigma*k now dividing everything
%by yt we get it/yt=delta*(kt/yt)--> 0.25=delta*4-->delta= 1/16

%Once we have computed sigma, we state the formula that relates all terms
%where the only unknown is z:
theta=0.67; h=0.31;x0=0;beta=0.99;delta=1/16;
f=@(z)((1-theta)*(4^(-theta))*(z*h)^theta)-delta;
x=fsolve(f,x0);
%Notice that in this case, the steady state capital (kss) is 4!
z=x;
%%
%PART B:
param=[beta delta theta 2*z h];
x0=[1 1 1 1]';
x=fsolve('steadystate',x0,optimset,param);
kss=x(1);
css=x(3);

%%
%PART C:
param=[beta delta theta 2*z h];
x1=[kss kss css css]';
    J=jacobian('steadystate',x1,param); % transition matrix
    MA=[J(1,1) J(1,3);
        J(2,1) J(2,3)];
    MB=[J(1,2) J(1,4);
        J(2,2) J(2,4)];
    MG=-inv(MA)*MB;
    [Mvec,Meig]=eig(MG);
    Meig2=diag(Meig);
    [Meig3,Mord]=sort(abs(Meig2));
    Meig4=Meig2(Mord);
    Mlambda=diag(Meig4);
    MM=Mvec(:,Mord);
    Mm=inv(MM);
    
    stab=-Mm(2,1)/Mm(2,2); % stabilizing constant
   
%Solution using the nonlinear model imposing stability at all points in time
k00=4*0.8; %where 0.8 is a shock 
nobs=100;
k0=k00;
c0=css+stab*(k0-kss);
vk4=[];vk4=[vk4;k0];
vc4=[];vc4=[vc4;c0];
    for i=1:nobs
        k1=kss+MG(1,1)*(k0-kss)+MG(1,2)*(c0-css);
        c1=css+stab*(k0-kss);
        c0=c1;
        k0=k1;
        vk4=[vk4;k0];
        vc4=[vc4;c0];
    end
y0=vk4.^(1-theta)*(2*z*h)^theta;
i0=y0-vc4;

%Graphs for the transitions:
subplot(2,2,1)
plot(vk4);
title('Solution using the linear approx');
ylabel('Capital stock');
xlabel('Periods');
subplot(2,2,2)
plot(vc4);
title('Solution using the linear approx');
ylabel('Consumption');
xlabel('Periods');  
subplot(2,2,3)
plot(y0);
title('Solution using the linear approx');
ylabel('Output');
xlabel('Periods');  
subplot(2,2,4)
plot(i0);
title('Solution using the linear approx');
ylabel('Investment=savings');
xlabel('Periods');  

%%
%PART D:
vk41_1=vk4(10,1)*0.8 %where 0.8 is the shock
param=[beta delta theta z h];
x0=[1 1 1 1]';
x=fsolve('steadystate',x0,optimset,param);
kss=x(1);
css=x(3);
param=[beta delta theta z h];
x1=[kss kss css css]';
J=jacobian('steadystate',x1,param); % transition matrix
    MA=[J(1,1) J(1,3);
        J(2,1) J(2,3)];
    MB=[J(1,2) J(1,4);
        J(2,2) J(2,4)];
    MG=-inv(MA)*MB;
    [Mvec,Meig]=eig(MG);
    Meig2=diag(Meig);
    [Meig3,Mord]=sort(abs(Meig2));
    Meig4=Meig2(Mord);
    Mlambda=diag(Meig4);
    MM=Mvec(:,Mord);
    Mm=inv(MM);
    stab1=-Mm(2,1)/Mm(2,2); % stabilizing constant

k0=vk41_1;
c0=css+stab1*(k0-kss);
vk41=[];vk41=[vk41;k0];
vc41=[];vc41=[vc41;c0];
    for i=1:nobs
        k1=kss+MG(1,1)*(k0-kss)+MG(1,2)*(c0-css);
        c1=css+stab1*(k0-kss);
        c0=c1;
        k0=k1;
        vk41=[vk41;k0];
        vc41=[vc41;c0];
    end
c_partd=zeros(111,1);c_partd(1:10,1)=vc4(1:10);c_partd(11:end,1)=vc41
k_partd=zeros(111,1);k_partd(1:10,1)=vk4(1:10);k_partd(11:end,1)=vk41
y01=k_partd.^(1-theta)*(z*h)^theta;
i01=y01-c_partd;

%Grpahs for the transitions:
subplot(2,2,1)
plot(k_partd);
title('Solution using the linear approx');
ylabel('Capital stock');
xlabel('Periods');
subplot(2,2,2)
plot(c_partd);
title('Solution using the linear approx');
ylabel('Consumption');
xlabel('Periods');  
subplot(2,2,3)
plot(y01);
title('Solution using the linear approx');
ylabel('Output');
xlabel('Periods');  
subplot(2,2,4)
plot(i01);
title('Solution using the linear approx');
ylabel('Investment=savings');
xlabel('Periods');   

%%
%PART E:
tauc=0.057;beta=0.99;delta=1/16;theta=0.67;z=1.0768;h=0.31;
param=[beta delta theta 2*z h tauc];
x0=[1 1 1 1]';
x=fsolve('steadystate_taxonc',x0,optimset,param);
kss=x(1);
css=x(3);
param=[beta delta theta 2*z h];
x1=[kss kss css css]';
    J=jacobian('steadystate',x1,param); % transition matrix
    MA=[J(1,1) J(1,3);
        J(2,1) J(2,3)];
    MB=[J(1,2) J(1,4);
        J(2,2) J(2,4)];
    MG=-inv(MA)*MB;
    [Mvec,Meig]=eig(MG);
    Meig2=diag(Meig);
    [Meig3,Mord]=sort(abs(Meig2));
    Meig4=Meig2(Mord);
    Mlambda=diag(Meig4);
    MM=Mvec(:,Mord);
    Mm=inv(MM);
    
    stab=-Mm(2,1)/Mm(2,2); % stabilizing constant
   
%Solution using the nonlinear model imposing stability at all points in time
k00=4*0.8; %where 0.8 is a shock 
nobs=100;
k0=k00;
c0=css+stab*(k0-kss);
vk4=[];vk4=[vk4;k0];
vc4=[];vc4=[vc4;c0];
    for i=1:nobs
        k1=kss+MG(1,1)*(k0-kss)+MG(1,2)*(c0-css);
        c1=css+stab*(k0-kss);
        c0=c1;
        k0=k1;
        vk4=[vk4;k0];
        vc4=[vc4;c0];
    end
y0=vk4.^(1-theta)*(2*z*h)^theta;
i0=y0-vc4;

%Graphs for exercise 1.e:
subplot(2,2,1)
plot(vk4);
title('Solution using the linear approx with taxonc');
ylabel('Capital stock');
xlabel('Periods');
subplot(2,2,2)
plot(vc4);
title('Solution using the linear approx with taxonc');
ylabel('Consumption');
xlabel('Periods');  
subplot(2,2,3)
plot(y0);
title('Solution using the linear approx');
ylabel('Output');
xlabel('Periods');  
subplot(2,2,4)
plot(i0);
title('Solution using the linear approx');
ylabel('Investment=savings');
xlabel('Periods');  

taucy=0,4287;beta=0.99;delta=1/16;theta=0.67;z=1.0768;h=0.31;
param=[beta delta theta 2*z h taucy];
x0=[1 1 1 1]';
x=fsolve('steadystate_taxony',x0,optimset,param);
kss=x(1);
css=x(3);
param=[beta delta theta 2*z h];
x1=[kss kss css css]';
    J=jacobian('steadystate',x1,param); % transition matrix
    MA=[J(1,1) J(1,3);
        J(2,1) J(2,3)];
    MB=[J(1,2) J(1,4);
        J(2,2) J(2,4)];
    MG=-inv(MA)*MB;
    [Mvec,Meig]=eig(MG);
    Meig2=diag(Meig);
    [Meig3,Mord]=sort(abs(Meig2));
    Meig4=Meig2(Mord);
    Mlambda=diag(Meig4);
    MM=Mvec(:,Mord);
    Mm=inv(MM);
    
    stab=-Mm(2,1)/Mm(2,2); % stabilizing constant
   
%Solution using the nonlinear model imposing stability at all points in time
k00=4*0.8; %where 0.8 is a shock 
nobs=100;
k0=k00;
c0=css+stab*(k0-kss);
vk4=[];vk4=[vk4;k0];
vc4=[];vc4=[vc4;c0];
    for i=1:nobs
        k1=kss+MG(1,1)*(k0-kss)+MG(1,2)*(c0-css);
        c1=css+stab*(k0-kss);
        c0=c1;
        k0=k1;
        vk4=[vk4;k0];
        vc4=[vc4;c0];
    end
y0=vk4.^(1-theta)*(2*z*h)^theta;
i0=y0-vc4;

%Graphs for exercise 1.e:
subplot(2,2,1)
plot(vk4);
title('Solution using the linear approx with taxony');
ylabel('Capital stock');
xlabel('Periods');
subplot(2,2,2)
plot(vc4);
title('Solution using the linear approx with taxony');
ylabel('Consumption');
xlabel('Periods');  
subplot(2,2,3)
plot(y0);
title('Solution using the linear approx');
ylabel('Output');
xlabel('Periods');  
subplot(2,2,4)
plot(i0);
title('Solution using the linear approx');
ylabel('Investment=savings');
xlabel('Periods');  
