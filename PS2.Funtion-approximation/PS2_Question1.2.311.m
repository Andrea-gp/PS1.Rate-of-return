%Question 1.1 (Taylor expansion)
syms x
f=x^0.321
a=1;n=5;%inter=[];
[expTaylor]=TE_approx(f,a,n)
A=[1,2,5,20]
for j=1:size(A,2)
   TE(j)=TE_approx(f,a,A(j)) 
   for i=0:0.1:4
       TA(j)=subs(TE(j),i)
   end
end

fplot([TE(1) TE(2) TE(3) TE(4) TE])
xlim([0 4]);ylim([-2 7])
title('Taylor Approximation')
legend('TE order 1','TE order 2','TE order 5','TE order 20')

%Here is the command for the Taylor approximation already implemented in
%MATLAB
% t1=taylor(f,x,'ExpansionPoint',1,'Order',1)
% t2=taylor(f,x,'ExpansionPoint',1,'Order',2)
% t5=taylor(f,x,'ExpansionPoint',1,'Order',5)
% t20=taylor(f,x,'ExpansionPoint',1,'Order',20)
%Plot for the Taylor approximation:
%fplot([t1 t2 t5 t20 f])
%xlim([0 4]);ylim([-2 7])

%%
%Question 1.2: The implemented function for TE is not working, here we just
%write it down 
%syms y
%g=(y+abs(y))/2
% T1=taylor(g,y,'ExpansionPoint',2,'Order',1)
% T2=taylor(g,y,'ExpansionPoint',2,'Order',2)
% T5=taylor(g,y,'ExpansionPoint',2,'Order',5)
% T20=taylor(g,y,'ExpansionPoint',2,'Order',20)

%Our suggestion to make it work:
%g=y+-y/2-->g1=y-y/2=0 and g2=y+y/2=y so we are just considering g=y, hence
%the function would be the same for whatever order, otherwise there is no derivative
%Here we apply directly the Taylor command:
syms y
g=(y+abs(y))/2
g1 = y
T1=taylor(g1,y,'ExpansionPoint',2,'Order',1)
T2=taylor(g1,y,'ExpansionPoint',2,'Order',2)
T5=taylor(g1,y,'ExpansionPoint',2,'Order',5)
T20=taylor(g1,y,'ExpansionPoint',2,'Order',20)
fplot([g g1])
xlim([-2 6])
ylim([-4 6])
title('Taylor approximation for the ramp function')
legend("Original function","TE approx")
%%
%Question 1.3.1:
s=3;n=5;m=10; %order of the polynomial
x = linspace(-1,1,20); %We are generating 20 points btw -1 and 1
%For function 1
f1=exp(1./x);
theta1_3 = polyfit(x,f1,s);
theta1_5 = polyfit(x,f1,n);
theta1_10 = polyfit(x,f1,m);
y1_3 = polyval(theta1_3,x);
y1_5 = polyval(theta1_5,x);
y1_10 = polyval(theta1_10,x);
e1_3=abs(f1-y1_3);
e1_5=abs(f1-y1_5);
e1_10=abs(f1-y1_10);
figure(1)
plot(x,f1,'LineWidth',2);
hold on
plot(x,y1_3,'LineWidth',2)
hold on
plot(x,y1_5,'LineWidth',2)
hold on
plot(x,y1_10,'LineWidth',2)
title('Evenly space interpolation for Function 1')
legend('Original function','Order 3','Order 5','Order 10')

figure(2)
plot(x,e1_3,'g','LineWidth',2);
hold on
plot(x,e1_5,'r','LineWidth',2)
hold on
plot(x,e1_10,'b','LineWidth',2)
title('Error terms Funtion 1')
legend("Error terms polynomial order 3","Error terms polynomial order 5","Error terms polynomial order 10")

%For Runge Function
f2=1./(1+25*x.^2);
theta2_3 = polyfit(x,f2,s);
theta2_5 = polyfit(x,f2,n);
theta2_10 = polyfit(x,f2,m);
y2_3 = polyval(theta2_3,x);
y2_5 = polyval(theta2_5,x);
y2_10 = polyval(theta2_10,x);
e2_3=abs(f2-y2_3);
e2_5=abs(f2-y2_5);
e2_10=abs(f2-y2_10);
figure(3)
plot(x,f2);
hold on
plot(x,y2_3)
hold on
plot(x,y2_5)
hold on
plot(x,y2_10)
title('Evenly space interpolation Runge Function')
legend('Original function','Order 3','Order 5','Order 10')

figure(4)
plot(x,e2_3,'g','LineWidth',2);
hold on
plot(x,e2_5,'r','LineWidth',2)
hold on
plot(x,e2_10,'b','LineWidth',2)
title('Error terms for Runge Function')
legend("Error terms polynomial order 3","Error terms polynomial order 5","Error terms polynomial order 10")

%For Ramp Function
f3=(x+abs(x))/2;
theta3_3 = polyfit(x,f3,s);
theta3_5 = polyfit(x,f3,n);
theta3_10 = polyfit(x,f3,m);
y3_3 = polyval(theta3_3,x);
y3_5 = polyval(theta3_5,x);
y3_10 = polyval(theta3_10,x);
e3_3=abs(f3-y3_3);
e3_5=abs(f3-y3_5);
e3_10=abs(f3-y3_10);
figure(5)
plot(x,f3);
hold on
plot(x,y3_3)
hold on
plot(x,y3_5)
hold on
plot(x,y3_10)
title('Evenly space interpolation Ramp Function')
legend('Original function','Order 3','Order 5','Order 10')

figure(6)
plot(x,e3_3,'g','LineWidth',2);
hold on
plot(x,e3_5,'r','LineWidth',2)
hold on
plot(x,e3_10,'b','LineWidth',2)
title('Error terms for Ramp Function')
legend("Error terms polynomial order 3","Error terms polynomial order 5","Error terms polynomial order 10")
