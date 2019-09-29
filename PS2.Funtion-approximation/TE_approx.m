function [expTaylor]=TE_approx(f,a,n)
inter=[];
syms x
for i=1:n
    deriv=diff(f,x,i)
    inter=[inter subs(deriv,a)*(x-a)^i/factorial(i)]
end
    expTaylor=sum(inter)
end