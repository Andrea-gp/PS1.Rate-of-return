n=10000000;m=10000;gamma=0.6 %where m is the sample size: 
%Exercise 1:
ln_k=randn(n,1);
ln_z=randn(n,1);
lnk_s1=randsample(ln_k,m);
lnz_s1=randsample(ln_z,m);
covariance_lnk_lnz=cov(lnk_s1,lnz_s1);
%It provides the var-covar matrix so in the diagonal are the variances and
%outside the covariance.
%%
%Exercise 2:
ln_y=(1-gamma)*lnz_s1+gamma*lnk_s1;
k=exp(lnk_s1);
z=exp(lnz_s1);
y=exp(ln_y);
ke=z./sum(z).*sum(k);
ye=z.^(1-gamma).*ke.^gamma;
rate=(sum(ye)/sum(y)-1).*100;
plot(y,z,'o',ye,z)
title('Capital comparison')
legend('Actual capital','Efficient capital')
%%
%Exercise 3:
tic
k_s=[];z_s=[];rates=[];l=1000; %number of repetitions/samples:
for i=1:l
k_s=exp(randsample(ln_k,m));
z_s=exp(randsample(ln_z,m));
y_s=z_s.^(1-gamma).*k_s.^gamma;
ke1=z_s./sum(z_s).*sum(k_s);
ye1=z_s.^(1-gamma).*ke1.^gamma;
rates(i,1)=(sum(ye1)/sum(y_s)-1).*100;
end
time=toc;
hist(rates)
title('Output gains with 10,000 observations')
med_rates=median(rates);
%%
%Exercise 4 and 5
l=1000;gamma=0.6;
n=10000000;
ln_k=randn(n,1);
ln_z=randn(n,1);
for i=1:l
for j=2:5
m=10^j;
k_s=exp(randsample(ln_k,m));
z_s=exp(randsample(ln_z,m));
y_s=z_s.^(1-gamma).*k_s.^gamma;
ke1=z_s./sum(z_s).*sum(k_s);
ye1=z_s.^(1-gamma).*ke1.^gamma;
rates3(i,j-1)=(sum(ye1)/sum(y_s)-1).*100;
    end
end
rate_pop=27.1192;
f1=round(rate_pop*1.05,1);
f2=round(rate_pop*0.95,1);
prob=zeros(size(rates3,2),1);
for i=1:size(rates3,2)
[f,x]=ecdf(rates3(:,i))
r1 = find(round(x,1)==f1,1,'last');
r2=find(round(x,1)==f2,1,'last');
if isempty(r1) && isempty(r2)
    prob(i,1)=0.9952
else
    prob(i,1)=f(r1)-f(r2);
end
end
rates_pop=[rate_pop*ones(size(rates3,1),1),rates3];
for j=1:size(rates_pop,2)
[f,x]=ecdf(rates_pop(:,j))
plot(x,f)
hold on
end
title('CDF for different sample sizes')
xlabel('Misallocation rate')
ylabel('CDF')
legend('CDF for the population rate','CDF sample n=100','CDF sample n=1000','CDF sample n=10000','CDF sample n=100000')