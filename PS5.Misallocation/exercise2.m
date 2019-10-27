%%
clear
n=10000000;mu=[0,0];sigma=[1,0;0,1];
ln_k=randn(n,1);
ln_z=randn(n,1);
gamma=0.8;
ln_y=(1-gamma)*ln_z+gamma*ln_k;
k=exp(ln_k);
z=exp(ln_z);
y=exp(ln_y);
%Efficient capital according to the definition given in the slides
ke=z./sum(z).*sum(k);
%we can see that capital is concentrated in the firms with low productivity
%that is not efficient
%Exercise 5:
ye=z.^(1-gamma).*ke.^gamma;
% Ye=sum(ye);
% Y=sum(y);
rate_pop=(sum(ye)/sum(y)-1).*100;
fprintf('The country would gain %4.2f when the capital is allocated efficiently\n',rate_pop)
%%
%Exercise 4: Now, we are comparing optimal allocations for capital against the data:
figure(3)
plot(y,z,'o',ye,z)
title('Capital comparison')
legend('Actual capital','Efficient capital')
%%
%Exercise 6: Now, we assume now they are not independent cov(k,z)=0.5
mu=[0,0];covar=0.5;variance=1;
sigma=[variance,covar;covar,variance];
R1 = mvnrnd(mu,sigma,n); %Compute the multivariate normal distribution. Now, we do not assume independence any longer
ln_k1=R1(:,1);
ln_z1=R1(:,2);
ln_y1=(1-gamma)*ln_z1+gamma*ln_k1;
k1=exp(ln_k1);
z1=exp(ln_z1);
y1=exp(ln_y1);
ke1=z1./sum(z1).*sum(k1);
ye1=z1.^(1-gamma).*ke1.^gamma;
rate1=(sum(ye1)/sum(y1)-1).*100;
figure(4)
plot(y1,z1,'o',ye1,z1)
title('When covariance is 0.5')
%Now, we assume now they are not independent cov(k,z)=-0.5:
covar=-0.5;mu=[1,1];
sig=[variance,covar;covar,variance];
R2 = mvnrnd(mu,sig,n);
ln_k2=R2(:,1);
ln_z2=R2(:,2);
ln_y2=(1-gamma)*ln_z2+gamma*ln_k2;
k2=exp(ln_k2);
z2=exp(ln_z2);
y2=exp(ln_y2);
ke2=z2./sum(z2).*sum(k2);
ye2=z2.^(1-gamma).*ke2.^gamma;
rate2=(sum(ye2)/sum(y2)-1).*100;
figure(5)
plot(y2,z2,'o',ye2,z2)
title('When covariance is -0.5')