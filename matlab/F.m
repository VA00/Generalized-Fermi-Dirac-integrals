function y=F(k,eta,theta) 
%Calculation using the integral representation
integrand=@(x)x.^k.*sqrt(1+theta*x/2)./(1+exp(x-eta));
y=integral(integrand,0,inf,'RelTol',4.4e-16,'Waypoints',[eta]);
 