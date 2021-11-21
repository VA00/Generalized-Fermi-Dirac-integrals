function y=dF10(k,eta,theta) 
s=@(z)1./(1+exp(-z));
integrand=@(x)x.^k.*sqrt(1+theta*x/2).*s(eta-x).*(1-s(eta-x));
y=integral(integrand,0,inf,'RelTol',4.4e-16,'Waypoints',[eta]);
end