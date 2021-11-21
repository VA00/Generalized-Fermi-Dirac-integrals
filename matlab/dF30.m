function y=dF30(k,eta,theta) 
s=@(z)1./(1+exp(-z));
integrand=@(x)x.^k.*sqrt(1+theta*x/2).*s(eta-x).*(1-s(eta-x)).*(1-6*s(eta-x)+6*s(eta-x).^2);
y=integral(integrand,0,inf,'RelTol',1.1e-16,'Waypoints',[eta]);
end