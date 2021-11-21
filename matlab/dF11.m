function y=dF11(k,eta,theta) 
s=@(z)1./(1+exp(-z));
integrand=@(x)1/4*x.^(k+1)./sqrt(1+theta*x/2).*s(eta-x).*(1-s(eta-x));
y=integral(integrand,0,inf,'RelTol',4.4e-16,'Waypoints',[eta]);
end