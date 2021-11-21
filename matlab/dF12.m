function y=dF12(k,eta,theta) 
s=@(z)1./(1+exp(-z));
integrand=@(x)-1/16*x.^(k+2)./(sqrt(1+theta*x/2).^3).*s(eta-x).*(1-s(eta-x));
y=integral(integrand,0,inf,'RelTol',4.4e-16,'Waypoints',[eta]);
end