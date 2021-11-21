function y=dF03(k,eta,theta) 
s=@(z)1./(1+exp(-z));
integrand=@(x)3/64*x.^(k+3)./(sqrt(1+theta*x/2).^5).*s(eta-x);
y=integral(integrand,0,inf,'RelTol',4.4e-16,'Waypoints',[eta]);
end