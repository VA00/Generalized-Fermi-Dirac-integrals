function y=F(k,eta,theta) 
%Calculation using the integral representation 
%A.Gil J.Segura N.M.Temme https://www.sciencedirect.com/science/article/pii/S0096300321007025
integrand=@(x)x.^k.*sqrt(1+theta*x/2)./(1+exp(x-eta));
y=integral(integrand,0,inf,'RelTol',4.4e-16,'Waypoints',[eta]);
 