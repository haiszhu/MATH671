function [ fxprime] = globalint( numg, xprime, sigma)
%Return the global part for the integration using trapezoidal rule
%   fTest: Density function on surface 
%   xprime: target point on surface
%   x,y,z: coordinates on surfac
%   W: element area on surface
%   h: grid size on u-v plane
%
%   10/29/2016 Hai

%%  
h = 2*pi/numg;
[U,V] = meshgrid( 0:h:2*pi, 0:h:2*pi);
[ xx, yy, zz] = Torus(U,V);
W = 2 + cos(U);

kernel = @(x,y,z) 1./(sqrt( (x-xprime(1)).^2 + (y-xprime(2)).^2 + (z-xprime(3)).^2)); 

force = sigma(xx,yy,zz);
F = kernel(xx, yy, zz);
integrand = W.*F.*force;
fxprime= h^2*sum(sum(integrand(1:end-1,1:end-1)));

end

%% torus geometry
function [ x, y, z] = Torus(U,V)
% 
x = (2+cos(U)).*cos(V);
y = (2+cos(U)).*sin(V);
z = sin(U);

end
