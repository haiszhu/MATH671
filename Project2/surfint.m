function [fxprime, flag] = surfint( numg, xprime, uvprime, sigma)
close all
flag = 0;
[ xprimex, xprimey, xprimez] = Torus(uvprime(1),uvprime(2)); xprime2 = [xprimex;xprimey;xprimez];
if norm(xprime(:)-xprime2) > 1e-14, flag = 1; end
h = 2*pi/numg;
[U,V] = meshgrid( 0:h:2*pi, 0:h:2*pi);

[ xglobalx, xglobaly, xglobalz] = Torus(U,V);

W = 2 + cos(U); % element area on torus, analytic

%% kernel function
 kernel = @(x,y,z) 1./(sqrt( (x-xprime(1)).^2 + (y-xprime(2)).^2 + (z-xprime(3)).^2)); 
 
%%%%%%%%%%%%%%%%%%%%% local integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 1 of local mollification: nearby points index computation 
% [ idxu, idxv, numrl, i, j, uvprimer] = lindex( uvprime, h, numg);
[ numrl, i, j, uvprimer] = lindex( uvprime, h);
ulocalr = uvprimer(1)-numrl*h : h : uvprimer(1)+numrl*h;
vlocalr = uvprimer(2)-numrl*h : h : uvprimer(2)+numrl*h;
[ Ulr, Vlr] = meshgrid( ulocalr, vlocalr);

%% step 2 of local mollification: local integration over u-v plane on regular grid 
[ xlocalrx, xlocalry, xlocalrz] = Torus(Ulr, Vlr);
wlocalr = 2 + cos(Ulr);
% F = kernel(xlocalrx, xlocalry, xlocalrz);
forcelr = sigma(xlocalrx,xlocalry,xlocalrz);
distlr = sqrt( (xlocalrx-xprime(1)).^2 + (xlocalry-xprime(2)).^2 + (xlocalrz-xprime(3)).^2);
if norm(xprime - [xlocalrx(numrl+1,numrl+1);xlocalry(numrl+1,numrl+1);xlocalrz(numrl+1,numrl+1)])<1e-14, xlocalrx(numrl+1,numrl+1) = xprime(1)+h; end   
[ fxprimelr] = intlocalr( distlr, numrl, xlocalrx, xlocalry, xlocalrz, wlocalr, h, kernel, forcelr);

%% step 3 of local mollification: local integration to add to global integration (on polar grid)
% step 4 of local mollification: substep 1 -- set up local nodes in polar grid
[ rp, thetap, dr, dtheta, ulocalp, vlocalp] = DistMatrixp( uvprime, h); 
% step 4 of local mollification: substep 2 -- local integration over on polar grid
[xlocalpx, xlocalpy, xlocalpz] = Torus(ulocalp,vlocalp); wlocalp = 2 + cos(ulocalp);
% F = kernel(xlocalpx, xlocalpy, xlocalpz);
forcelp = sigma(xlocalpx,xlocalpy,xlocalpz);
[ fxprimelp] = intlocalp( xprime, xlocalpx, xlocalpy, xlocalpz, wlocalp, rp, h, dr, dtheta, kernel, forcelp);


%%%%%%%%%%%%%%%%%%%%% global integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% global integration over u-v plane on regular grid (check with exponential convergence rate)
forceg = sigma(xglobalx,xglobaly,xglobalz);
if norm(xprime - [xglobalx(j,i);xglobaly(j,i);xglobalz(j,i)])<1e-14, xglobalx(j,i) = xprime(1)+h; end   
% just to avoid the exception of x = x', doesn't change overall algorithm
[fxprimegr] = intglobalr( xglobalx, xglobaly, xglobalz, W, h, kernel, forceg);   % fx' on global regular grid

%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fxprimer = fxprimegr - fxprimelr;
fxprime = fxprimer + fxprimelp;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Integral evaluation on polar grid locally
function fxprimelp = intlocalp( xprime, xlocalpx, xlocalpy, xlocalpz, wlocalp, rp, h, dr, dtheta, fTest, force)

% Hai 10/22/16

% mollifier function
dist = sqrt( (xlocalpx-xprime(1)).^2 + (xlocalpy-xprime(2)).^2 + (xlocalpz-xprime(3)).^2);
psi = mypsi(dist/sqrt(h));

F = fTest(xlocalpx(:,2:end), xlocalpy(:,2:end), xlocalpz(:,2:end)).*rp(:,2:end);
F = [ones(size(F(:,1))),F];
integrand = F.*wlocalp.*psi;
fxprimelp = dr*dtheta*( sum(sum(integrand(2:end-1,2:end-1))) + 1/2*sum(integrand(1,2:end-1)) + 1/2*sum(integrand(end,2:end-1))...
    + 1/2*sum(integrand(2:end-1,1)) + 1/2*sum(integrand(2:end-1,end)) ...
    + 1/4*(integrand(1,1)+integrand(1,end)+integrand(end,1)+integrand(end,end)));

end


%% set up local integration nodes in polar coordinates
function [ rp, thetap, dr, dtheta, ulocalp, vlocalp] = DistMatrixp(uvprime, h)

nump = 4*round(ceil(sqrt(h)/h));  
dr = 2*sqrt(h)/nump;
dtheta = 2*pi/nump;   

r = 0:dr:2*sqrt(h);  
theta = 0:dtheta:2*pi;   
[rp,thetap] = meshgrid(r,theta);  
ulocalp = uvprime(1)+rp.*cos(thetap);
vlocalp = uvprime(2)+rp.*sin(thetap);

end


%% Integral Evaluation on regular grid globally (periodic)
function [ fxprimegr] = intglobalr( x, y, z, W, h, kernel, force)
%Return the global part for the integration using trapezoidal rule
%   fTest: Density function on surface 
%   xprime: target point on surface
%   x,y,z: coordinates on surfac
%   W: element area on surface
%   h: grid size on u-v plane
%
%   10/22/2016 Hai

%%  
% x(j,i) = x(j,i) + h;    % this is also bad in global integration
F = kernel(x, y, z);
integrand = W.*F.*force;
fxprimegr= h^2*sum(sum(integrand(1:end-1,1:end-1)));

% dist = sqrt( (x-xprime(1)).^2 + (y-xprime(2)).^2 + (z-xprime(3)).^2);
% integrandpsi = integrand.*(1-mypsi(dist/sqrt(h)));
% fxpsi = h^2*sum(sum(integrandpsi(1:end-1,1:end-1)));

end

%% Integral Evaluation on regular grid locally 
function fxprimelr = intlocalr( dist, numrl, x, y, z, W, h, kernel, force)
%Return the local part for the integration using trapezoidal rule
%   fTest: Density function on surface 
%   xprime: target point on surface
%   x,y,z: coordinates on surfac
%   W: element area on surface
%   h: grid size on u-v plane
%
%   10/22/2016 Hai

% mollifier function
psi = mypsi(dist/sqrt(h));
F = kernel(x, y, z);
integrand = W.*F.*psi.*force;
% integrand(numrl+1,numrl+1) = 0;
fxprimelr= h^2*sum(sum(integrand));

end


%% local index computation
function [ numrl, i, j, uvprimer] = lindex( uvprime, h)
% numg = global mesh num
% numrl = local mesh num on regular grid
% uvprimer = round version of uvprime
idxuvp = round(uvprime/h)+1;
i = idxuvp(1); j = idxuvp(2); 
uvprimer = [h*(i-1), h*(j-1)];
numrl = 2*ceil(sqrt(h)/h);
% idxu = mod( (i-numrl:i+numrl)-1, numg)+1;
% idxv = mod( (j-numrl:j+numrl)-1, numg)+1;
end


%% torus geometry
function [ x, y, z] = Torus(U,V)
% 
x = (2+cos(U)).*cos(V);
y = (2+cos(U)).*sin(V);
z = sin(U);

end

%% mollifier function
function Psi = mypsi(dist)

Psi = NaN(size(dist));

H = @(r) exp(-2*exp(-1./r)./(1-r));
Psi = H(dist);
Psi(dist>=1)=0;

end
