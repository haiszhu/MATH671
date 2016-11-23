function [fx,Fx] = main(n,u,v,X,Y,Z,sigma)
% n = number of grid on uv plane, n^2 is the number of grid on surface
% u = surface target point u-coord, input should be vector
% v = surface target point v-coodr, input should be vector
% 
% X = target point x-coord in 3D cube, input should be vector
% Y = y-coord
% Z = z-coord
% sigma = force function
%
% fx = surface target point potential
% Fx = target point potential in 3D cube


%% 2d torus surface
u = u(:); v = v(:);
fx = NaN(numel(u),1);
for k = 1:numel(u)
    [x,y,z]=Torus(u(k),v(k));xprimes = [x;y;z]; 
    uvprime = [u(k),v(k)];
    [fx(k),flag] = surfint(n,xprimes,uvprime,sigma);
end

%% 3d domain
X = X(:); Y = Y(:); Z = Z(:);
Fx = NaN(numel(X),1);
for k = 1:numel(X)
    xprime = [X(k);Y(k);Z(k)];
    dist = Z(k)^2 + (sqrt(X(k)^2+Y(k)^2)-2)^2;
    
    if abs(dist-1)>1e-14,
        [Fx(k)] = globalint(n,xprime,sigma);
    else
        vx = atan(Y(k)/X(k)); 
        if X(k)<0,
            vx = vx + pi; 
        else if Y(k)<0, 
                 vx = vx + 2*pi; 
             end
        end
        ux = asin(Z(k));
        if sqrt(X(k)^2+Y(k)^2)-2 < 0
            ux = ux + pi;
        else if Z(k)<0,
                ux = ux + 2*pi;
            end
        end        
        uvprime = [ux,vx];
        [Fx(k)] = surfint(n,xprime,uvprime,sigma);
    end

end

end

%% torus geometry
function [ x, y, z] = Torus(U,V)
% 
x = (2+cos(U)).*cos(V);
y = (2+cos(U)).*sin(V);
z = sin(U);

end

