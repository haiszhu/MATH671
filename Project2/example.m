function [fx,Fx] = example()
%% This function mainly serves as an example of computing surface integration on Torus and 3D domain
% test over surface and 3-d domain
% the following would be examples of input format

% input of surface integration
n = 2^7; u=[0,pi/3,pi]; v=[pi/3,0,pi];

% input should be a vector contains x, y, z coord
[X,Y,Z] = Torus(u,v);
% X=X+0.01; Y=Y+0.01; Z=Z+0.01;

% input should be force function format
sigma = @(x,y,z) ones(size(x));


[fx, Fx] = main(n,u,v,X,Y,Z,sigma);

end

%% torus geometry
function [ x, y, z] = Torus(U,V)
% 
x = (2+cos(U)).*cos(V);
y = (2+cos(U)).*sin(V);
z = sin(U);

end
