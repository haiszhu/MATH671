% main function to test project curvature flow project
% Hai 09/10/2016

clear all
close all

%% Initial geometry set up
s = @(t) [(4+cos(3*t)).*cos(t), (4+cos(3*t)).*sin(t)];

%% Part 1 curve shape evolution
N = 400;    % time discretization
M = 256;     % space discretization
tfinal = 4; % final evolution time

tt = 0:tfinal;
count = numel(tt)-1;
xx = 2*pi*(0:1/M:(M-1)/M)';
x = NaN( M, 2, count);
x(:,:,1) = s(xx);
figure()
    plot(real(x(:,1,1)),real(x(:,2,1)))
for k = 1:count 
    x(:,:,k+1) = curvflow( N/count, M, x(:,:,k), tt(k), tt(k+1));  
    figure()
    plot(real(x(:,1,k+1)),real(x(:,2,k+1)))
end
