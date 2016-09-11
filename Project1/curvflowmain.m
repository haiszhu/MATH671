% main function to test project curvature flow project
% Hai 09/10/2016

clear all
close all

%% Initial geometry set up
s = @(t) [(4+cos(3*t)).*cos(t), (4+cos(3*t)).*sin(t)];

%% Part 1 curve shape evolution
N = 400;    % time discretization
M = 32;     % space discretization
tfinal = 4; % final evolution time

tt = 0:tfinal;
count = numel(tt)-1;
xx = 2*pi*(0:1/M:(M-1)/M)';
x = NaN( M, 2, count);
x(:,:,1) = s(xx);
figure()
    subplot(3,2,1)
    plot(real(x(:,1,1)),real(x(:,2,1)))
    title(sprintf('time t = %f',tt(1)))
for k = 1:count 
    x(:,:,k+1) = curvflow( N/count, M, x(:,:,k), tt(k), tt(k+1));  
    subplot(3,2,k+1)
    plot(real(x(:,1,k+1)),real(x(:,2,k+1)))
    title(sprintf('time t = %f',tt(k+1)))
end

%% Part 2 CPU time per time-step
M = 2.^(7:12);
tcpu = NaN(numel(M),1);
for  k=1:numel(M)
   xx = 2*pi*(0:1/M(k):(M(k)-1)/M(k))';
   x = s(xx);
   tic
   x = curvflow( 1, M(k), x, 0, tfinal/N); 
   tcpu(k) = toc;
end
figure()
myloglog(M/M(1),tcpu,'c')   % M has been rescaled by M(1)
mytable(M,tcpu,'c')

%% Part 3 Order of accuracy
N = 2.^(2:6)*100;
M = 64;     % space discretization
tfinal = 4; % final evolution time

count = numel(N);
xx = 2*pi*(0:1/M:(M-1)/M)';
L = NaN(count,1);
for k = 1:count
    [x,xp] = curvflow( N(k), M, s(xx), 0, tfinal);
    L(k) = sum(sqrt(diag((xp*xp'))));
end
h = 4./N(1:end-1); E = L-L(end); E = abs(E(1:end-1));
figure()
myloglog(h,E,'e')
mytable(h,E,'e')



