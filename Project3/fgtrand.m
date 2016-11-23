function t = fgtrand( N)
% This is fgt for random points
%
% Source points position and Gaussian summation value at target points are stored by N/4 by 4 matrix.
% Each column corresponding to position of points inside each box.

if nargin < 1, N = 1*1e+4; end
p = 10;                     % number of fourier modes
L = 7;                      % lambda = L/(p*sqrt(delta))
delta = 1/64;                % support of Gaussian
lambda = L/(p*sqrt(delta));


%% Source point function (source strength, can be played with)
f = @(x,y) x.^2 + y.^2;


%% partition domain into uniform boxes of size sqrt(delta)
Cbx = sqrt(delta)/2:sqrt(delta):1-sqrt(delta)/2; % center of boxes based on delta
Nbox = numel(Cbx)^2;

N = Nbox*N;
%% form random points location
% Here I form random points inside each box. So I sort of already know
% which point is inside which box...
x = NaN(N/Nbox,Nbox);
y = NaN(N/Nbox,Nbox);
for i = 1:sqrt(Nbox)
    for j = 1:sqrt(Nbox)
        x(:,(i-1)*sqrt(Nbox)+j) = rand(N/Nbox,1)/2 + (i-1)*sqrt(delta); 
        y(:,(i-1)*sqrt(Nbox)+j) = rand(N/Nbox,1)/2 + (j-1)*sqrt(delta);            
    end
end


%% Some parameter
CBx = reshape(repmat(Cbx,sqrt(Nbox),1),1,Nbox);     % x-coord of center of boxes 
CBy = repmat(Cbx,1,sqrt(Nbox));                     % y-coord of center of boxes

mid = (2*p+1)*p+p+1;                                % by symmetry, only need half of s2w, w2l, l2t
freq = -p:p;                                        % frequency range
freqk = reshape(repmat(freq,2*p+1,1),1,(2*p+1)^2);  
freqk=freqk(1:mid);                                 % frequency k (only need about half)
freqs = repmat(freq,1,2*p+1); 
freqs=freqs(1:mid);                                 % frequency s ( (k,s) come in pairs in computation )


tic
%% Sources to Wave Expansions
w = NaN(mid,Nbox);                                  % store s2w for half frequency
fy = f(x,y);
for i = 1:mid
    w(i,:) = sum( fy.*exp( 1i*lambda*( freqs(i)*(ones(N/Nbox,1)*CBy-y)+freqk(i)*(ones(N/Nbox,1)*CBx-x))));
end

%% Wave to Local expansions
v = zeros(mid,Nbox);                                % store w2l for half frequency
for i = 1:mid
    for j = 1:Nbox
        v(i,j) =  sum( w(i,:).*exp(1i*lambda* (freqs(i)*(CBy(j)-CBy) + freqk(i)*(CBx(j)-CBx))));
    end
end

%% Local Expansion to Targets
Gks = (L/(2*p*sqrt(pi)))^2*exp(-lambda^2*(freqk.^2+freqs.^2)*delta/4);
F = zeros(N/Nbox,Nbox);                                   % store Gaussian summation for 4 boxes
for i = 1:N/Nbox
    for j = 1:Nbox
        F(i,j) = 2*sum( real( Gks(1:end-1)'.*v(1:end-1,j).*exp(1i*lambda*( freqs(1:end-1)'*(y(i,j)-CBy(j)) + freqk(1:end-1)'*(x(i,j)-CBx(j)))) ) )...
            +Gks(end)*v(end,j);
    end
end
t = toc




%% compare with direct sum
tnum = 3;                          % randomly pick tnum*tnum number of testing points
r1 = randi(N/Nbox, tnum,1);
r2 = randi(Nbox,tnum,1);
targetx = x(r1,r2);
targety = y(r1,r2);
Fexa = zeros(size(targetx));
[L,S] = size(Fexa);
for l = 1:L
    for s = 1:S
        for i = 1:N/Nbox
            for j = 1:Nbox
                Fexa(l,s) = Fexa(l,s) + f(x(i,j),y(i,j))*exp(-( (targetx(l,s)-x(i,j))^2+(targety(l,s)-y(i,j))^2 )/delta);
            end
        end
    end
end
err = abs((Fexa - F(r1,r2))./Fexa)
toc


end










