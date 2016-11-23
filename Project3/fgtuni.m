function fgtuni()
% This is fgt for uniformly located grid point. x = h/2:h:1-h/2, y = h/2:h:1-h/2;
%
%   Source points position and Gaussian summation value at target points are stored in sqrt(N)*sqrt(N) matrix 
%   And no source strength function f(x)

N = 1*1e+4;         % number of points (both source and target)
p = 10;             % number of fourier modes
L = 7;              % lambda = L/(p*sqrt(delta))
delta = 1/4;        % support of Gaussian
lambda = L/(p*sqrt(delta));

h = 1/sqrt(N);      % mesh along each dimension
x = h/2:h:1-h/2;    % form uniform located grid points (source and targes)


%% partition domain into uniform boxes of size sqrt(delta)
Cbx = sqrt(delta)/2:sqrt(delta):1-sqrt(delta)/2; % center of boxes based on delta
Cdx = sqrt(delta)/2:sqrt(delta):1-sqrt(delta)/2;
Nbox = numel(Cbx)^2;                             % number of boxes

%% Sources to Wave Expansions
% loop over frequency k, s for each box
w = NaN(2*p+1,2*p+1,Nbox);
for k = -p:p
    for s = -p:p
        for i = 1:sqrt(Nbox)
            % i th column of boxes
            for j = 1:sqrt(Nbox)
                % j th row of boxes
                nindex = [(i-1),i]*sqrt(N)/sqrt(Nbox);      % source points in box ith column, jth row
                mindex = [(j-1),j]*sqrt(N)/sqrt(Nbox);         
                w(k+p+1,s+p+1,(i-1)*sqrt(Nbox)+j) = ...
                    sum(sum(exp(1i*lambda*s*(Cbx(j)-x(mindex(1)+1:mindex(2)))).'*...
                    exp(1i*lambda*k*(Cbx(i)-x(nindex(1)+1:nindex(2))))));
            end
        end
    end
end

%% Wave to Local expansions (could use some optimization)
% compute W2L for each target box (i,j), each frequency (k,s), influenced
% by each source box (si,sj)
v = zeros(2*p+1,2*p+1,Nbox);
for i = 1:sqrt(Nbox)
    % i th column of boxes
    for j = 1:sqrt(Nbox)
        % j th row of boxes
        % Cdx(i)-- x coord, Cdx(j)-- y coord
        for k = -p:p
            for s = -p:p
                for si = 1:sqrt(Nbox)
                    % source box
                    for sj = 1:sqrt(Nbox)
                        v(k+p+1,s+p+1,(i-1)*sqrt(Nbox)+j) = v(k+p+1,s+p+1,(i-1)*sqrt(Nbox)+j)+...
                            w(k+p+1,s+p+1,(si-1)*sqrt(Nbox)+sj)*...
                            exp(1i*lambda*( k*(Cdx(i)-Cbx(si)) + s*(Cdx(j)-Cbx(sj)) ));
                    end
                end
            end
        end
    end
end

%% Local Expansion to Targets
[kk,ss] = meshgrid(-p:p,-p:p);
Gks = (L/(2*p*sqrt(pi)))^2*exp(-lambda^2*(kk.^2+ss.^2)*delta/4);    % G_k hat coefficient that depends on frequency 

F = zeros(sqrt(N),sqrt(N));
for i = 1:sqrt(Nbox)
    % i th column of boxes for calculation of wk
    for j = 1:sqrt(Nbox)
        % j th row of boxes for calculation of wk
        nindex = [(i-1),i]*sqrt(N)/sqrt(Nbox);
        mindex = [(j-1),j]*sqrt(N)/sqrt(Nbox);       
        % sum over kernel potential of |k|<p frequency
        for k = -p:p
            for s = -p:p
                F(mindex(1)+1:mindex(2),nindex(1)+1:nindex(2)) = F(mindex(1)+1:mindex(2),nindex(1)+1:nindex(2))+...
                    Gks(k+p+1,s+p+1)*v(k+p+1,s+p+1,(i-1)*sqrt(Nbox)+j)*...   
                    exp(1i*lambda*s*(x(mindex(1)+1:mindex(2))-Cbx(j))).'*...
                    exp(1i*lambda*k*(x(nindex(1)+1:nindex(2))-Cbx(i)));
            end
        end        
    end
end




%% testing
%target = [0.625,0.875];
r = randi([1 sqrt(N)], 1, 2);   % randomly pick a point to compute Gaussian sum directly
target = [x(r(1)),x(r(2))];
Fexa = 0;
for n = 1:sqrt(N)
   for m = 1:sqrt(N)
        Fexa = Fexa + exp(-( (target(1)-x(n))^2 + (target(2)-x(m))^2 )/delta);
   end
end
err = (Fexa - F(r(2),r(1)))/Fexa

end


