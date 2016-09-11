function myloglog(h,E,s)
%
% Produce loglog plot of CPU time vs. Number of points.
% Estimate order of accuracy by doing a linear least squares fit.
%
% Based on http://www.amath.washington.edu/~rjl/fdmbook/  (2007)
% Hai 09/11/16

h = h(:);            % make sure it's a column vector
E = E(:);            % make sure it's a column vector
ntest = length(h);
clf
loglog(h,E,'o-')
axis([.5*min(h) 1.5*max(h)  .5*min(E) 1.5*max(E)])
if s == 'c'
    title('log-log plot of CPU time vs. M')
else
    title('log-log plot of h vs. error for Euler')
end

% Estimate order of accuracy from least squares fit:
Ap = ones(ntest,2);
Ap(:,2) = log(h);
bp = log(E);
Kp = Ap\bp;
K = Kp(1);
p = Kp(2);
if s == 'c' 
    disp(' ')
    disp(sprintf('Least squares fit gives T(M) = %g * M^%g',exp(K),p))
    disp(' ')
else
    disp(' ')
    disp(sprintf('Least squares fit gives E(h) = %g * h^%g',exp(K),p))
    disp(' ')
end

% add graph of this line to loglog plot:
hold on
err1 = exp(K)*h.^p;
loglog(h,err1,'r')
legend('errors', 'least squares fit','Location','SouthEast')
hold off