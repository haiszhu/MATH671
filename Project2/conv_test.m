% singular integration convergence test
N = 2.^(4:11);
val = NaN(size(N));
fx = val;
fxr = val;
fxpsi = val;
for k = 1:length(N)
    n = N(k);
    uvprime = [2*pi-0.1,2*pi-0.1];
    [fx(k),val(k),fxr(k)]  = singular(n,uvprime);
end
err = (val - val(end));
error_table(2*pi./N(1:end-1),abs(err(1:end-1)))
err = (fx - fx(end));
error_table(2*pi./N(1:end-1),abs(err(1:end-1)))

