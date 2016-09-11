% discrete fourier transform test

fx = @(x) 4*cos(x)+cos(x).*cos(3*x);
m = 4;
t = 2*pi*(0:1/m:(m-1)/m)';

f = fx(t);
fhat = 1/m*ftransform(m)*f;
fphat = 1i*(repmat((-m/2:m/2-1)',1,size(fhat,2))).*fhat;
xp = ftransform(m)\(m*fphat);