function V = ftransform( n)
% FTRANSFORM - Discrete Fourier Transform Circulant Matrix
%
% V = ftransform( n) returns discrete fourier tranform matrix
% Inputs: n = number of point in physical space
% Outputs: V = transformation matrix
% Discrete Fourier Transform matrix n dimension
% Hai 09/09/2016

V = ones(n);
w = exp(-1i*(-n/2:n/2-1)'*2*pi/n);
for k = 2:n
    V(:,k) = V(:,k-1).*w;
end

