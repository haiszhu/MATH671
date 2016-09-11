function V = iftransform( n)
% FTRANSFORM - Discrete Inverse Fourier Transform Circulant Matrix
%
% V = ftransform( n) returns discrete inverse fourier tranform matrix
% Inputs: n = number of point in fourier space
% Outputs: V = transformation matrix
% Discrete Inverse Fourier Transform matrix n dimension
% Hai 09/10/2016

V = ones(n);
w = exp(1i*(-n/2:n/2-1)*2*pi/n);
for k = 2:n
    V(k,:) = V(k-1,:).*w;
end