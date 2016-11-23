function mainvcycle()
%
% function mainvcycle()
%
% multigrid for -u_xx-u_yy = f on [0,1]x[0,1], problem setup and
% initialization. The coarse grid finite difference 5-point stencil is
% generated directly.



N = 64;                % grid size (power of 2)
alpha = 1; beta = 1;    % ratio between dx dy and h
vcycles = 15;           % number of V-cycles
smoothnum = 3;          % number of smoothing iteration, default value is 5 if no value specified
matrix = 'interp';      % matrix = 'direct', default is 'direct' if no method is specified

%% generate mesh
h = 1/N; dx = alpha*h; dy = beta*h;
x = 0:dx:1; y = 0:dy:1;
[X,Y] = meshgrid(x,y);
[ny,nx] = size(X);

%% get exact reference solution and rhs
uexact = ufun(X,Y);
f = ffun(X,Y);

%% set up boundary data
u = zeros(ny,nx);
u(1:ny, 1) = uexact(1:ny, 1);
u(1:ny,nx) = uexact(1:ny,nx);
u( 1,1:nx) = uexact( 1,1:nx);
u(ny,1:nx) = uexact(ny,1:nx);

%% finite difference matrix
A = fd5point(nx,ny,alpha,beta);

%% multigrid V-cycle
rnorm = zeros(vcycles+1,1);
rnorm(1) = norm(residual(A,f,u,h))*h;
enorm = zeros(vcycles+1,1);
enorm(1) = norm(u-uexact)*h;
for i=2:vcycles+1 

	u = mgv(A, f, u, h, alpha, beta, smoothnum, matrix );
    rnorm(i)  = norm( residual(A,f,u,h) )*h;
    enorm(i)  = norm( u-uexact )*h;
	
	disp(['cycle ', num2str(i-1), ' residual norm = ', num2str(rnorm(i))]);

end

% plot residual vs. number of vcycles
subplot(1,2,1);
semilogy(0:vcycles,rnorm,'r*-');
xlabel('Vcycle');
ylabel('Residual Norm');
title('Multigrid Method');

% plot solution error vs. number of sweeps
subplot(1,2,2);
semilogy(0:vcycles,enorm,'r*-');
xlabel('Vcycle');
ylabel('Error Norm');
title('Multigrid Method');


end

function u = mgv(A, f, u, h, alpha, beta, smoothnum, matrix)
%
% function u = mgv(A, f, u, h, smoothnum, alpha, beta)
%
% Multigrid V-cycle (recursively call itself)
%
% Input
%   A         :    finite difference matrix for 2d laplace 
%               (everything has been taken care of, except h^2)
%   f         :    - u_xx - u_yy = f (in matrix form)
%   u         :    - u_xx - u_yy = f (in matrix form)
%   h         :    mesh size
%   alpha     :    ratio between dx and h
%   beta      :    ratio between dy and h
%   smoothnum :    number of smoothing steps
%   matrix    :    choose to form fd matrix directly or by interpolation
%
% Output
%   u         :    solution to A*u = f after one V-cycle
%
% 
%  11/19/2016 Hai Zhu
%

% intialize
if nargin < 8, matrix = 'direct'; end
if nargin < 7, smoothnum = 5; end

[ny,nx] = size(f);

% set up coarse grid
if nx == 3 || ny == 3
    temp = A\(h^2*reshape(f(2:ny-1,2:nx-1),(nx-2)*(ny-2),1));
    u(2:ny-1,2:nx-1) = reshape(temp,ny-2,nx-2); 
else
    %% smoothing
    u = djrelax(A, f, u, h, smoothnum );
    
    
    %% residual correction scheme
    res = residual(A, f, u, h );
    
    % restriction of residual
    nxc = (nx-1)/2 + 1;
    nyc = (ny-1)/2 + 1;
    fc = zeros(nyc,nxc);
    temp = restrict_mat(nx,ny)*reshape(res(2:end-1,2:end-1),(nx-2)*(ny-2),1);
    fc(2:end-1,2:end-1) = reshape(temp,nyc-2,nxc-2);
    
    uc = zeros(nyc,nxc);    % initial guess on coarse grid
    hc = 2*h;
    
    switch matrix
        
        case 'direct'
            Ac = fd5point(nxc,nyc,alpha,beta);
            
        case 'interp'
            Ac = fd5pointc(A,nx,ny);
     
    end
    
    uc = mgv(Ac, fc, uc, hc, alpha, beta, smoothnum, matrix);
    
    %% prolong
    u = u + c2f( uc );
    
    %% post smoothing
    u = djrelax(A, f, u, h , smoothnum);
      
end


end



function Ac = fd5pointc(A,nx,ny)
%
% function Ac = fd5pointc(A,nx,ny)
%
% Finite difference scheme 5 point stencil for coarse grid
%
% Input
%   A          :    coefficient matrix for finite difference on fine grid
%   nx         :    number of grid along x direction
%   ny         :    number of grid along y direction
%       (the reason we require nx and ny here is because A is always a square matrix, 
%       it is a (nx-2)*(ny-2) by (nx-2)*(ny-2) matrix)
%
% Output
%   Ac         :    coefficient matrix on coarse grid
% 
%  11/20/2016 Hai Zhu
%

Iinject = inject_mat(nx,ny);                % get injection matrix
Iprolon = (4*restrict_mat(nx,ny))';         % get prolongation matrix

Ac = Iinject*2*A*Iprolon;                   % finally form finite difference matrix on coarse grid 

end


function Iinject = inject_mat(nx,ny)
%
% function Iinject = inject_mat(nx,ny)
%
% Prolongation matrix
%
% Input
%   nx          :    number of grid along x direction
%   ny          :    number of grid along y direction
%
% Output
%   Iinject     :    direct injection matrix
%
% 
%  11/20/2016 Hai Zhu
%

Nx = nx-2; Ny = ny-2;           % number of points inside computational domain

onevecy = ones(Ny,1);          
T = spdiags(onevecy, 0, Ny,Ny);
Imat = T(2:2:Ny-1,:);

onevecx = ones(Nx,1);
Ix = spdiags(onevecx, 0,Nx,Nx);
Ix = Ix(2:2:Nx-1,:);

Iinject = kron(Ix,Imat);

end



function A = fd5point(nx,ny,alpha,beta)
%
% function A = fd5point(nx,ny,alpha,beta)
%
% Finite difference scheme 5 point stencil
%
% Input
%   nx         :    number of grid along x direction
%   ny         :    number of grid along y direction
%   alpha      :    ratio between dx and h
%   beta       :    ratio between dy and h
%
% Output
%   A          :    coefficient matrix without effect of h
% 
%  11/19/2016 Hai Zhu
%

%% locally affected by ajacent rowwise elements (horizontally)
Nx = nx - 2; 
Ix = speye(Nx);
ex = ones(Nx,1)/alpha^2;
Tx = spdiags([ex -2*ex ex],[-1 0 1],Nx,Nx);

%% locally affected by ajacent columnwise elements (vertically)
Ny = ny-2;
Iy = speye(Ny);
ey = ones(Ny,1)/beta^2;
Ty = spdiags([ey -2*ey ey],[-1 0 1],Ny,Ny);

%% kron product to combine ajacent 4 points
A = -(kron(Ix,Ty) + kron(Tx,Iy));

end


function u = c2f( uc )
%
% function uc = c2f( u )
%
% Prolongation
%
% Input
%   uc         :    solution on coarse grid
%
% Output
%   u          :    solution on fine grid by prolongation
%
% 
%  11/19/2016 Hai Zhu
%

%% get prolongation matrix
[ny,nx] = size( uc );
Ih2h = prolong_mat(nx,ny);

%% assemble solution on fine grid
nxf = 2*(nx-1)+1;  nyf = 2*(ny-1)+1;        % grid size on fine grid
u = zeros(nyf,nxf);
temp = Ih2h*reshape(uc(2:ny-1,2:nx-1),(nx-2)*(ny-2),1);
u(2:end-1,2:end-1) = reshape(temp,nyf-2,nxf-2);


end

function Ih2h = prolong_mat(nx,ny)
%
% function Ih2h = prolong_mat(nx,ny)
%
% Prolongation matrix
%
% Input
%   nx          :    number of grid along x direction
%   ny          :    number of grid along y direction
%
% Output
%   Ih2h        :    I^h_{2h} by interpolation
%
% 
%  11/19/2016 Hai Zhu
%

nxf = 2*(nx-1)+1;       % fine grid size
nyf = 2*(ny-1)+1;
Nx = nxf-2;
Ny = nyf-2;

%% form subblocks of restriction ( weighted u for each column ) 
onevecy = ones(Ny,1);
T = spdiags([1/16*onevecy 1/8*onevecy 1/16*onevecy], [-1,0,1], Ny,Ny);
Imat = T(2:2:Ny-1,:);

%% form global pattern restriction matrix ( weighted u for all columns )
onevecx = ones(Nx,1);
Ix = spdiags([onevecx 2*onevecx onevecx], [-1,0,1],Nx,Nx);
Ix = Ix(2:2:Nx-1,:);

%% combine global pattern with local subblock
Ih2h = 4*kron(Ix,Imat)';

end


function uc = f2c( u )
%
% function uc = f2c( u )
%
% Restriction
%
% Input
%   u          :    solution on fine grid
%
% Output
%   uc         :    solution on coarse grid by restriction
%
% 
%  11/19/2016 Hai Zhu
%

%% get restriction matrix
[nx,ny] = size(u);                              % fine grid size
I2hh = restrict_mat(nx,ny);                     % restriction matrix


%% assemble solution on coarse grid
nxc = (nx-1)/2+1; nyc = (nx-1)/2+1;             % coarse grid size
uc = zeros(nyc,nxc);
temp = I2hh*reshape(u(2:ny-1,2:nx-1),(nx-2)*(ny-2),1);
uc(2:end-1,2:end-1) = reshape(temp,nyc-2,nxc-2);

end

function I2hh = restrict_mat(nx,ny)
%
% function I2hh = restrict_mat(nx,ny)
%
% Restriction matrix
%
% Input
%   nx          :    number of grid along x direction
%   ny          :    number of grid along y direction
%
% Output
%   I2hh        :    I^{2h}_h by taking full weight
%
% 
%  11/19/2016 Hai Zhu
%

Nx = nx-2; Ny = ny-2;           % number of points inside computational domain

%% form subblocks of restriction ( weighted u for each column ) 
onevecy = ones(Ny,1);           % each column has Ny grids inside domian
T = spdiags([1/16*onevecy 1/8*onevecy 1/16*onevecy], [-1,0,1], Ny,Ny);
Imat = T(2:2:Ny-1,:);

%% form global pattern restriction matrix ( weighted u for all columns )
onevecx = ones(Nx,1);
Ix = spdiags([onevecx 2*onevecx onevecx], [-1,0,1],Nx,Nx);
Ix = Ix(2:2:Nx-1,:);

%% combine global pattern with local subblock
I2hh = kron(Ix,Imat);

end


function u = djrelax(A, f, u, h , smoothnum )
%
% function u = djrelax(A, f, u, h )
%
% Relaxation using damped Jacobi iteration
%
% Input
%   A           :    finite difference matrix for 2d laplace 
%                   (everything has been taken care of, except h^2)
%   f           :    - u_xx - u_yy = f (in matrix form)
%   u           :    - u_xx - u_yy = f (in matrix form)
%   h           :    mesh size
%   smoothnum   :    number of relaxation
%
% Output
%   u           :    iterative solution after smoothnum number of iteration
%
% 
%  11/19/2016 Hai Zhu
%

%% factorize matrix A
[ny,nx] = size(f);          % computational domain size
h2 = h*h;                   % A needs to take in account of h^2 
Nx = nx - 2; Ny = ny - 2;   % number of points inside computational domain

D = diag(diag(A));          % sparse diagonal of A
U = -triu(A,1);             % A = D - L -U
L = -tril(A,-1);            % sparse lower triagnular matrix


%% actual smoothing and assembly
temp = reshape(u(2:end-1,2:end-1),Nx*Ny,1);     % all u inside computational domain 
ftemp = reshape(f(2:end-1,2:end-1),Nx*Ny,1);    % all rhs inside computational domain
w = 4/5;                                        % chosen to minimize high frequency eignvalue
for i = 1:smoothnum                             
    temp = ((1-w)*speye(Nx*Ny)+w*(D\(L+U)))*temp + w*(D\(h2*ftemp));
end
u(2:end-1,2:end-1) = reshape(temp,Ny,Nx);       % reshape to original matrix format

end


function res = residual( A, f, u, h )
%
% function res = residual( A, f, u, h )
%
% Compute residual r = b - A*u
%
% Input
%   A           :    finite difference matrix for 2d laplace 
%                   (everything has been taken care of, except h^2)
%   f           :    - u_xx - u_yy = f (in matrix form)
%   u           :    - u_xx - u_yy = f (in matrix form)
%   h           :    mesh size
%
% Output
%   res         :    residual r = b - A*u
%
% 
%  11/19/2016 Hai Zhu
%

[ny,nx] = size(f);                              % computational domain size
h2 = h*h;                                       % A needs to take in account of h^2 
Nx = nx - 2;                                    % number of points inside computational domain along x direction
Ny = ny - 2;                                    % number of points inside computational domain along y direction
res = zeros(ny,nx);                             % initialize residual on grid

temp = reshape(u(2:end-1,2:end-1),Nx*Ny,1);     % all u inside computational domain
ftemp = reshape(f(2:end-1,2:end-1),Nx*Ny,1);    % all rhs inside computational domain 
temp = reshape(ftemp - A*temp/h2,Ny,Nx);        % residual inside computational domain
res(2:end-1,2:end-1) = temp;

end

function u = ufun( X, Y )
%
% function u = ufun( X, Y )
%
% Compute exact solution for 2d laplace on grid to compare with
%
% Input
%   X   :    x-coordinates 2d format generated by meshgrid
%   Y   :    y-coordinates 
%
% Output
%   u   :    exact solution on 
%
% 
%  11/19/2016 Hai Zhu
%

u = (X.^2-X.^4).*(Y.^4-Y.^2);

end 

function f = ffun( X, Y )
%
% function f = ffun( X, Y )
%
% Compute rhs for 2d laplace ( - u_xx - u_yy = f ) on grid 
%
% Input
%   X   :    x-coordinates 2d format generated by meshgrid
%   Y   :    y-coordinates 
%
% Output
%   f   :    right hand side of Au = f
%
% 
%  11/19/2016 Hai Zhu
%

f = 2*( (1-6*X.^2).*Y.^2.*(1-Y.^2) + (1-6*Y.^2).*X.^2.*(1-X.^2) );

end



