%% ======================================================================
%  POISSON EQUATION  --  JACOBI ITERATIVE METHOD
%  ----------------------------------------------------------------------
%  Solves   u_xx + u_yy = f(x,y)   on the rectangle  [ax,bx] x [ay,by]
%  with the source term f AND the Dirichlet boundary conditions entered
%  DIRECTLY as function handles in the USER INPUT block. No hand
%  discretisation and no manual calculation are required.
%
%  SIGN CONVENTION: the code solves  u_xx + u_yy = f(x,y).
%  If your problem is written as  -(u_xx + u_yy) = g(x,y), simply enter
%  f = @(x,y) -g(x,y). Nothing else changes.
%
%  Method: 5-point finite differences on an Nx-by-Ny grid.
%  Jacobi update (ALL right-hand-side values taken from iteration k):
%
%     u(i,j)^(k+1) = [ cx*( u(i+1,j)^(k) + u(i-1,j)^(k) )
%                    + cy*( u(i,j+1)^(k) + u(i,j-1)^(k) )
%                    - f(x_i, y_j) ] / ( 2*(cx+cy) )
%
%     cx = 1/hx^2 ,   cy = 1/hy^2
%
%  Default example: f = -2*pi^2*sin(pi*x)*sin(pi*y) with u = 0 on the
%  whole boundary of the unit square. Exact solution (for checking):
%     u(x,y) = sin(pi*x) * sin(pi*y)
%
%  Requires MATLAB R2016b or newer (local function inside a script).
%  ======================================================================
clear; clc; close all;

%% --------------------- USER INPUT (edit only here) -------------------
ax = 0;   bx = 1;                 % x-interval [ax, bx]
ay = 0;   by = 1;                 % y-interval [ay, by]
Nx = 41;                          % grid points in x (>= 3)
Ny = 41;                          % grid points in y (>= 3)

% Source term of  u_xx + u_yy = 8x^2y^2. Use element-wise operators
% (.*, ./, .^). A constant such as f = @(x,y) 4 also works.
f = @(x,y) -2*pi^2 * sin(pi*x) .* sin(pi*y);

% Dirichlet boundary conditions, written exactly as given in the problem.
g_left   = @(y) 0*y;              % u(ax, y)
g_right  = @(y) 0*y;              % u(bx, y)
g_bottom = @(x) 0*x;              % u(x, ay)
g_top    = @(x) 0*x;              % u(x, by)

tol        = 1e-8;                % stop when max|change| < tol
maxIter    = 50000;               % iteration cap (choose any value)
printEvery = 500;                 % progress print interval (0 = silent)
% To run EXACTLY N iterations with no early stop: set tol = 0; maxIter = N;
%% ---------------------------------------------------------------------

if Nx < 3 || Ny < 3, error('Nx and Ny must both be at least 3.'); end

% ----- grid -----
x  = linspace(ax,bx,Nx);   hx = x(2) - x(1);
y  = linspace(ay,by,Ny);   hy = y(2) - y(1);
[X,Y] = meshgrid(x,y);     % U(j,i) approximates u( x(i), y(j) )

cx = 1/hx^2;   cy = 1/hy^2;   denom = 2*(cx + cy);

% ----- source term evaluated once on the whole grid -----
F = f(X,Y);
if isscalar(F), F = F * ones(Ny,Nx); end       % allow constant f

% ----- initial array: boundary values imposed, zeros in the interior -----
U = zeros(Ny,Nx);
U(1 ,:)  = evalBC(g_bottom, x);
U(Ny,:)  = evalBC(g_top   , x);
U(: ,1)  = evalBC(g_left  , y).';
U(: ,Nx) = evalBC(g_right , y).';

fprintf('Poisson / Jacobi.  Grid %d x %d,  hx = %.4g,  hy = %.4g\n', ...
        Nx, Ny, hx, hy);

% ----- Jacobi iteration -----
errHist = zeros(maxIter,1);
for k = 1:maxIter
    Uold = U;                                  % freeze iteration k
    for j = 2:Ny-1
        for i = 2:Nx-1
            U(j,i) = ( cx*(Uold(j,i+1) + Uold(j,i-1)) ...
                     + cy*(Uold(j+1,i) + Uold(j-1,i)) ...
                     - F(j,i) ) / denom;
        end
    end
    errHist(k) = max(abs(U(:) - Uold(:)));     % infinity-norm change
    if printEvery > 0 && mod(k,printEvery) == 0
        fprintf('   iteration %6d   max change = %.3e\n', k, errHist(k));
    end
    if errHist(k) < tol, break; end
end
errHist = errHist(1:k);

fprintf('Finished: %d iterations, final max change = %.3e\n', k, errHist(end));
if errHist(end) >= tol && tol > 0
    warning('Tolerance not reached -- increase maxIter or loosen tol.');
end

% ----- results -----
figure;
surf(X,Y,U,'EdgeColor','none'); colorbar
xlabel('x'); ylabel('y'); zlabel('u(x,y)');
title(sprintf('Poisson equation, Jacobi (%d iterations)', k));

figure;
contourf(X,Y,U,20); colorbar; axis equal tight
xlabel('x'); ylabel('y');
title('Poisson equation, Jacobi: solution contours');

figure;
semilogy(1:k, errHist, 'LineWidth', 1.2); grid on
xlabel('iteration k'); ylabel('max |U^{(k)} - U^{(k-1)}|');
title('Jacobi convergence history');

% Solution value at any point (xq,yq) inside the domain:
%     uq = interp2(X,Y,U, xq, yq)
% Check against the exact solution of the default example:
%     Uex = sin(pi*X).*sin(pi*Y);
%     fprintf('Max error vs exact: %.3e\n', max(abs(U(:)-Uex(:))));

%% ------------------- helper (leave at end of file) -------------------
function v = evalBC(g, s)
% Evaluates a boundary-condition handle on the vector s and expands
% scalar outputs (e.g. g = @(x) 5) to the correct size.
    v = g(s);
    if isscalar(v), v = repmat(v, size(s)); end
    v = reshape(v, size(s));
end