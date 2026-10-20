%==========================================================================
% LR.m   : LAPLACE EQUATION   u_xx + u_yy = 0    on [ax,bx] x [ay,by]
% METHOD : RELAXATION METHOD (SOR)   (5-point finite differences)
%          Each node first gets its Gauss-Seidel value uGS, then is
%          over-relaxed:  U(i,j) <- (1-omega)*U(i,j) + omega*uGS,
%          with 0 < omega < 2  (omega = 1 reproduces Gauss-Seidel).
% BCs    : Dirichlet on all four edges, typed directly as formulas.
%
% HOW TO USE : edit ONLY the "USER INPUT" block. Boundary functions are
% evaluated point-by-point, so ordinary scalar syntax (*, ^) is fine.
% Leave omega = 0 and the code computes the OPTIMAL relaxation factor
% for your grid by itself:  omega_opt = 2/(1 + sqrt(1 - rhoJ^2)),
% rhoJ = ( cx*cos(pi/(Nx-1)) + cy*cos(pi/(Ny-1)) )/(cx+cy)  being the
% Jacobi spectral radius of the discrete problem (Young's formula).
%
% NOTE   : omega-related statements use elementwise operators (.*, ./, .^)
% so the Code Analyzer can verify they are scalar (for scalars they equal
% *, /, ^ exactly). Plot commands are written ONE PER LINE so that no
% statement depends on a comma separator.
%==========================================================================
clc; clear; close all;

%============================ USER INPUT ==================================
ax = 0;  bx = 1;                  % x-interval :  ax <= x <= bx
ay = 0;  by = 1;                  % y-interval :  ay <= y <= by

Nx = 41;                          % grid points in x (boundaries included)
Ny = 41;                          % grid points in y (boundaries included)

gBottom = @(x) 0;                 % u(x,ay)   bottom edge
gTop    = @(x) sin(pi*x);         % u(x,by)   top edge
gLeft   = @(y) 0;                 % u(ax,y)   left edge
gRight  = @(y) 0;                 % u(bx,y)   right edge
% (Demo problem: exact solution u = sin(pi*x)*sinh(pi*y)/sinh(pi).)

maxIter = 20000;                  % iterate up to ANY number of sweeps
tol     = 1e-8;                   % stop when max|U_new - U_old| < tol
                                  % set tol = 0 to force exactly maxIter sweeps
omega   = 0;                      % relaxation factor, 0 < omega < 2
                                  % leave omega = 0 => optimal value is
                                  % computed automatically below
%==========================================================================

% ---- grid and spacings (computed automatically) -------------------------
x  = linspace(ax,bx,Nx);   hx = x(2)-x(1);
y  = linspace(ay,by,Ny);   hy = y(2)-y(1);

cx = 1/hx^2;   cy = 1/hy^2;   denom = 2*(cx+cy);

% ---- optimal relaxation factor if the user did not fix one --------------
if omega <= 0
    rhoJ  = ( cx.*cos(pi./(Nx-1)) + cy.*cos(pi./(Ny-1)) ) ./ (cx + cy);
    omega = 2./(1 + sqrt(1 - rhoJ.^2));
end
fprintf('Relaxation factor in use: omega = %.4f\n', omega);

% ---- initial array with the boundary conditions imposed -----------------
U = zeros(Ny,Nx);                        % U(i,j) approximates u(x(j),y(i))
U(1 ,:) = arrayfun(gBottom,x);           % edge y = ay
U(Ny,:) = arrayfun(gTop   ,x);           % edge y = by
U(: ,1) = arrayfun(gLeft  ,y).';         % edge x = ax
U(:,Nx) = arrayfun(gRight ,y).';         % edge x = bx
% (Corners take the last assignment; with consistent data this is immaterial.)

% ---- relaxation (SOR) iteration ------------------------------------------
errHist = zeros(maxIter,1);
for k = 1:maxIter
    Uold = U;                            % stored ONLY to measure the change
    for i = 2:Ny-1
        for j = 2:Nx-1
            uGS    = ( cx.*(U(i,j-1)+U(i,j+1)) ...
                     + cy.*(U(i-1,j)+U(i+1,j)) ) ./ denom; % Gauss-Seidel value
            U(i,j) = (1 - omega).*U(i,j) + omega.*uGS;     % over-relaxation
        end
    end
    errHist(k) = max(abs(U(:)-Uold(:)));
    if errHist(k) < tol
        break
    end
end
nIter   = k;
errHist = errHist(1:nIter);

% ---- automatic verification: discrete u_xx+u_yy at interior nodes -------
R = cx*( U(2:end-1,1:end-2) - 2*U(2:end-1,2:end-1) + U(2:end-1,3:end) ) ...
  + cy*( U(1:end-2,2:end-1) - 2*U(2:end-1,2:end-1) + U(3:end,2:end-1) );

fprintf('RELAXATION (SOR), Laplace equation\n');
if errHist(nIter) < tol
    fprintf('Converged: tolerance %.1e reached after %d sweeps.\n', tol, nIter);
else
    fprintf('Stopped at the iteration limit maxIter = %d sweeps.\n', maxIter);
end
fprintf('Last change per sweep   max|dU|      = %.3e\n', errHist(nIter));
fprintf('PDE check               max|Lap(u)|  = %.3e (should be near 0)\n', ...
        max(abs(R(:))));

% ---- plots (one command per line: nothing depends on a comma) ------------
figure
surf(x,y,U)
shading interp
colorbar
xlabel('x')
ylabel('y')
zlabel('u(x,y)')
title(sprintf('Laplace equation - Relaxation, \\omega = %.3f (%d iterations)', ...
      omega, nIter))

figure
contourf(x,y,U,20)
colorbar
axis equal tight
xlabel('x')
ylabel('y')
title('Contours of u')

figure
semilogy(1:nIter, errHist, 'LineWidth', 1.2)
grid on
xlabel('iteration k')
ylabel('max|U^{(k)} - U^{(k-1)}|')
title('Relaxation (SOR) convergence history')