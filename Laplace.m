%% Poisson equation:  u_xx + u_yy = -2(x^2 + y^2)
%  Square 0 <= x,y <= 3,  u = 0 on the boundary,  mesh length h = 1
%  Solved by Jacobi and Gauss-Seidel iteration
clc; clear;

h    = 1;
x    = 0:h:3;                  % 4 grid lines in x
y    = 0:h:3;                  % 4 grid lines in y
N    = numel(x);
tol  = 1e-8;
maxit= 500;

% source term f(x,y) = -2(x^2 + y^2)
f = @(X,Y) -2*(X.^2 + Y.^2);

[X,Y] = meshgrid(x,y);         % row index -> y, col index -> x
F     = f(X,Y);

%% ---------------- Jacobi ----------------
U = zeros(N,N);                % includes boundaries (all zero)
fprintf('=== JACOBI ===\n');
for it = 1:maxit
    Uold = U;
    for i = 2:N-1              % y index
        for j = 2:N-1          % x index
            U(i,j) = 0.25*( Uold(i-1,j) + Uold(i+1,j) ...
                + Uold(i,j-1) + Uold(i,j+1) ...
                - h^2*F(i,j) );
        end
    end
    err = max(abs(U(:)-Uold(:)));
    if it <= 5 || err < tol
        fprintf('iter %3d :  u1=%8.5f  u2=%8.5f  u3=%8.5f  u4=%8.5f   err=%.2e\n',...
            it, U(3,2), U(3,3), U(2,2), U(2,3), err);
    end
    if err < tol, break; end
end
U_jac = U;
fprintf('Jacobi converged in %d iterations\n\n', it);

%% ------------- Gauss-Seidel -------------
U = zeros(N,N);
fprintf('=== GAUSS-SEIDEL ===\n');
for it = 1:maxit
    Uold = U;
    for i = N-1:-1:2           % sweep top row first (u1,u2 then u3,u4)
        for j = 2:N-1
            U(i,j) = 0.25*( U(i-1,j) + U(i+1,j) ...
                + U(i,j-1) + U(i,j+1) ...
                - h^2*F(i,j) );
        end
    end
    err = max(abs(U(:)-Uold(:)));
    if it <= 5 || err < tol
        fprintf('iter %3d :  u1=%8.5f  u2=%8.5f  u3=%8.5f  u4=%8.5f   err=%.2e\n',...
            it, U(3,2), U(3,3), U(2,2), U(2,3), err);
    end
    if err < tol, break; end
end
U_gs = U;
fprintf('Gauss-Seidel converged in %d iterations\n\n', it);

%% ---------------- Results ----------------
disp('Jacobi solution (rows: y = 3,2,1,0):');
disp(flipud(U_jac));
disp('Gauss-Seidel solution (rows: y = 3,2,1,0):');
disp(flipud(U_gs));

fprintf('Exact interior values:  u1 = u4 = %.6f,  u2 = %.6f,  u3 = %.6f\n',...
    17/3, 23/3, 14/3);