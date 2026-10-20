clc; clear;

A = [ 4 -1 -1  0;
     -1  4  0 -1;
     -1  0  4 -1;
      0 -1 -1  4];
B = [100; 100; 200; 200];

U = zeros(4,1);        % initial guess
D = diag(diag(A));     % = 4*eye(4)
R = A - D;

for k = 1:3            % Jacobi iteration
    U = D \ (B - R*U);
    fprintf('--- Iteration %d ---\n', k)
    disp(U')
end