clc;
clear;
close all;

% Given data
L = 2;
n = 4;                 % Number of spatial intervals
h = L/n;               % Spatial step: h = 0.5
k = 0.05;              % Time step
T = 0.1;
r = k/h^2;             % r = 0.2

x = (0:h:L)';          % Spatial grid
u = sin(pi*x/2);       % Initial condition

% Enforce boundary conditions
u(1) = 0;
u(end) = 0;

% Number of interior points
m = n - 1;

% Crank-Nicolson matrices
A = (1+r)*eye(m) ...
    - (r/2)*diag(ones(m-1,1),1) ...
    - (r/2)*diag(ones(m-1,1),-1);

B = (1-r)*eye(m) ...
    + (r/2)*diag(ones(m-1,1),1) ...
    + (r/2)*diag(ones(m-1,1),-1);

% Interior initial values
U = u(2:end-1);

% Number of time steps
Nt = round(T/k);

% Crank-Nicolson time iteration
for j = 1:Nt
    U = A \ (B*U);
end

% Complete solution including boundaries
u_final = [0; U; 0];

% Display the result
fprintf('Approximate solution at t = %.2f:\n\n',T);
fprintf('      x           u(x,t)\n');
fprintf('-------------------------\n');

for i = 1:length(x)
    fprintf('%8.2f      %10.6f\n',x(i),u_final(i));
end

% Plot
figure;
plot(x,u_final,'o-','LineWidth',1.5,'MarkerSize',7);
xlabel('x');
ylabel('u(x,0.1)');
title('Crank-Nicolson Approximation at t = 0.1');
grid on;