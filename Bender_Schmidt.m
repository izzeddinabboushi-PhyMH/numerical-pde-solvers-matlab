%% Bender-Schmidt (Schmidt) method for u_t = c^2 u_xx
% BCs : u(0,t) = u(4,t) = 0
% IC  : u(x,0) = x.*(4-x)
% Given: c^2 = 0.5, h = 1, alpha = 0.5  =>  k = alpha*h^2/c^2 = 1
clc; clear;

c2    = 0.5;               % c^2
h     = 1;                 % space step
alpha = 0.5;               % alpha = k*c^2/h^2
k     = alpha*h^2/c2;      % time step -> 1

x = 0:h:4;                 % nodes i = 0,...,4  (domain is [0,4])
t = 0:k:5*k;               % levels j = 0,...,5
N = numel(x);
M = numel(t);

u = zeros(M,N);
u(1,:) = x.*(4-x);         % initial condition
u(1,1) = 0;  u(1,N) = 0;   % enforce BCs at t=0 (consistent here)

% Bender-Schmidt relation (valid since alpha = 1/2):
%   u(i,j+1) = ( u(i-1,j) + u(i+1,j) ) / 2
for j = 1:M-1
    for i = 2:N-1
        u(j+1,i) = 0.5*( u(j,i-1) + u(j,i+1) );
    end
end

% Print table (rows: t_j, columns: x_i)
fprintf('%8s', 't \ x');  fprintf('%9.1f', x);  fprintf('\n');
for j = 1:M
    fprintf('%8.2f', t(j));  fprintf('%9.4f', u(j,:));  fprintf('\n');
end