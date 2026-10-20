function Hyperbolic_NPDE_Solver
%HYPERBOLIC_NPDE_SOLVER  1-D wave equation: four numerical methods.
%
%   Solves   u_tt = c^2*u_xx  on 0 < x < L, t > 0,  with Dirichlet
%   conditions u(0,t) = u(L,t) = 0 and initial data u(x,0) = f(x),
%   u_t(x,0) = g(x), and compares four finite-difference methods with
%   the exact Fourier-series solution.  With r = c*dt/dx and
%   dxx(v) = v(i+1)-2v(i)+v(i-1):
%
%   M1 solveExplicitCTCS      U(i,n+1) = 2U(i,n)-U(i,n-1)+r^2*dxx(U)
%                             3 levels, O(dx^2+dt^2), stable iff r<=1
%   M2 solveImplicitTheta     (I-th*r^2*dxx)U(:,n+1) = ... tridiagonal
%                             solve, unconditionally stable, th >= 1/4
%   M3 solveLaxWendroff       Lax-Wendroff on p = u_t, q = c*u_x, then
%                             u by the trapezoidal rule; dissipative
%   M4 solveCharacteristics   U(i,n+1) = U(i+1,n)+U(i-1,n)-U(i,n-1)
%                             at r = 1; exact when g = 0
%
%   Examples: 1 fundamental string mode (g = 0); 2 superposed modes
%   with g not 0; 3 PHYSICS/EM, a standing wave in a vacuum cavity
%   between two perfect electric conductors, with H recovered from
%   Faraday's law and the stored energy monitored.  A CFL experiment
%   and a grid-refinement study follow.
%
%   Units: x [m], t [s], c [m/s], u [m] (string) or [V/m] (E field).
%   Usage: hyperbolic_npde_solver.  No toolboxes; MATLAB or Octave.

close all;
fprintf('\n%s\n', repmat('=', 1, 68));
fprintf('  HYPERBOLIC PDEs -- 1-D WAVE EQUATION\n');
fprintf('  Four numerical methods verified against exact solutions\n');
fprintf('%s\n', repmat('=', 1, 68));

runExample1();
runExample2();
runExample3EM();
runCFLExperiment();
runConvergenceStudy();
fprintf('\nAll three examples completed successfully.\n\n');
end

%% ---- EXAMPLE 1: fundamental standing mode of a string -----------
function runExample1
%RUNEXAMPLE1  f = sin(pi*x/L), g = 0, c = 2 m/s, L = 1 m.
%   Exact: u(x,t) = sin(pi*x/L)*cos(c*pi*t/L), period T1 = 2L/c = 1 s.
c = 2.0;  L = 1.0;  T = 1.0;  Nx = 100;  rTarget = 0.90;
f   = @(x) sin(pi*x/L);
g   = @(x) zeros(size(x));
uex = @(x,t) sin(pi*x/L) .* cos(c*pi*t/L);

printHeader('EXAMPLE 1', 'Fundamental standing mode of a string');
fprintf('  c = %.4f m/s, L = %.4f m, T = %.4f s\n', c, L, T);

fprintf('\n  (a) Coarse grid reproducing the hand computation\n');
[Uh, xh, th] = solveExplicitCTCS(c, L, 3/32, 4, 0.50, f, g);
printLevelTable(Uh, buildExactField(uex, xh, th), xh, th, 0:2);

fprintf('\n  (b) Four methods on Nx = %d\n', Nx);
compareAllMethods(c, L, T, Nx, rTarget, f, g, uex, 'Example 1');

[U, x, t] = solveExplicitCTCS(c, L, T, Nx, rTarget, f, g);
Uex = buildExactField(uex, x, t);
plotProfiles(U, Uex, x, t, [0.125 0.375 0.625 0.875], ...
    'Example 1: fundamental mode', 'u  [m]');
plotErrorHistory(U, Uex, x, t, 'Example 1');
end

%% ---- EXAMPLE 2: superposed modes, non-zero initial velocity -----
function runExample2
%RUNEXAMPLE2  f = sin(pi*x/L)+0.5*sin(3*pi*x/L),
%   g = (c*pi/2L)*sin(2*pi*x/L).  With omega_n = c*n*pi/L the only
%   non-zero coefficients are A1 = 1, A3 = 1/2 and B2 = 1/4, so the
%   Fourier series terminates and the exact solution is a finite sum.
c = 1.0;  L = 1.0;  T = 2.0;  Nx = 120;  rTarget = 0.80;
w   = @(n) c*n*pi/L;
f   = @(x) sin(pi*x/L) + 0.5*sin(3*pi*x/L);
g   = @(x) (c*pi/(2*L)) * sin(2*pi*x/L);
uex = @(x,t) 1.00*sin(1*pi*x/L).*cos(w(1)*t) ...
           + 0.25*sin(2*pi*x/L).*sin(w(2)*t) ...
           + 0.50*sin(3*pi*x/L).*cos(w(3)*t);

printHeader('EXAMPLE 2', 'Superposed modes, non-zero initial velocity');
fprintf('  c = %.4f m/s, L = %.4f m, T = %.4f s\n', c, L, T);

fprintf('\n  (a) Coarse grid reproducing the hand computation\n');
[Uh, xh, th] = solveExplicitCTCS(c, L, 1/6, 6, 0.50, f, g);
printLevelTable(Uh, buildExactField(uex, xh, th), xh, th, 0:2);

fprintf('\n  (b) Four methods on Nx = %d\n', Nx);
compareAllMethods(c, L, T, Nx, rTarget, f, g, uex, 'Example 2');

[U, x, t] = solveExplicitCTCS(c, L, T, Nx, rTarget, f, g);
Uex = buildExactField(uex, x, t);
plotProfiles(U, Uex, x, t, [0.00 0.25 0.50 0.75], ...
    'Example 2: superposed modes', 'u  [m]');
plotErrorHistory(U, Uex, x, t, 'Example 2');
end

%% ---- EXAMPLE 3: physics/EM, a 1-D PEC cavity resonator ----------
function runExample3EM
%RUNEXAMPLE3EM  Standing wave between two perfect electric conductors
%   a distance d apart.  For E = E_y(x,t)*yhat, H = H_z(x,t)*zhat
%   Maxwell's curl equations reduce to
%       dE_y/dx = -mu0*dH_z/dt ,   dH_z/dx = -eps0*dE_y/dt ,
%   hence  d2E_y/dt2 = c^2*d2E_y/dx2  with c = 1/sqrt(mu0*eps0).  A
%   PEC forces the tangential E to vanish, so E_y(0,t) = E_y(d,t) = 0
%   are physics, not modelling convenience, and f_n = n*c/(2*d).  With
%   E_y(x,0) = E0*sin(pi*x/d) and dE_y/dt(x,0) = 0 the exact fields
%   are E_y = E0*sin(pi*x/d)*cos(w1*t) and
%   H_z = -(E0/eta0)*cos(pi*x/d)*sin(w1*t): a quarter period apart,
%   with the stored energy eps0*E0^2*d/4 constant.
eps0 = 8.8541878128e-12;                 % [F/m]
mu0  = 1.25663706212e-6;                 % [H/m]
c    = 1/sqrt(mu0*eps0);                 % [m/s]
eta0 = sqrt(mu0/eps0);                   % [ohm]
d  = 0.05;                               % PEC separation [m]
E0 = 1.0;                                % field amplitude [V/m]
f1 = c/(2*d);                            % fundamental resonance [Hz]
T1 = 1/f1;                               % period [s]
w1 = pi*c/d;                             % angular frequency [rad/s]
Nx = 200;  rTarget = 0.90;  T = 2*T1;

f   = @(x) E0*sin(pi*x/d);
g   = @(x) zeros(size(x));
uex = @(x,t)  E0*sin(pi*x/d).*cos(w1*t);
Hex = @(x,t) -(E0/eta0)*cos(pi*x/d).*sin(w1*t);

printHeader('EXAMPLE 3', 'PHYSICS/EM: standing wave in a 1-D PEC cavity');
fprintf('  c    = %.9e m/s   (1/sqrt(mu0*eps0))\n', c);
fprintf('  eta0 = %.6f ohm      (sqrt(mu0/eps0))\n', eta0);
fprintf('  d    = %.4f m, E0 = %.2f V/m\n', d, E0);
fprintf('  f1   = %.6f GHz, T1 = %.6f ps\n', f1/1e9, T1*1e12);

fprintf('\n  (a) Coarse grid reproducing the hand computation\n');
[Eh, xh, tCh] = solveExplicitCTCS(c, d, 5.5e-11, 5, 0.80, f, g);
printLevelTable(Eh, buildExactField(uex, xh, tCh), xh, tCh, 0:2);

fprintf('\n  (b) Four methods on Nx = %d\n', Nx);
compareAllMethods(c, d, T, Nx, rTarget, f, g, uex, 'Example 3');

% ---- physics checks on the explicit run -------------------------
[E, x, t, r] = solveExplicitCTCS(c, d, T, Nx, rTarget, f, g);
Eexact = buildExactField(uex, x, t);
[H, xH, tH] = recoverMagneticField(E, x, t, mu0);
Hexact = buildExactField(Hex, xH, tH);
fMeas  = measureFrequency(t, E(round((Nx+1)/2), :));
W      = electromagneticEnergy(E, H, x, eps0, mu0);
Wex    = 0.25*eps0*E0^2*d;

fprintf('\n  Physics verification\n');
fprintf('  analytical resonance f1        : %.9f GHz\n', f1/1e9);
fprintf('  measured  resonance (zero xing): %.9f GHz\n', fMeas/1e9);
fprintf('  relative frequency error       : %.3e\n', abs(fMeas-f1)/f1);
fprintf('  max|H_z| numerical / analytical: %.6e / %.6e A/m\n', ...
    max(abs(H(:))), E0/eta0);
fprintf('  max error in H_z               : %.3e A/m\n', ...
    max(max(abs(H-Hexact))));
fprintf('  E/H peak ratio (should be eta0): %.4f ohm\n', ...
    max(abs(E(:)))/max(abs(H(:))));
fprintf('  stored energy: exact %.6e J/m^2, drift %.3e\n', ...
    Wex, (max(W)-min(W))/Wex);
fprintf('  CFL number used r = %.4f, dt = %.6f ps\n', r, ...
    (t(2)-t(1))*1e12);

plotProfiles(E, Eexact, x, t, T1*[0 0.125 0.25 0.375], ...
    'Example 3: PEC cavity, E_y', 'E_y  [V/m]');
im = round((Nx+1)/2);
plotFieldHistory(t*1e12, [Eexact(im,:); E(im,:)], ...
    tH*1e12, 1e3*[Hexact(1,:); H(1,:)]);
end

%% ---- STUDY A: the CFL stability experiment ----------------------
function runCFLExperiment
%RUNCFLEXPERIMENT  One physical problem, three CFL numbers.  The PDE,
%   the data and the exact solution are those of Example 1; only dt
%   changes, so every difference belongs to the discretisation.
c = 2.0;  L = 1.0;  T = 1.0;  Nx = 50;
rList = [0.50 1.00 1.05];
f   = @(x) sin(pi*x/L);
g   = @(x) zeros(size(x));
uex = @(x,t) sin(pi*x/L).*cos(c*pi*t/L);

printHeader('STUDY A', 'CFL experiment: r < 1, r = 1, r > 1');
fprintf('\n  %-8s %-12s %-12s %-14s %-14s %-9s\n', 'r', 'dt [s]', ...
    't_end [s]', 'max|U| [m]', 'Linf err [m]', 'verdict');
fprintf('  %s\n', repmat('-', 1, 74));
for k = 1:numel(rList)
    [U, x, t, r] = solveExplicitCTCS(c, L, T, Nx, rList(k), f, g);
    [~, ~, eInf] = calculateErrors(U, buildExactField(uex, x, t), x);
    if r > 1 + 1e-12
        verdict = 'UNSTABLE';
    elseif abs(r-1) < 1e-12
        verdict = 'exact';
    else
        verdict = 'stable';
    end
    fprintf('  %-8.3f %-12.6f %-12.6f %-14.4e %-14.4e %-9s\n', ...
        r, t(2)-t(1), t(end), max(abs(U(:,end))), eInf(end), verdict);
end
r = 1.05;
fprintf('\n  r = 1.05: predicted amplification of the sawtooth\n');
fprintf('  mode (theta = pi) per step, |xi| = 2r^2-1+2r*sqrt(r^2-1)');
fprintf(' = %.6f\n', 2*r^2 - 1 + 2*r*sqrt(r^2-1));
end

%% ---- STUDY B: observed order of accuracy ------------------------
function runConvergenceStudy
%RUNCONVERGENCESTUDY  Order p = log2( E(dx)/E(dx/2) ) for each method,
%   with r held fixed so that dt and dx refine together.  A non-zero g
%   is added so that Method 4 is not exact and a genuine order can be
%   measured for all four.  A second-order method must give p -> 2.
c = 2.0;  L = 1.0;  T = 0.5;  rTarget = 0.80;
gridList = [25 50 100 200 400];
f   = @(x) sin(pi*x/L);
g   = @(x) (c*pi/(2*L))*sin(2*pi*x/L);
uex = @(x,t) sin(pi*x/L).*cos(c*pi*t/L) ...
           + 0.25*sin(2*pi*x/L).*sin(2*c*pi*t/L);

printHeader('STUDY B', 'Grid-refinement study at fixed CFL number');
names = {'M1 explicit CTCS', 'M2 implicit theta=1/4', ...
         'M3 Lax-Wendroff', 'M4 characteristics'};
err = zeros(4, numel(gridList));
for k = 1:numel(gridList)
    Nx = gridList(k);
    err(1,k) = peakError(@() solveExplicitCTCS(c,L,T,Nx,rTarget,f,g), uex);
    err(2,k) = peakError(@() solveImplicitTheta(c,L,T,Nx,rTarget,f,g,0.25), uex);
    err(3,k) = peakError(@() solveLaxWendroff(c,L,T,Nx,rTarget,f,g), uex);
    err(4,k) = peakError(@() solveCharacteristics(c,L,T,Nx,f,g), uex);
end
for m = 1:4
    fprintf('\n  %s\n', names{m});
    fprintf('  %-8s %-12s %-16s %-10s\n', 'Nx', 'dx [m]', ...
        'Linf error', 'order p');
    fprintf('  %s\n', repmat('-', 1, 50));
    for k = 1:numel(gridList)
        dx = L/gridList(k);
        if k == 1
            fprintf('  %-8d %-12.6f %-16.6e %-10s\n', ...
                gridList(k), dx, err(m,k), '--');
        else
            fprintf('  %-8d %-12.6f %-16.6e %-10.3f\n', gridList(k), ...
                dx, err(m,k), log2(err(m,k-1)/err(m,k)));
        end
    end
end
figure('Name', 'Grid-refinement study', 'NumberTitle', 'off', ...
       'Color', 'w');
dxList = L ./ gridList;
marks = {'o-', 's-', 'd-', '^-'};
for m = 1:4
    loglog(dxList, max(err(m,:), eps), marks{m}, 'LineWidth', 1.6, ...
        'MarkerSize', 6); hold on;
end
loglog(dxList, err(1,1)*(dxList/dxList(1)).^2, 'k--', 'LineWidth', 1.2);
hold off; grid on;
xlabel('\Deltax  [m]'); ylabel('L_\infty error');
title('Observed convergence at fixed CFL number');
legend([names, {'slope 2'}], 'Location', 'northwest');
end

%% ---- METHOD 1: explicit centred scheme (CTCS, leap-frog) --------
function [U, x, t, r] = solveExplicitCTCS(c, L, T, Nx, rTarget, f, g)
%SOLVEEXPLICITCTCS  U(i,n+1) = 2U-U(old) + r^2*dxx(U), r = c*dt/dx.
%   Three-level, matrix-free, O(dx^2+dt^2), stable iff r <= 1.
[x, t, dx, dt, r, Nt] = makeGrid(c, L, T, Nx, rTarget);       %#ok<ASGLU>
if r > 1 + 1e-12
    warning('npde:cfl', 'r = %.4f > 1: CTCS is unstable.', r);
end
U      = zeros(Nx+1, Nt+1);
U(:,1) = f(x);                              % level n = 0
U([1 end], 1) = 0;                          % Dirichlet corners
in = 2:Nx;                                  % interior node indices
U(in,2) = U(in,1) + dt*g(x(in)) ...
        + 0.5*r^2*( U(in+1,1) - 2*U(in,1) + U(in-1,1) );
last = Nt + 1;
for n = 2:Nt
    U(in,n+1) = 2*U(in,n) - U(in,n-1) ...
              + r^2*( U(in+1,n) - 2*U(in,n) + U(in-1,n) );
    U([1 end], n+1) = 0;                    % re-impose the BCs
    if ~all(isfinite(U(:,n+1))) || max(abs(U(:,n+1))) > 1e12
        last = n + 1;
        fprintf('  *** blow-up guard at t = %.4e s ***\n', t(n+1));
        break;
    end
end
U = U(:, 1:last);  t = t(1:last);
end

%% ---- METHOD 2: implicit theta-weighted scheme -------------------
function [U, x, t, r] = solveImplicitTheta(c, L, T, Nx, rTarget, ...
                                          f, g, theta)
%SOLVEIMPLICITTHETA  Weighted three-level scheme,
%   ( I - theta*r^2*dxx ) U(:,n+1) = 2U(:,n) - U(:,n-1)
%       + (1-2*theta)*r^2*dxx( U(:,n) ) + theta*r^2*dxx( U(:,n-1) ).
%   Unconditionally stable for theta >= 1/4; theta = 0 recovers
%   Method 1; theta = 1/4 is Newmark average acceleration.
[x, t, dx, dt, r, Nt] = makeGrid(c, L, T, Nx, rTarget);       %#ok<ASGLU>
m  = Nx - 1;                                % number of unknowns
lo = -theta*r^2 * ones(m,1);                % sub-diagonal
di = (1 + 2*theta*r^2) * ones(m,1);         % main diagonal
up = -theta*r^2 * ones(m,1);                % super-diagonal
U      = zeros(Nx+1, Nt+1);
U(:,1) = f(x);
U([1 end], 1) = 0;
gv = g(x);
% start-up level, from the same ghost-node elimination
rhs = U(2:Nx,1) + dt*gv(2:Nx) ...
    + 0.5*(1-2*theta)*r^2*dxx(U(:,1), Nx) ...
    - theta*r^2*dt*dxx(gv, Nx);
U(2:Nx,2) = thomas(lo, di, up, rhs);
U([1 end], 2) = 0;
for n = 2:Nt
    rhs = 2*U(2:Nx,n) - U(2:Nx,n-1) ...
        + (1-2*theta)*r^2*dxx(U(:,n), Nx) ...
        + theta*r^2*dxx(U(:,n-1), Nx);
    U(2:Nx,n+1) = thomas(lo, di, up, rhs);
    U([1 end], n+1) = 0;
end
end

%% ---- METHOD 3: Lax-Wendroff on the first-order system -----------
function [U, x, t, r] = solveLaxWendroff(c, L, T, Nx, rTarget, f, g)
%SOLVELAXWENDROFF  p = u_t, q = c*u_x, then Lax-Wendroff on
%   p_t = c*q_x, q_t = c*p_x, and u recovered by the trapezoidal rule.
%   Two levels, O(dx^2+dt^2), dissipative, stable iff r <= 1.  At a
%   Dirichlet wall u vanishes for all t, so p is odd and q is even
%   about the wall; the ghost values below encode exactly that
%   symmetry and keep the boundary second-order accurate.
[x, t, dx, dt, r, Nt] = makeGrid(c, L, T, Nx, rTarget);
if r > 1 + 1e-12
    warning('npde:cfl', 'r = %.4f > 1: Lax-Wendroff is unstable.', r);
end
U      = zeros(Nx+1, Nt+1);
U(:,1) = f(x);
U([1 end], 1) = 0;
p = g(x);  p([1 end]) = 0;                   % p = u_t
fv = U(:,1);
q = zeros(Nx+1, 1);                          % q = c*u_x
q(2:Nx) = c*( fv(3:Nx+1) - fv(1:Nx-1) ) / (2*dx);
q(1)    = c*( fv(2) + fv(2) ) / (2*dx);      % odd mirror f(-1)=-f(1)
q(Nx+1) = c*( -fv(Nx) - fv(Nx) ) / (2*dx);
for n = 1:Nt
    pE = [-p(2); p; -p(Nx)];                 % p odd  about the walls
    qE = [ q(2); q; q(Nx)];                  % q even about the walls
    j  = 2:Nx+2;                             % extended interior
    pNew = pE(j) + 0.5*r*( qE(j+1) - qE(j-1) ) ...
         + 0.5*r^2*( pE(j+1) - 2*pE(j) + pE(j-1) );
    qNew = qE(j) + 0.5*r*( pE(j+1) - pE(j-1) ) ...
         + 0.5*r^2*( qE(j+1) - 2*qE(j) + qE(j-1) );
    U(:,n+1) = U(:,n) + 0.5*dt*( p + pNew );
    p = pNew;  q = qNew;
    p([1 end]) = 0;                          % exact at a PEC/fixed end
    U([1 end], n+1) = 0;
end
end

%% ---- METHOD 4: method of characteristics (d'Alembert) -----------
function [U, x, t, r] = solveCharacteristics(c, L, T, Nx, f, g)
%SOLVECHARACTERISTICS  Unit-CFL marching on the characteristic grid.
%   With dt = dx/c exactly, the CTCS stencil collapses to
%       U(i,n+1) = U(i+1,n) + U(i-1,n) - U(i,n-1),
%   d'Alembert's formula on the grid: the centre node disappears and
%   no dispersion error is committed.  The start-up level is the
%   characteristic average U(i,1) = ( f(i+1)+f(i-1) )/2 + dt*g(i).
%   dt is not a free parameter: r = 1 always.
[x, t, dx, dt, r, Nt] = makeGrid(c, L, T, Nx, 1.0);           %#ok<ASGLU>
U      = zeros(Nx+1, Nt+1);
U(:,1) = f(x);
U([1 end], 1) = 0;
in = 2:Nx;
U(in,2) = 0.5*( U(in+1,1) + U(in-1,1) ) + dt*g(x(in));
U([1 end], 2) = 0;
for n = 2:Nt
    U(in,n+1) = U(in+1,n) + U(in-1,n) - U(in,n-1);
    U([1 end], n+1) = 0;
end
end

%% ---- electromagnetic post-processing (Example 3) ----------------
function [H, xH, tH] = recoverMagneticField(E, x, t, mu0)
%RECOVERMAGNETICFIELD  Faraday's law on the staggered Yee grid,
%       H(i+1/2,n+1/2) = H(i+1/2,n-1/2)
%                      - dt/(mu0*dx) * ( E(i+1,n) - E(i,n) ),
%   which is exactly Yee's FDTD update.  H lives at the cell centres
%   xH and the half time levels tH, in time quadrature with E.  H_z is
%   odd in time about t = 0, so H(-1/2) = -H(+1/2) and the first half
%   step carries the factor 1/2.
dx = x(2) - x(1);
dt = t(2) - t(1);
xH = 0.5*( x(1:end-1) + x(2:end) );
tH = t(1:end-1) + 0.5*dt;
H  = zeros(numel(xH), numel(tH));
H(:,1) = -0.5*dt/(mu0*dx) * ( E(2:end,1) - E(1:end-1,1) );
for n = 2:numel(tH)
    H(:,n) = H(:,n-1) - dt/(mu0*dx) * ( E(2:end,n) - E(1:end-1,n) );
end
end

function W = electromagneticEnergy(E, H, x, eps0, mu0)
%ELECTROMAGNETICENERGY  Stored energy per unit area [J/m^2],
%   W = integral( eps0*E^2/2 + mu0*H^2/2 ) dx, with the magnetic term
%   averaged over the two half time levels straddling t^n.
dx = x(2) - x(1);
nT = size(H,2);
W  = zeros(1, nT-1);
for n = 2:nT
    we = 0.5*eps0*dx*( sum(E(2:end-1,n).^2) ...
                     + 0.5*(E(1,n)^2 + E(end,n)^2) );
    wh = 0.5*mu0 *dx* sum( 0.5*(H(:,n-1).^2 + H(:,n).^2) );
    W(n-1) = we + wh;
end
end

function fMeas = measureFrequency(t, y)
%MEASUREFREQUENCY  Frequency from linearly interpolated zero
%   crossings: consecutive crossings are half a period apart.
s = sign(y);
k = find(s(1:end-1).*s(2:end) < 0);
tz = t(k) - y(k).*(t(k+1)-t(k))./(y(k+1)-y(k));
fMeas = 1 / (2*mean(diff(tz)));
end

%% ---- shared numerical utilities ---------------------------------
function [x, t, dx, dt, r, Nt] = makeGrid(c, L, T, Nx, rTarget)
%MAKEGRID  Uniform space-time grid; dt is derived from the CFL number,
%   so stability is an input choice and never an accident.
dx = L / Nx;
dt = rTarget * dx / c;
Nt = max(1, round(T/dt));
x  = transpose(linspace(0, L, Nx+1));
t  = (0:Nt) * dt;
r  = c * dt / dx;
end

function d2 = dxx(v, Nx)
%DXX  Centred second difference on the interior nodes, delta_x^2 v.
d2 = v(3:Nx+1) - 2*v(2:Nx) + v(1:Nx-1);
end

function y = thomas(lo, di, up, rhs)
%THOMAS  Tridiagonal solve (Thomas algorithm), O(m) work.  Solves
%   lo(i)*y(i-1) + di(i)*y(i) + up(i)*y(i+1) = rhs(i).
m = numel(rhs);
cp = zeros(m,1);  dp = zeros(m,1);  y = zeros(m,1);
cp(1) = up(1)/di(1);
dp(1) = rhs(1)/di(1);
for i = 2:m
    den   = di(i) - lo(i)*cp(i-1);
    cp(i) = up(i)/den;
    dp(i) = ( rhs(i) - lo(i)*dp(i-1) )/den;
end
y(m) = dp(m);
for i = m-1:-1:1
    y(i) = dp(i) - cp(i)*y(i+1);
end
end

function Uex = buildExactField(uexact, x, t)
%BUILDEXACTFIELD  Sample the analytical solution on the whole grid, so
%   that it can be compared entrywise with the numerical field.
[tG, xG] = meshgrid(t, x);
Uex = uexact(xG, tG);
end

function [errField, errL2, errInf] = calculateErrors(U, Uex, x)
%CALCULATEERRORS  Pointwise, discrete L2 and L-infinity errors.  The
%   factor dx makes errL2 approximate the continuous L2(0,L) norm, so
%   errors from different grids are directly comparable.
dx       = x(2) - x(1);
errField = abs(U - Uex);
errL2    = sqrt( dx * sum(errField.^2, 1) );
errInf   = max(errField, [], 1);
end

function e = peakError(solverThunk, uexact)
%PEAKERROR  Worst L-infinity error over the whole run of one method.
[U, x, t] = solverThunk();
[~, ~, eInf] = calculateErrors(U, buildExactField(uexact, x, t), x);
e = max(eInf);
end

function compareAllMethods(c, L, T, Nx, rTarget, f, g, uexact, tag)
%COMPAREALLMETHODS  Run the four methods on one problem and tabulate
%   the CFL number actually used, the error norms and the cost.
names = {'M1 explicit CTCS  ', 'M2 implicit th=1/4', ...
         'M3 Lax-Wendroff   ', 'M4 characteristics'};
fprintf('  %-18s %-7s %-9s %-13s %-13s\n', 'method', 'r', 'steps', ...
    'max L2 err', 'max Linf err');
fprintf('  %s\n', repmat('-', 1, 64));
for m = 1:4
    switch m
      case 1, [U,x,t,r] = solveExplicitCTCS(c,L,T,Nx,rTarget,f,g);
      case 2, [U,x,t,r] = solveImplicitTheta(c,L,T,Nx,rTarget,f,g,0.25);
      case 3, [U,x,t,r] = solveLaxWendroff(c,L,T,Nx,rTarget,f,g);
      case 4, [U,x,t,r] = solveCharacteristics(c,L,T,Nx,f,g);
    end
    [~, eL2, eInf] = calculateErrors(U, buildExactField(uexact,x,t), x);
    fprintf('  %-18s %-7.4f %-9d %-13.4e %-13.4e\n', ...
        names{m}, r, numel(t)-1, max(eL2), max(eInf));
end
fprintf('  (%s, Nx = %d, requested r = %.2f)\n', tag, Nx, rTarget);
end

%% ---- reporting and plotting -------------------------------------
function printHeader(label, description)
%PRINTHEADER  Console banner for one example or study.
fprintf('\n%s\n', repmat('-', 1, 68));
fprintf(' %s: %s\n', label, description);
fprintf('%s\n', repmat('-', 1, 68));
end

function printLevelTable(U, Uex, x, t, levels)
%PRINTLEVELTABLE  Full space profile at the first few time levels,
%   which is what a hand computation produces.
for n = levels
    fprintf('  n = %d,  t = %.6e\n', n, t(n+1));
    fprintf('  %-6s %-16s %-16s %-12s\n', 'i', 'numerical', ...
        'analytical', '|error|');
    for i = 1:numel(x)
        fprintf('  %-6d %-16.8f %-16.8f %-12.3e\n', i-1, ...
            U(i,n+1), Uex(i,n+1), abs(U(i,n+1)-Uex(i,n+1)));
    end
end
end

function plotProfiles(U, Uex, x, t, times, titleText, yLabel)
%PLOTPROFILES  Numerical and analytical profiles at selected times.
%   The exact solution is a continuous line and the numerical one is
%   thinned markers, so agreement is visible rather than hidden.
stride = max(1, round(numel(x)/25));
figure('Name', titleText, 'NumberTitle', 'off', 'Color', 'w');
for k = 1:numel(times)
    [~, n] = min(abs(t - times(k)));
    subplot(2, 2, k);
    plot(x, Uex(:,n), '-', 'LineWidth', 2, 'Color', [.15 .15 .15]);
    hold on;
    plot(x(1:stride:end), U(1:stride:end,n), 'o', 'MarkerSize', 5, ...
        'LineWidth', 1.2, 'Color', [.85 .33 .10]);
    hold off; grid on;
    xlabel('x  [m]'); ylabel(yLabel);
    title(sprintf('t = %.4e s', t(n)));
    legend('analytical', 'numerical', 'Location', 'best');
end
end

function plotErrorHistory(U, Uex, x, t, tagText)
%PLOTERRORHISTORY  Discrete L2 and L-infinity error versus time.
[~, eL2, eInf] = calculateErrors(U, Uex, x);
figure('Name', [tagText ': error vs time'], 'NumberTitle', 'off', ...
       'Color', 'w');
plot(t, eL2, '-', 'LineWidth', 1.8, 'Color', [0 .45 .74]); hold on;
plot(t, eInf, '--', 'LineWidth', 1.8, 'Color', [.85 .33 .10]);
hold off; grid on;
xlabel('t  [s]'); ylabel('error norm');
title([tagText ': growth of the discretisation error']);
legend('discrete L_2', 'L_\infty', 'Location', 'northwest');
end

function plotFieldHistory(tE, Erows, tH, Hrows)
%PLOTFIELDHISTORY  E at the cavity centre and H at the PEC wall, each
%   passed as [analytical; numerical].  The quarter-period phase shift
%   between them is the signature of a standing wave: energy sloshes
%   between the electric and magnetic reservoirs twice per cycle.
lbl = {'E_y(d/2,t)  [V/m]', 'H_z(\Deltax/2,t)  [mA/m]'};
ttl = {'Electric field at the cavity centre', ...
       'Magnetic field beside the PEC wall (quarter-period lag)'};
figure('Name', 'Example 3: E and H histories', 'NumberTitle', ...
       'off', 'Color', 'w');
for k = 1:2
    if k == 1, tv = tE; rows = Erows; else, tv = tH; rows = Hrows; end
    subplot(2,1,k);
    plot(tv, rows(1,:), '-', 'LineWidth', 2, 'Color', [.15 .15 .15]);
    hold on;
    plot(tv, rows(2,:), '--', 'LineWidth', 1.4, 'Color', [.85 .33 .10]);
    hold off; grid on;
    xlabel('t  [ps]'); ylabel(lbl{k}); title(ttl{k});
    legend('analytical', 'numerical', 'Location', 'northeast');
end
end