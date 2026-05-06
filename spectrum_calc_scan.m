%% Calculate spectrum

% Choose energy and wavenumber bounds to calculate spectrum over:
kbound = [-15, 15];
Ebound = [0, 1.2];

% Discretization parameters. Nk, Ne; initial size of (k, E) grid to
% calculate over. L; number of adaptive grid steps. Number of grid points
% in each direction is doubled each step so that the final grid size is
% Nk*2^(L-1) by Ne*2^(L-1). N is number of discretization points for the
% ODE solver to use at each point.
Nk = 600;
Ne = 1200;
N = 256;
L = 2;
x = linspace(-1, 1, N+1);

% Define one or more functions for \Omega(y). For discontinuous profiles
% define two functions, one for each continuous interval:

%p1 = spline([x(1), x(N/2+1), x(end)], [0, -1, 0, 1, 0]);
% p2 = spline([x(N/2+1), x(3*N/4+1), x(N+1)], [0, 0.5, 0, 1, 0]);
%m1 = @(s) ppval(p1, s);
% m2 = @(s) ppval(p2, s);
m1 = @(s) 2*s + 1;
m2 = @(s) 2*s - 1;


B = @cold_plasma_tm_sc;
tic;

% Choose initial error e, final error ef, and specify the locations of
% discontinuities in terms of the coordinate of x. Then calculate spectrum:
e = 0.5;
ef = 10^(-5);
discont = [N/2+1]; % if none set to []
[Ae] = grid_adapt(kbound, Ebound, Nk, Ne, e, ef, L, B, {m1, m2}, x, discont);
% [Ae] = residue_map(kbound, Ebound, Nk, Ne, B, {m1, m2}, x, discont);
% Functions "residue_map.m" and "uniform_grid.m" may also be used but do
% not perform as well.
t = toc/60

%% Plot Spectrum:

%figure, imshow(Ae, [0, .1]);
k = linspace(kbound(1), kbound(2), Nk*2^L+1);
E = linspace(Ebound(1), Ebound(2), Ne*2^L+1);
plot_grid(Ae, ones(size(Ae)), ones(size(Ae)), ones(size(Ae)), k, E, 10^(-5));
axis([kbound(1), kbound(2), Ebound(1), Ebound(2)])