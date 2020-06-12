rng(0);

% parameters
p = asteroid_params();
ub = 10;
dt = 0.1;

n_steps = 100;
u_old = randn(3, n_steps);

Qf = eye(6);
Q = eye(6);
R = eye(3);

s0 = [0; 0; 10000; 0; 0; 0];
sf = [10; 0; 10000; 0; 0; 0];

[s, u] = scp_solution(Q, R, Qf, ub, sf, s0, n_steps, dt, p);