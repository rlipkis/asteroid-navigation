clear all

% parameters
p = asteroid_params();
ub = 0.01; % was 1e-3
dt = 10;

n_steps = 500;
u_old = randn(3, n_steps);

Q = zeros(6);
Qf = diag([1 1 1 10 10 10]);
R = 1e3*eye(3);

s0 = [0; 0; 5000; 0; 0; 0];
sf = [-2666; 0; 1000; 0; 0; 0];
    
n = 6;
m = 3; 

s_ref = zeros(n*n_steps, 1);
s_ref(1:n) = s0;
u_ref = zeros(m*n_steps, 1);

% initial forward pass
for i=1:(n_steps-1)
    si = s_ref((i-1)*n+1:i*n);
    ui = u_ref((i-1)*m+1:i*m);
    s_ref(i*n+1:(i+1)*n) = si + dynamics(si, ui, p)*dt;
end

iter = 1;
tol = 1e-1;
err = Inf;
while err > tol   
    fprintf('Iteration %i: ', iter);
    [s, u] = scp(s_ref, u_ref, ub, Q, R, Qf, sf, s0, n_steps, dt, p);
    err = norm(s - s_ref, 'inf') + norm(u - u_ref, 'inf');
    s_ref = s;
    u_ref = u;
    iter = iter + 1;
end
display_sol(s, u)