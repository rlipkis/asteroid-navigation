clear all

% parameters
p = asteroid_params();
ub = 0.01; % was 1e-3
dt = 10;

n_steps = 500; %200; % 500;

Q = zeros(6);
Qf = diag([1 1 1 10 10 10]);
R = 1e3*eye(3);

% r0 = p.r2 + (p.R2 + 5)*[4; 0; 1]/sqrt(17);
% rf = p.r1 + (p.R1 + 5)*[0; -1; 0];
% s0 = [r0; 0; 0; 0];
% sf = [rf; 0; 0; 0]; % [0; 0; 5000; 0; 0; 0];
r0 = p.r2 + (p.R2 + 5)*[4; 0; 1]/sqrt(17);
s0 = [r0; 0; 0; 0];
rf = [0; 0; 5000; 0; 0; 0];

n = 6;
m = 3; 

s_ref = zeros(n*n_steps, 1);
s_ref(1:n) = s0;
u_ref = zeros(m*n_steps, 1);
% guess for ascent
u_ref(3:m:end) = 4e-4; % was 1e-3
u_ref(2:m:end) = 0; %-1e-3; % 0
u_ref(1:m:end) = 0; %-4e-4; % 0

% initial forward pass
for i=1:(n_steps-1)
    si = s_ref((i-1)*n+1:i*n);
    ui = u_ref((i-1)*m+1:i*m);
    s_ref(i*n+1:(i+1)*n) = si + dynamics(si, ui, p)*dt;
end
% display_sol(s_ref, u_ref, dt)
% return

% iteration
iter = 1;
tol = 1e-1;
err = Inf;
conv = [];
while err > tol   
    fprintf('Iteration %i: ', iter);
    [s, u] = scp(s_ref, u_ref, ub, Q, R, Qf, sf, s0, n_steps, dt, p);
    err = norm(s - s_ref, 'inf') + norm(u - u_ref, 'inf');
    if isnan(err)
        s = s_ref;
        u = u_ref;
        break
    end
    s_ref = s;
    u_ref = u;
    conv = [conv, err];
    iter = iter + 1;
end
display_sol(s, u, dt)
plot(conv)