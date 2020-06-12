rng(0);

% parameters
p = asteroid_params();
ub = 0.001;
dt = 1;

n_steps = 200;
u_old = randn(3, n_steps);

Qf = eye(6);
Q = eye(6);
R = eye(3);

s0 = [0; 0; 5000; 0; 0; 0];
sf = [10; 0; 5000; 0; 0; 0];
    
n = 6;
m = 3; 

s_ref = zeros(n*n_steps, 1);
s_ref(1:n) = s0;
u_ref = 1e-3*randn(m*n_steps, 1);

% initial forward pass
for i=1:(n_steps-1)
    si = s_ref((i-1)*n+1:i*n);
    ui = u_ref((i-1)*m+1:i*m);
    s_ref(i*n+1:(i+1)*n) = si + dynamics(si, ui, p)*dt;
end

% % debugging
% r = zeros(n_steps, 3);
% for i = 1:n_steps
%      r(i,:) = s_ref((i-1)*n+1:(i-1)*n+3);
% end
% plot3(r(:,1), r(:,2), r(:,3))
    
[s, u] = scp(s_ref, u_ref, ub, Q, R, Qf, sf, s0, n_steps, dt, p);
ctr = 0;

while (norm(s - s_ref, 'inf') + norm(u - u_ref, 'inf') > 0.5)
    ctr = ctr + 1;
    fprintf('Iteration %i: ', ctr);
    s_ref = s;
    u_ref = u;
    [s, u] = scp(s_ref, u_ref, ub, Q, R, Qf, sf, s0, n_steps, dt, p);
end