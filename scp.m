function [s, u] = scp(s_ref, u_ref, u_ub, Q, R, Qf, sf, s0, n_steps, dt, p)
    %% indexing 
    z_ref = [s_ref; u_ref];
    n = 6;
    m = 3;
    n_var = (n + m)*n_steps;
    u_shift = n*n_steps;
    s_start = @(i) (i-1)*n + 1;
    s_end = @(i) i*n;
    u_start = @(i) u_shift + (i-1)*m + 1;
    u_end = @(i) u_shift + i*m;
    
    %% bounds
    ub = 10000*ones(n_var, 1);
    lb = -ub;
    ub(u_start(1):end) = u_ub;
    lb(u_start(1):end) = -u_ub;

    %% target state
    z0 = zeros(n_var, 1);
    for i=1:n_steps
        z0(s_start(i):s_end(i)) = sf;
    end

    %% build cost matrix
    M = zeros(n_var);
    for i=1:n_steps
        if i < n_steps
            M(s_start(i):s_end(i), s_start(i):s_end(i)) = Q;
        else
            M(s_start(i):s_end(i), s_start(i):s_end(i)) = Qf;
        end
        M(u_start(i): u_end(i), u_start(i): u_end(i)) = R;
    end

    %% linear dynamics constraints
    n_constr = n*n_steps;

    % constraints will have the form C*z == d
    C = zeros(n_constr, n_var);
    d = zeros(n_constr, 1);

    for i = 1:n_steps-1
        sk = z_ref(s_start(i):s_end(i));
        uk = z_ref(u_start(i):u_end(i));
        [A, B] = linear_dynamics(sk, p);
        fk = dynamics(sk, uk, p);

        C(s_start(i+1):s_end(i+1), s_start(i+1):s_end(i+1)) = eye(n);
        C(s_start(i+1):s_end(i+1), s_start(i):s_end(i)) = -eye(n) - A*dt; 
        C(s_start(i+1):s_end(i+1), u_start(i):u_end(i)) = -B*dt;
        d(s_start(i+1):s_end(i+1)) = (fk - A*sk - B*uk)*dt;
    end
    
    % initial condition constraint
    C(s_start(1):s_end(1), s_start(1):s_end(1)) = eye(n);
    d(s_start(1):s_end(1)) = s0;
    
    % trust radius
    rho = 10.0;
    
    %% optimization
    cvx_begin quiet

        variable z(n_var);
        cost = quad_form(z - z0, M);

        minimize(cost)
        subject to
        % dynamics constraints
        C*z == d;
        % bounds
        lb <= z <= ub;
        % trust region
        abs(z - z_ref) <= rho;
    cvx_end

    var = z;
    s = var(1:n*n_steps);
    u = var(u_shift+1:end);
    
    % printing
    fprintf('%s / ', cvx_status)
    fprintf('current cost: %f \n', (var - z0)'*M*(var - z0))
end
