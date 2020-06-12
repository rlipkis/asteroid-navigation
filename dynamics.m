function f = dynamics(s_ref, u_ref, p)
    % unpack
    r = s_ref(1:3);
    v = s_ref(4:6);
    
    % dynamics
    dr1 = r - p.r1;
    dr2 = r - p.r2;
    
    f = zeros(6,1);
    f(1:3) = v;
    f(4:6) = -p.mu1*dr1/norm(dr1)^3 - p.mu2*dr2/norm(dr2)^3 + u_ref;    
end