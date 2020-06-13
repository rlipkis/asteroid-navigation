function [A, B] = linear_dynamics(s_ref, p)
    % unpack
    r = s_ref(1:3);
    v = s_ref(4:6);
    
    % interior correction
    [corr1, corr2] = correction(r, p);
    
    % dynamics
    A = zeros(6);
    A(1:3,4:6) = eye(3);
    A = A + corr1*A_build(r - p.r1, p.mu1);
    A = A + corr2*A_build(r - p.r2, p.mu2);
    
    B = zeros(6,3);
    B(4:6,:) = eye(3);
end

function A = A_build(dr, mu)
    % graviational component to dynamics matrix
    A = zeros(6);
    A(4,1) = mu*(2*dr(1)^2 - dr(2)^2 - dr(3)^2)/norm(dr)^5;
    A(4,2) = mu*(3*dr(1)*dr(2))/norm(dr)^5;
    A(4,3) = mu*(3*dr(1)*dr(3))/norm(dr)^5;  
    A(5,1) = A(4,1);
    A(5,2) = mu*(-dr(1)^2 + 2*dr(2)^2 - dr(3)^2)/norm(dr)^5;
    A(5,3) = mu*(3*dr(2)*dr(3))/norm(dr)^5;  
    A(6,1) = A(4,3);
    A(6,2) = A(5,3);
    A(6,3) = mu*(-dr(1)^2 - dr(2)^2 + 2*dr(3)^2)/norm(dr)^5;
end