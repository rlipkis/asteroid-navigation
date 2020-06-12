x0 = [1;1;1];
xc = [0.5;0.5;0.5];

cvx_begin
    variable x(3)
    cost = quad_form(x - x0, eye(3))    
    minimize(cost)  
cvx_end